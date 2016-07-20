#!/usr/bin/python

import os
from astropy.time import Time
import psycopg2
import sys
import socket
from optparse import OptionParser
import subprocess
import datetime
import time

#********************************
#Script that downloads GPU box files, runs cotter, updates the mwa_qc database with the uvfits
#locations, and deletes the GPU box files. It performs a check for existing GPU box files
#and will not download them if they already exist. Option for downloading uvfits files and 
#bypassing cotter. Files are downloaded to whichever nodes have sufficient free space. Download
#and cotter processes are parallelized with grid engine for efficiency.

#Script written by Nichole Barry and Ruby Byrne, July 2016.

def main():

	#Parse the command line inputs. 
	#Set version and subversion to a number to indicate which cotter settings to use (defined in 
	#run_cotter in the cotter_args dict). Required unless uvfits_download_check is set, which 
	#forces uvfits download only and bypasses cotter.
	parser=OptionParser()
	parser.add_option("-v", "--version", dest="version", \
		help="Version for cotter arguments defined in run_cotter function")
	parser.add_option("-s", "--subversion", dest="subversion", \
		help="Subversion for cotter arguments defined in run_cotter function")
	parser.add_option("-o", "--obsfile_name", dest="obsfile_name", \
		help="Path to a file with the observation IDs for processing")
	parser.add_option("-u", "--uvfits_download_check", dest="uvfits_download_check", \
		help="Download only the uvfits from the ngas server")
	parser.add_option("-c", "--db_comment", dest="db_comment", \
		help="Optional comment to be placed in the uvfits database.",default='')
	(options, args)=parser.parse_args()
	version=options.version
	subversion=options.subversion
	obsfile_name=options.obsfile_name
	uvfits_download_check=options.uvfits_download_check
	db_comment=options.db_comment

	if (version is None) and (subversion is None) and (uvfits_download_check is None):
		print "ERROR: version, subversion, and uvfits_download_check were not set."
		print "To run cotter, set a version and subversion (defined in script)."
		print "To download uvfits files only and bypass cotter, set -u 1 on the command line."
		sys.exit(1)

	obs_per_chunk = 10 #number of obsids to run in parallel

	#get obsids to download:
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()
	nonredundant_obsids = list(set(obsids))
	if len(obsids) != len(nonredundant_obsids):
		print "WARNING: Obs list contains redundant entries."
	obsids = nonredundant_obsids	  
    
	#separate obsids into chunks:
	obs_chunks = []
	obs_start = 0
	while obs_start <= len(obsids)-1:
		obs_chunks.append(obsids[obs_start:obs_start+obs_per_chunk])
		obs_start += obs_per_chunk
	chunk_submitted = [False for i in range(len(obs_chunks))]

	#find which nodes have enough space for downloads:
	all_nodes = ["eor-02", "eor-03", "eor-04", "eor-05", "eor-06", "eor-07", "eor-08", "eor-10", "eor-11", "eor-12", "eor-13", "eor-14"]
	all_nodes = ["/nfs/" + nodename + "/r1/" for nodename in nodes]
	free_nodes = filespace(all_nodes)
	if len(free_nodes) == 0:
		print "ERROR: No file space found."
		sys.exit(1)

	children = []
	nodes = [1,2,3]
	for process, node in enumerate(nodes):
		pid = os.fork()
		if pid:
			children.append(pid)
		else:
			node
			os._exit(0)
	print children
	for child in children:
		os.waitpid(child, 0)


	for i in range(len(obs_chunks)):
		obs_chunk = obs_chunks[i]
		t = Time([int(obs) for obs in obs_chunk], format="gps", scale="utc")
		chunk_jds = t.jd
		chunk_jds = [int(date) for date in chunk_jds]
		save_directories = ["EoRuvfits/jd" + str(jd) + "v"+ str(version) + "_" + str(subversion) + "/" for jd in chunk_jds]

		#check to see if GPU box files already exist:
		download = [True for i in range(len(obs_chunk))] #indicates which obsids need to be downloaded
		save_paths = []
		for i in range(len(obs_chunk)):
			gpu_loc_path = find_gpubox(obs_chunk[i], save_directories[i], all_nodes)
			if gpu_loc_path == False:
				save_paths.append(node + save_directories)
			else:
				save_paths.append(gpu_loc_path)
				download[i] = False

		#Find the path to python using a child process
		stdoutpointer = subprocess.Popen(["which","python"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout_data, stderr_data = stdoutpointer.communicate()

		#Check to see if a path was found for python
		if stderr_data:
			print 'ERROR: The command "which python" did not return the path of the python installation you want to use.'
			print 'Please add the path to python.'
			sys.exit(1)

		python_path = stdout_data

		#Check to make sure the log_files directory exists on the node
		if not os.path.exists(node + 'EoRuvfits/log_files/'):
			os.makedirs(node + 'EoRuvfits/log_files/')

		#Check to make sure the obsid directory exists on the node for each obsid
		#Otherwise, script runs ahead of Grid Engine, and requires directory before GE can make it.
		for i in range(len(obs_chunk)):
			if not os.path.exists(save_paths[i] + obs_chunk[i]):
				os.makedirs(save_paths[i] + obs_chunk[i])

		#initialize
		task_jobid = False

		#Download the files (a uvfits or gpuboxes depending on uvfits_download_check)
		if any(download) or uvfits_download_check:
			(task_jobid, download_script_path) = download_files(save_paths, obs_chunk, uvfits_download_check, python_path, node, download)
			download_script_paths.append(download_script_path)


		#If metafits does not exist in the same location as the gpubox files, set up logic to create it
		metafits_logic = []
		for i in range(len(obs_chunk)):
			if not os.path.isfile(save_paths[i] + obs_chunk[i] +'/' + obs_chunk[i] + '.metafits'):
				metafits_logic.append(True)
			else:
				metafits_logic.append(False)
				print "Using metafits file found for obsid " + obs_chunk[i] + " located in " + save_paths[i]

		if any(metafits_logic):
			#Make a metafits file for the obsids, will bypass if all the metafits exists.    
			(task_jobid, metafits_script_path) = make_metafits(obs_chunk, save_paths,task_jobid,python_path,node,metafits_logic)
			metafits_script_paths.append(metafits_script_path)
    
		#Run cotter if gpubox files were downloaded
		if not uvfits_download_check:
			(task_jobid, cotter_version, cotter_script_path) = run_cotter(version,subversion,save_paths,obs_chunk,task_jobid,node)
			cotter_script_paths.append(cotter_script_path)
    
		chunk_submitted[chunk_index] = True

		#Grab the last Grid Engine jobid to watch while the program sleeps
		final_task_jobids.append(task_jobid)
   
		#sleep while periodically checking to see if a job finishes in grid engine:
		completed_node = []
		sleep_time = 20 #check after this number of seconds
		job_finish_array = [False for obsid in obs_chunk]

		#Check each of tasks in the task array for the last submitted job
		for task_array_index in range(len(obs_chunk)):
			#Talk to Grid Engine about the last submitted job for one of the tasks
			qsub_command = 'qacct -j ' + str(task_jobid) + ' -t ' + str(task_array_index)
			stdoutpointer = subprocess.Popen(qsub_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout_data, stderr_data = stdoutpointer.communicate()

			#If the command that only works when the job is done does not throw an error, then the job finished
				if not stderr_data:
					exit_start = stdout_data.rfind('exit_status') + 11
					exit_end = stdout_data.rfind('ru_wallclock')
					exit_status = stdout_data[exit_start:exit_end].split()[0]
					job_finish_array[task_array_index] = exit_status

		#If all of the tasks are done, then break the sleeper loop
		if all(job_finish_array):
			break
		free_nodes = filespace(completed_node)

		#Remove the temporary bash script for download, metafits, and cotter in Grid Engine
		os.remove(download_script_path)
		os.remove(metafits_script_path)
		os.remove(cotter_script_path)



		#Fill the database with the new uvfits
		fill_database(obs_chunk,version,subversion,save_paths,cotter_version,db_comment,uvfits_download_check)

		#Check that all gpubox files were successfully downloaded


		#Delete the gpubox files
		delete_gpubox(obs_chunk,save_paths)



  




    while chunk_submitted.count(False) > 0:
        
        final_task_jobids = []
	download_script_paths = []
	metafits_script_paths = []
	cotter_script_paths = []

        for node in free_nodes[0:chunk_submitted.count(False)]: 
            chunk_index = chunk_submitted.index(False)
            obs_chunk = obs_chunks[chunk_index]
            t = Time([int(obs) for obs in obs_chunk], format="gps", scale="utc")
            chunk_jds = t.jd
            chunk_jds = [int(date) for date in chunk_jds]
            save_directories = ["EoRuvfits/jd" + str(jd) + "v"+ str(version) + "_" + str(subversion) + "/" for jd in chunk_jds]
            
            #check to see if GPU box files already exist:
            download = [True for i in range(len(obs_chunk))] #indicates which obsids need to be downloaded
            save_paths = []
            for i in range(len(obs_chunk)):
		gpu_loc_path = find_gpubox(obs_chunk[i], save_directories[i], nodes)
		if gpu_loc_path == False:
		    save_paths.append(node + save_directories)
		else:
		    save_paths.append(gpu_loc_path)
		    download[i] = False           
            
            #submit obs_chunk (list), save_paths (list), node (string), version (int), subversion(int), download (list) to grid engine for download and cotter
            #Find the path to python using a child process
            stdoutpointer = subprocess.Popen(["which","python"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout_data, stderr_data = stdoutpointer.communicate()

            #Check to see if a path was found for python
            if stderr_data:
                print 'ERROR: The command "which python" did not return the path of the python installation you want to use.'
                print 'Please add the path to python.'
                sys.exit(1)

            python_path = stdout_data
            
            #Check to make sure the log_files directory exists on the node
            if not os.path.exists(node + 'EoRuvfits/log_files/'):
                os.makedirs(node + 'EoRuvfits/log_files/')

            #Check to make sure the obsid directory exists on the node for each obsid
	    #Otherwise, script runs ahead of Grid Engine, and requires directory before GE can make it.
	    for i in range(len(obs_chunk)):
                if not os.path.exists(save_paths[i] + obs_chunk[i]):
                    os.makedirs(save_paths[i] + obs_chunk[i])

	    #initialize
	    task_jobid = False

            #Download the files (a uvfits or gpuboxes depending on uvfits_download_check)
	    if any(download) or uvfits_download_check:
                (task_jobid, download_script_path) = download_files(save_paths, obs_chunk, uvfits_download_check, python_path, node, download)
		download_script_paths.append(download_script_path)


    	    #If metafits does not exist in the same location as the gpubox files, set up logic to create it
    	    metafits_logic = []
    	    for i in range(len(obs_chunk)):
        	if not os.path.isfile(save_paths[i] + obs_chunk[i] +'/' + obs_chunk[i] + '.metafits'):
	            metafits_logic.append(True)
	        else:
	    	    metafits_logic.append(False)
	    	    print "Using metafits file found for obsid " + obs_chunk[i] + " located in " + save_paths[i]

    	    if any(metafits_logic):
            	#Make a metafits file for the obsids, will bypass if all the metafits exists.    
            	(task_jobid, metafits_script_path) = make_metafits(obs_chunk, save_paths,task_jobid,python_path,node,metafits_logic)
		metafits_script_paths.append(metafits_script_path)
            
            #Run cotter if gpubox files were downloaded
            if not uvfits_download_check:
                (task_jobid, cotter_version, cotter_script_path) = run_cotter(version,subversion,save_paths,obs_chunk,task_jobid,node)
		cotter_script_paths.append(cotter_script_path)
            
            chunk_submitted[chunk_index] = True

	    #Grab the last Grid Engine jobid to watch while the program sleeps
	    final_task_jobids.append(task_jobid)
	   
        #sleep while periodically checking to see if a job finishes in grid engine:
        completed_node = []
        sleep_time = 20 #check after this number of seconds
        job_finish_array = [False for obsid in obs_chunk]

        while len(completed_node) == 0:
	    #wait
            time.sleep(sleep_time)

	    #Check each of tasks in the task array for the last submitted job
	    for task_array_index in range(len(obs_chunk)):
                #Talk to Grid Engine about the last submitted job for one of the tasks
		qsub_command = 'qacct -j ' + str(task_jobid) + ' -t ' + str(task_array_index)
	        stdoutpointer = subprocess.Popen(qsub_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	        stdout_data, stderr_data = stdoutpointer.communicate()

		#If the command that only works when the job is done does not throw an error, then the job finished
	        if not stderr_data:
	            exit_start = stdout_data.rfind('exit_status') + 11
	            exit_end = stdout_data.rfind('ru_wallclock')
	            exit_status = stdout_data[exit_start:exit_end].split()[0]
		    job_finish_array[task_array_index] = exit_status

	    #If all of the tasks are done, then break the sleeper loop
	    if all(job_finish_array):
	        break
	    free_nodes = filespace(completed_node)

        #Remove the temporary bash script for download, metafits, and cotter in Grid Engine
        os.remove(download_script_path)
        os.remove(metafits_script_path)
        os.remove(cotter_script_path)



        #Fill the database with the new uvfits
        fill_database(obs_chunk,version,subversion,save_paths,cotter_version,db_comment,uvfits_download_check)

	#Check that all gpubox files were successfully downloaded
	
        
        #Delete the gpubox files
        delete_gpubox(obs_chunk,save_paths)
#********************************

#********************************
#Module that searches for saved GPU box files
def find_gpubox(obsid, save_directory, all_nodes):
	for gpu_loc_node in all_nodes:
		gpu_loc_path = gpu_loc_node + save_directory
	if os.path.isdir(gpu_loc_path + obsid): #checks to see if the directory exists
		directory_contents = os.listdir(gpu_loc_path + obsid)
		gpubox00 = 0
		gpubox01 = 0
		flags = 0
		metafits = 0
		for file in directory_contents: #counts how many of each type of file exists
			if file.endswith("_00.fits"):
				gpubox00 += 1
			if file.endswith("_01.fits"):
				gpubox01 += 1
			if file.endswith("_flags.zip"):
				flags += 1
			if file.endswith("_metafits_ppds.fits"):
				metafits += 1
		if gpubox00 >= 24 and gpubox01 >= 24 and flags >= 1 and metafits >= 1:
			print "Using saved GPU box files for obsid " + obsid + " located in " + gpu_loc_path
			if gpubox00 != 24 or gpubox01 != 24 or flags != 1 or metafits != 1:
		     		print "WARNING: Directory contains extra GPU box files."
			return gpu_loc_path   
	return False
#********************************   

#********************************
#Module that takes a list of file paths and returns those that have more than 
#a specified amount of free disk space; prints a warning if the free disk space
#is less than another specified amount
def filespace(nodes):
	enough_space = 2 #disks with less than this amount in TB of free space will not be used
	space_warning = 4 #disks with less than this amount in TB of free space will return a warning
	free_nodes = []
	for node in nodes:
		if os.path.isdir(node):
			stat = os.statvfs(node)
			free_space = (stat.f_bavail * stat.f_frsize)/1024.**4
			if free_space > enough_space:
				free_nodes.append(node)
				if free_space < space_warning:
					print "WARNING: Limited disk space in " + node
			else:
				print "WARNING: No disk space in " + node
		else:
			print "WARNING: Disk " + node + " not found."
	return free_nodes
#********************************

#********************************
def download_files(save_paths, obs_chunk, uvfits_download_check, python_path, node, download):

	#Check to see that the MWA_Tools is in the path so obsdownload.py can be found
	#Find the path of MWA_Tools by looking in the system path variable 
	mwa_tools_path=""
	for parsed_path in os.environ['PATH'].split(':'):
		if "MWA_Tools" in parsed_path:
			mwa_tools_path = parsed_path

	#If the MWA_Tools path doesn't exist, throw an error.
	if not mwa_tools_path:
		print 'ERROR: MWA_Tools is not in the path, obsdownload.py not found!'
		print 'Please add the path to MWA_Tools to your system path.'
		sys.exit(1)
			
	#Setup the path to the download script and log files
	obsdownload_path = mwa_tools_path[0:mwa_tools_path.find("MWA_Tools")+9] + '/scripts/obsdownload.py'
	log_path = (save_paths[0])[0:(save_paths[0]).rfind("jd")] + "log_files/"
			
	#Pick out the obsids that need to be downloaded given the download logic created in main.
	obs_chunk_download=[]
	save_paths_download=[]
	for i in range(len(download)):
		if download[i]:
			obs_chunk_download.append(obs_chunk[i])
			save_paths_download.append(save_paths[i])

	#Setup the u tag argument, with or without the uvfits only option	
	if not uvfits_download_check:	 
		u_arg = ''
	else:
		u_arg = ' -u 1'

	#Write a bash script so that Grid Engine can run a task array for the downloads.
	download_script_path = save_paths_download[0] + 'download_commands_file_chunk'+obs_chunk_download[0]+'.sh'
	download_commands_file = open(download_script_path, 'w')
	#Write the contents of the file and the necessary arguments
	download_commands_file.write('#!/bin/bash\n\n' + \
		'#$ -S /bin/bash\n\n' + \
		'save_paths_download=(0 '+" ".join(save_paths_download) + ')\n' + \
		'obs_chunk_download=(0 ' +  " ".join(obs_chunk_download) + ')\n' + \
		'ls ${save_paths_download[$SGE_TASK_ID]} > /dev/null \n' + \
		python_path.strip('\n') + ' ' + obsdownload_path + ' -d ${save_paths_download[$SGE_TASK_ID]} -o ${obs_chunk_download[$SGE_TASK_ID]} ' + u_arg)
	#Close the file
	download_commands_file.close() 

	#Make the file executable 
	os.chmod(download_script_path, 0775)

	#Setup the the bulk of the qsub command
	qsub_command = "qsub -l h_vmem=1G,h_stack=512,h_rt=08:00:00,h=" + node.split('/')[2] + ".mit.edu -pe chost 1 -e " + \
		log_path + " -o " + log_path + " -N obs_download -t 1:" + str(len(obs_chunk_download))

	#Run the qsub command
	stdoutpointer = subprocess.Popen((qsub_command + ' ' + download_script_path).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()
		
	#If Grid Engine had an error, report and exit
	if stderr_data:
		print "WARNING: Grid Engine threw an error trying to run obsdownload.py."
		print stderr_data
		sys.exit(1)

	#Get the jobid of the download task array to hold future jobs till it finishes.	  
	task_jobid = (stdout_data.split())[2].rsplit('.')[0]
	
	return (task_jobid, download_script_path)
#********************************           

#********************************
#Module for running cotter given input version and subversion, and for finding 
#the version of cotter itself. Will make a metafits if it does not exist in the
#save_path directory
def run_cotter(version,subversion,save_paths,obs_chunk,task_jobid,node):

	#Perform check to make sure essential information is known for the module
	if not obs_chunk:
		print "ERROR: obs_chunk not defined in run_cotter."
		sys.exit(1)
	if not version:
		print "ERROR: version not defined in run_cotter"
		sys.exit(1)
	if not subversion:
		print "ERROR: subversion not defined in run_cotter"
		sys.exit(1)
	if not save_paths:
		print "ERROR: save_paths not defined in run_cotter"
		sys.exit(1)
	if not node:
		print "ERROR: node not defined in run_cotter"
		sys.exit(1)
	
	#Warn the user if uvfits were automatically downloaded for at least one observation in the 
	#chunk. Will delete this automatic uvfits and rerun cotter with the specifications.
	uvfits_logic = []
	for i in range(len(save_paths)):
		if os.path.isfile(save_paths[i] + obs_chunk[i] + '/' + obs_chunk[i] + '.uvfits'):
			uvfits_logic.append(True)
			os.remove(save_paths[i] + obs_chunk[i] + '/' + obs_chunk[i] + '.uvfits')

	if any(uvfits_logic):
		print "WARNING: At least one uvfits was automatically downloaded, which will be deleted and rerun with your cotter specifications."
		print "Please set -u to download uvfits only and bypass cotter."

	#A dictionary of all version and subversion cotter arguments. FEEL FREE TO ADD MORE VERSIONS,
	#but please add a comment below.
	#3,3 was used to test compressed fits
	#3,4 was a rerun of 3,1 with a newer version of cotter before that version was recorded
	#4,0 went back to old settings for an industrial run
	#4,1 was the same as 4,0 but for running on compressed gpubox files
	#5,0 was a test to phase all obs to zenith (phasing needs to be added per obs currently)
	cotter_args = { \
		"0,0": "-timeavg 4 -freqavg 2 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"1,0": "-timeavg 4 -freqavg 2 -flagedges 2", \
		"2,0": "-timeavg 4 -freqavg 2 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"2,1": "-timeavg 4 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"2,2": "-timeavg 4 -freqavg 2 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"2,3": "-timeavg 1 -freqavg 2 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"3,0": "-timeavg 4 -freqavg 2 -flagedges 2 -usepcentre -initflag 2 -noflagautos", \
		"3,1": "-timeavg 4 -freqavg 1 -edgewidth 80 -usepcentre -initflag 0 -noflagautos", \
		"3,2": "-timeavg 1 -freqavg 2 -edgewidth 80 -usepcentre -initflag 0 -noflagautos", \
		"3,3": "-timeavg 1 -freqavg 2 -edgewidth 80 -usepcentre -initflag 0 -noflagautos", \
		"3,4": "-timeavg 4 -freqavg 1 -edgewidth 80 -usepcentre -initflag 0 -noflagautos", \
		"4,0": "-timeres 2 -freqres 80 -edgewidth 80 -usepcentre -initflag 2 -noflagautos", \
		"4,1": "-timeres 2 -freqres 80 -edgewidth 80 -usepcentre -initflag 2 -noflagautos", \
		"5,0": "-timeres 2 -freqres 80 -edgewidth 80 -initflag 2 -noflagautos" \
		}

	#Check that the version and subversion supplied exist in the argument dictionary
	if not cotter_args.get(version + ',' + subversion):
		print 'ERROR: Cotter version ' + version + ', subversion ' + subversion + ' does not exist in the python dictionary in module run_cotter.'
		print 'Please choose another version and subversion or update the run_cotter module.'
		sys.exit(1)

	#Find the path to cotter using a child process
	stdoutpointer = subprocess.Popen(["which","cotter"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()

	#Check to see if a path was found for cotter
	if stderr_data:
		print 'ERROR: The command "which cotter" did not return the path of the cotter installation you want to use.'
		print 'Please add the path to cotter.'
		sys.exit(1)

	#Path setups
	cotter_path = stdout_data.strip('\n')
	metafits_path = [save_paths[i] + obs_chunk[i] + '/' + obs_chunk[i] + '.metafits' for i in range(len(obs_chunk))]
	uvfits_path = [save_paths[i] + obs_chunk[i] + '/' + obs_chunk[i] + '.uvfits' for i in range(len(obs_chunk))]
	gpubox_path = [save_paths[i] + obs_chunk[i] + '/' + obs_chunk[i] + '*gpubox*.fits' for i in range(len(obs_chunk))]

	#Find out the version of the found cotter using a child process
	stdoutpointer = subprocess.Popen((cotter_path + " --version").split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()

	#String manipulation to extract the version and version date
	cotter_version_string = stdout_data.splitlines()[0]
	cotter_version = cotter_version_string[cotter_version_string.rfind("version"):cotter_version_string.rfind(".")]
  
	#Setup the path to the log files
	log_path = (save_paths[0])[0:(save_paths[0]).rfind("jd")] + "log_files/"

	#Write a bash script so that Grid Engine can run a task array for the downloads.
	cotter_script_path = save_paths[0] + 'cotter_commands_file_chunk'+obs_chunk[0]+'.sh'
	cotter_commands_file = open(cotter_script_path, 'w')
	#Write the contents of the file and the necessary arguments
	cotter_commands_file.write('#!/bin/bash\n\n' + \
		'#$ -S /bin/bash\n\n' + \
		'uvfits_path=(0 '+" ".join(uvfits_path) + ')\n' + \
		'gpubox_path=(0 ' +  " ".join(gpubox_path) + ')\n' +  \
		'metafits_path=(0 ' + " ".join(metafits_path) + ')\n' + \
		'ls ${save_paths[$SGE_TASK_ID]} > /dev/null\n' + \
		cotter_path + ' ' + cotter_args[str(version)+','+str(subversion)] + ' -m ${metafits_path[$SGE_TASK_ID]} -o ${uvfits_path[$SGE_TASK_ID]} ${gpubox_path[$SGE_TASK_ID]}')
	#Close the file
	cotter_commands_file.close() 

	#Make the file executable 
	os.chmod(cotter_script_path, 0775)

	#Setup the bulk of the Grid Engine command, depending on if there is a task to wait on
	if task_jobid:
		qsub_command = "qsub -V -b y -hold_jid " + task_jobid + " -l h_vmem=5G,h_stack=512,h_rt=08:00:00,h=" + node.split('/')[2] \
			 + ".mit.edu -pe chost 1 -e " + log_path + " -o " + log_path +" -N cotter -t 1:" + str(len(obs_chunk)) + " " + cotter_script_path 
	else:
		qsub_command = "qsub -V -b y -l h_vmem=1G,h_stack=512,h_rt=08:00:00,h=" + node.split('/')[2] \
			 + ".mit.edu -pe chost 5 -e " + log_path + " -o " + log_path +" -N cotter -t 1:" + str(len(obs_chunk)) + " " + cotter_script_path 

	#Run cotter with the correct arguments, the path to the metafits, the uvfits output path, and the gpubox file paths
	stdoutpointer = subprocess.Popen((qsub_command + ' ' + cotter_script_path).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()  

	#If there is data in the standard error output from cotter, let the user know
	if stderr_data:
		print 'WARNING: Grid Engine has thrown an error in trying to run cotter.'
		print stderr_data
  
	#Get the jobid of the cotter task array for checking later.	  
	task_jobid = (stdout_data.split())[2].rsplit('.')[0]

	#return the version of cotter used for fill_database and the qsub job id
	return (task_jobid, cotter_version, cotter_script_path)
#********************************

#********************************
def make_metafits(obs_chunk, save_paths, task_jobid, python_path, node, metafits_logic):

	#Elements of obsids to make metafits for
	obs_elements = [i for i, x in enumerate(metafits_logic) if x]

	#Find the path of MWA_Tools by looking in the system path variable 
	mwa_tools_path=""
	for parsed_path in os.environ['PATH'].split(':'):
		if "MWA_Tools" in parsed_path:
			mwa_tools_path = parsed_path

	#If the MWA_Tools path doesn't exist, throw an error.
	if not mwa_tools_path:
		print 'ERROR: MWA_Tools is not in the path, make_metafits.py not found!'
		print 'Please add the path to MWA_Tools to your system path.'
		sys.exit(1)

	#Setup the path to make_metafits.py and to the metafits file
	make_metafits_path = mwa_tools_path[0:mwa_tools_path.find("MWA_Tools")+9] + '/scripts/make_metafits.py'
	metafits_path = [save_paths[obs_elements[i]] + obs_chunk[obs_elements[i]] + '/' + obs_chunk[obs_elements[i]] + '.metafits' for i in range(len(obs_elements))]

	#Setup the log path for Grid Engine
	log_path = (save_paths[0])[0:(save_paths[0]).rfind("jd")] + "log_files/"

	#Write a bash script so that Grid Engine can run a task array for the metafits.
	metafits_script_path = save_paths[obs_elements[0]] + 'metafits_commands_file_chunk'+obs_chunk[obs_elements[0]]+'.sh'
	metafits_commands_file = open(metafits_script_path, 'w')
	#Write the contents of the file and the necessary arguments
	metafits_commands_file.write('#!/bin/bash\n\n' + \
		'#$ -S /bin/bash\n\n' + \
		'save_paths_metafits=(0 '+" ".join(metafits_path) + ')\n' + \
		'obs_chunk_metafits=(0 ' +  " ".join([obs_chunk[obs_element] for obs_element in obs_elements]) + ')\n' + \
		'ls ${save_paths_metafits[$SGE_TASK_ID]} > /dev/null\n' + \
		python_path.strip('\n') + ' ' + make_metafits_path + ' -o ${save_paths_metafits[$SGE_TASK_ID]} --gps ${obs_chunk_metafits[$SGE_TASK_ID]} ')
	#Close the file
	metafits_commands_file.close() 

	#Make the file executable 
	os.chmod(metafits_script_path, 0775)

	#Setup the bulk of the Grid Engine command, depending on if there is a task to wait on
	if task_jobid: 
		qsub_command = "qsub -hold_jid " + task_jobid + " -l h_vmem=1G,h_stack=512,h_rt=08:00:00,h=" + node.split('/')[2] + \
			".mit.edu -V -pe chost 1 -e " + log_path + " -o " + log_path +" -N make_metafits -t 1:" + str(len(obs_elements))
	else:
		qsub_command = "qsub -l h_vmem=1G,h_stack=512,h_rt=08:00:00,h=" + node.split('/')[2] + \
			".mit.edu -V -pe chost 1 -e " + log_path + " -o " + log_path +" -N make_metafits -t 1:" + str(len(obs_elements))

	#Run the metafits script to create a metafits for the obsid in the uvfits folder
	stdoutpointer = subprocess.Popen((qsub_command + ' ' + metafits_script_path).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()

	#If there is data in the standard error output from the metafits script, let the user know
	if stderr_data:
		print 'WARNING: Grid Engine threw an error trying to run make_metafits.py'
		print stderr_data

	#Get the jobid of the metafits task array to hold future jobs till it finishes.	  
	task_jobid = (stdout_data.split())[2].rsplit('.')[0]
		
	return (task_jobid, metafits_script_path)
#********************************

#********************************
#Module designed to fill the uvfits table on the mwa_qc database located on eor-00.
#Fills the uvfits table with the uvfits version and subversion for the obsid run, as well as
#the path to the uvfits, the cotter version that generated it, a timestamp for when it was
#run, an optional comment, and the bottom and top of the frequency range the uvfits file
#spans.
def fill_database(obs_chunk,version,subversion,save_paths,cotter_version,db_comment,uvfits_download_check):

	#Connect to the database on eor-00 using the mwa username.
	try:
		conn = psycopg2.connect(database='mwa_qc',user='mwa',password='BowTie',host='eor-00.mit.edu')
	except:
		print 'Could not connect to mwa database.'
		sys.exit(1)

	#Module only works on the MIT cluster.
	if not 'mit.edu' in socket.gethostname():
		print 'Sorry, this script is currently only supported on eor-xx.mit.edu machines.'
		sys.exit(1)  
		
	cur = conn.cursor()

	#Perform check to make sure essential information is known about the uvfits file
	if not obs_chunk:
		return "WARNING: obs_chunk not defined in fill_database. Database not updated"
	if not version and not uvfits_download_check:
		return "WARNING: version not defined in fill_database. Database not updated"
	if not subversion and not uvfits_download_check:
		return "WARNING: subversion not defined in fill_database. Database not updated"
	if not save_paths:
		return "WARNING: save_paths not defined in fill_database. Database not updated"
	if not cotter_version and not uvfits_download_check:
		return "WARNING: cotter_version not defined in fill_database. Database not updated"

	iteration = 0
	for obsid in obs_chunk:
		save_path = save_paths[iteration] + obsid + '/'
		iteration = iteration + 1

		#Check to make sure the uvfits and metafits specified exist
		if not os.path.isfile(save_path + obsid + '.uvfits'):
			print "ERROR: " + save_path + obsid + ".uvfits does not exist! Database not updated"
			sys.exit(1)
		if not os.path.isfile(save_path + obsid + '.metafits'):
			print "WARNING: " + save_path + obsid + ".metafits does not exist! Database not updated"
			return

		#Open up the metafits file that was made with the uvfits file (assumes they are in the same location)
		metafits_file = save_path + obsid + '.metafits'
		hdu_list_metafits = fits.open(metafits_file)
		header_metafits = hdu_list_metafits[0].header

		#Get the frequency center and bandwidth (assumes contiguous frequency channels)
		freq_cent = header_metafits['FREQCENT']
		bandwidth = header_metafits['BANDWDTH']
		top_freq_mhz = "{0:.2f}".format(freq_cent + bandwidth/2)
		bottom_freq_mhz = "{0:.2f}".format(freq_cent - bandwidth/2)

		#Close the metafits file
		hdu_list_metafits.close()

		#Get the time for the timestamp in UTC
		timestamp = str(datetime.datetime.utcnow())

		#Check to make sure that a similar uvfits file does not already exist, throw warning if it does.
		cur.execute("SELECT uvfits.obsid FROM uvfits WHERE (obsid,version,subversion,cotter_version,bottom_freq_mhz,top_freq_mhz)=(%s,%s,%s,%s,%s,%s);", \
			(obsid,version,subversion,cotter_version,bottom_freq_mhz,top_freq_mhz))
		if cur.fetchall():
			print "WARNING: A uvfits file of this obsid, version, subversion, cotter version, and frequency range already exists."

		#Create the database row, and fill it with the inputs. 
		cur.execute("INSERT INTO uvfits(obsid,version,subversion,path,cotter_version,timestamp,comment,bottom_freq_mhz,top_freq_mhz) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
			(obsid,version,subversion,save_path,cotter_version,timestamp,db_comment,bottom_freq_mhz,top_freq_mhz))

	#Commit all the cur.execute, and close the connection.
	conn.commit()
	cur.close()
	conn.close()

	#Print completion message
	print "Filled the uvfits table in the mwa_qc database with the current uvfits information." 
#********************************

#********************************
#Module for deleting gpubox files after the uvfits creation. Will check to see if a uvfits file
#exists before deletion
def delete_gpubox(obs_chunk,save_paths):

	iteration = 0
	for obsid in obs_chunk:
		save_path = save_paths[iteration]
		iteration = iteration + 1

		#Perform check to make sure essential information is known about the uvfits file
		if not obsid:
			return "WARNING: obsid not defined in delete_gpubox. Gpubox files not deleted"
		if not save_path:
			return "WARNING: save_path not defined in delete_gpubox. Gpubox files not deleted"
	
		#If the uvfits file does not exist with the gpubox files, do not delete the gpubox files
		if not os.path.isfile(save_path + str(obsid) + '.uvfits'):
			return "WARNING: uvfits file does not exist in the directory with the gpubox files. Gpubox files not deleted"

		#If the gpubox files do not exist, exit module
		if not os.path.isfile(save_path + str(obsid) + '*gpubox*.fits'):
			return "WARNING: there are not gpubox files to delete in " + save_path + " for obsid " + str(obsid) 

		#Remove gpubox files
		os.remove(save_path + str(obsid) + '*gpubox*.fits')

		print "Gpubox files in " + save_path + " for obsid " + str(obsid) + " have been deleted."
#********************************

if __name__ == '__main__':
	main()
