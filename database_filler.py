#!/usr/bin/python

import subprocess
import os
from astropy.time import Time
import sys
import psycopg2
import socket
from astropy.io import fits
import datetime
import time

def main():

	attempt_metafits = True #Set to false to fill database for obsids that already have metafits files (bypasses gridengine)

	version = 5
	subversion = 1
	obsfile_name = "/nfs/eor-00/h1/rbyrne/sidelobe_survey_obsIDs.txt"
	db_comment = "Diffuse survey EoR 2015"

	#Get obsids
	obsfile = open(obsfile_name, "r")
	all_obsids = [line.split( ) for line in obsfile.readlines()]
	all_obsids = [obs[0] for obs in all_obsids]
	obsfile.close()
	nonredundant_obsids = list(set(all_obsids))
	if len(all_obsids) != len(nonredundant_obsids):
		print "WARNING: Obs list contains redundant entries."
		all_obsids = nonredundant_obsids
	
	#Process the obsids in chunks	
	start_obs = 0	
	obs_per_chunk = 5
	while start_obs < len(all_obsids):
	
		obsids = all_obsids[start_obs:start_obs+obs_per_chunk]
		start_obs += obs_per_chunk
		
		print "Processing the following chunk: " 
		print obsids
		
		t = Time([int(obsid) for obsid in obsids], format="gps", scale="utc")
		jds = t.jd
		jds = [int(jd) for jd in jds]
		save_directories = ["EoRuvfits/jd" + str(jd) + "v"+ str(version) + "_" + str(subversion) + "/" for jd in jds]
	
		#Find file locations	
		all_nodes = ["eor-00", "eor-02", "eor-03", "eor-04", "eor-05", "eor-06", "eor-07", "eor-08", "eor-09", "eor-10", "eor-11", "eor-12", "eor-13", "eor-14"]
		node_paths = ["/nfs/" + nodename + "/" + disk + "1/" for nodename in all_nodes for disk in ["r","d","h"]]
		save_paths = []
		metafits_logic = []
		obsids_not_found = []
		for obsindex, obsid in enumerate(obsids):
			obsid_found = False
			for node in node_paths:
				if os.path.isdir(node + save_directories[obsindex] + obsid):
					if os.path.isfile(node + save_directories[obsindex] + obsid + "/" + obsid + ".uvfits"):
						save_paths.append(node + save_directories[obsindex])
						obsid_found = True						
						if os.path.isfile(node + save_directories[obsindex] + obsid + "/" + obsid + ".metafits"):
							metafits_logic.append(False) #Metafits file exists, will not be created
						else:
							metafits_logic.append(True) #Metafits file does not exist, will be created
						break
							
			if not obsid_found:
				print "uvfits file for obsid " + obsid + "not found."
				obsids_not_found.append(obsid)
			
		#Ignore obsids that did not have uvfits files
		for obsid in obsids_not_found:
			del obsids[obsids.index(obsid)]
		
		if attempt_metafits:
			if metafits_logic.count(True) > 0:
				for metafits_attempts in range(2): #Making metafits files sometimes failes; give it a second try if it fails the first time
					#Make metafits
					task_jobid = make_metafits(obsids, metafits_logic, save_paths)
		
					#Wait for metafits to finish being made
					stderr_data = True
					while stderr_data:
						time.sleep(30)
						#Talk to Grid Engine about the last submitted job for one of the tasks
						qsub_command = 'qacct -j ' + str(task_jobid)
						stdoutpointer = subprocess.Popen(qsub_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
						stdout_data, stderr_data = stdoutpointer.communicate()
						#The job is finished when stderr_data = False
						
					#Check to see if all metafits files were made
					for i, obsid in enumerate(obsids):
						if os.path.isfile(save_paths[i] + obsid + "/" + obsid + ".metafits"):
							metafits_logic[i] = False
						else:
							print "WARNING: Failed to make a metafits file for obsid " + obsid
					if metafits_logic.count(True) == 0:
						break
			else:
				print "Metafits files already exist, skipping make_metafits module."

		#Fill database
		fill_database(obsids, save_paths, version, subversion, db_comment)
		
		
def make_metafits(obsids, metafits_logic, save_paths):

	#Find the path to python using a child process
	stdoutpointer = subprocess.Popen(["which","python"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()
	if stderr_data:
		print 'ERROR: The command "which python" did not return the path of the python installation you want to use.'
		print 'Please add the path to python.'
		sys.exit(1)
	python_path = stdout_data

	#Elements of obsids to make metafits for
	obs_elements = [i for i, x in enumerate(metafits_logic) if x]
	print "Making metafits files for the following obsids:" 
	print [obsids[i] for i in obs_elements]

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
	metafits_path = [save_paths[obs_elements[i]] + obsids[obs_elements[i]] + '/' + obsids[obs_elements[i]] + '.metafits' for i in range(len(obs_elements))]

	#Setup the log path for Grid Engine
	log_path = (save_paths[0])[0:(save_paths[0]).rfind("jd")] + "log_files/"

	#Write a bash script so that Grid Engine can run a task array for the metafits.
	metafits_script_path = save_paths[obs_elements[0]] + 'metafits_commands_file_chunk'+obsids[obs_elements[0]]+'.sh'
	metafits_commands_file = open(metafits_script_path, 'w')
	
	#Write the contents of the file and the necessary arguments
	metafits_commands_file.write('#!/bin/bash\n\n' + \
		'#$ -S /bin/bash\n\n' + \
		'save_paths_metafits=(0 '+" ".join(metafits_path) + ')\n' + \
		'obsids_metafits=(0 ' +  " ".join([obsids[obs_element] for obs_element in obs_elements]) + ')\n' + \
		'ls ${save_paths_metafits[$SGE_TASK_ID]} > /dev/null\n' + \
		python_path.strip('\n') + ' ' + make_metafits_path + ' -o ${save_paths_metafits[$SGE_TASK_ID]} --gps ${obsids_metafits[$SGE_TASK_ID]} ')
		
	#Close the file
	metafits_commands_file.close() 

	#Make the file executable 
	os.chmod(metafits_script_path, 0775)

	qsub_command = "qsub -l h_vmem=1G,h_stack=512,h_rt=01:00:00 -V -pe chost 1 -e " + log_path + " -o " + log_path +" -N make_metafits -t 1:" + str(len(obs_elements))

	#Run the metafits script to create a metafits for the obsid in the uvfits folder
	stdoutpointer = subprocess.Popen((qsub_command + ' ' + metafits_script_path).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()
	task_jobid = (stdout_data.split())[2].rsplit('.')[0]
	
	#If there is data in the standard error output from the metafits script, let the user know
	if stderr_data:
		print 'WARNING: Grid Engine threw an error trying to run make_metafits.py'
		print stderr_data
		
	return task_jobid
	
		
def fill_database(obsids, save_paths, version, subversion, db_comment):

	#Find the cotter version	
	stdoutpointer = subprocess.Popen(["which","cotter"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()	
	cotter_path = stdout_data.strip('\n')
	stdoutpointer = subprocess.Popen((cotter_path + " --version").split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_data, stderr_data = stdoutpointer.communicate()
	cotter_version_string = stdout_data.splitlines()[0]
	cotter_version = cotter_version_string[cotter_version_string.rfind("version"):cotter_version_string.rfind(".")]

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

	for iteration, obsid in enumerate(obsids):
		print "Filling the mwa_qc database for obsid " + obsid
		save_path = save_paths[iteration] + obsid + '/'

		#Check to make sure the uvfits and metafits specified exist
		if not os.path.isfile(save_path + obsid + '.uvfits'):
			print "ERROR: " + save_path + obsid + ".uvfits does not exist! Database not updated"
			continue
		if not os.path.isfile(save_path + obsid + '.metafits'):
			print "WARNING: " + save_path + obsid + ".metafits does not exist! Database not updated"
			continue

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
			print "WARNING: A uvfits file for obsid " + obsid + ", version " + str(version) + ", subversion " + str(subversion) + \
				", cotter " + cotter_version + ", and frequency range " + bottom_freq_mhz + "-" + top_freq_mhz + " already exists."

		#Create the database row, and fill it with the inputs. 
		cur.execute("INSERT INTO uvfits(obsid,version,subversion,path,cotter_version,timestamp,comment,bottom_freq_mhz,top_freq_mhz) \
			VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);", \
			(obsid,version,subversion,save_path + obsid + '.uvfits',cotter_version,timestamp,db_comment,bottom_freq_mhz,top_freq_mhz))

	#Commit all the cur.execute, and close the connection.
	conn.commit()
	cur.close()
	conn.close()

	#Print completion message
	print "Filled the uvfits table in the mwa_qc database with the current uvfits information."

	
	
if __name__ == '__main__':
	main()
