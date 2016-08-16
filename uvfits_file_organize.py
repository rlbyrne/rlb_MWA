#!/usr/bin/python

import os
from astropy.time import Time
import sys

def main():

	original_path = "/nfs/eor-05/r1/EoRuvfits/transfer"
	destination_node = "/nfs/eor-05/r1"
	version = 5 
	subversion = 1
	
	if not os.path.isdir(original_path):
		print "Source directory does not exit."
		sys.exit(1)
	directory_contents = os.listdir(original_path)
	print len(directory_contents)
	for filename in directory_contents:
		if filename.endswith(".uvfits"):
			obsid = filename[0:10]
			t = Time(int(obsid), format="gps", scale="utc")
			jd = t.jd
			jd = int(jd)
			save_directory = "/EoRuvfits/jd" + str(jd) + "v"+ str(version) + "_" + str(subversion) + "/"
			if not os.path.isdir(destination_node + save_directory + obsid):
				if not os.path.isdir(destination_node + save_directory):
					os.system("mkdir " + destination_node + save_directory)
				os.system("mkdir " + destination_node + save_directory + obsid)
			os.system("mv " + original_path + "/" + filename + " " + destination_node + save_directory + obsid)
			print "moving file " + filename + " to " + destination_node + save_directory + obsid
	
	
	
if __name__ == '__main__':
	main()
