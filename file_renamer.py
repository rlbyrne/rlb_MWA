#!/usr/bin/python

import os
import sys
import time

def main():

	filepath_start = "/nfs/mwa-08/d1/DiffuseSurvey2015"
	run_name = "fhd_rlb_diffuse_survey_threeobs_nodiffuse"
	replacement_string = "sidelobe_survey_obsIDs_firstobs"
	initial_string = "diffuse_survey_threeobs_noneor0"

	filepaths = [filepath_start + "/" + run_name + directory for directory in ["/Healpix","/ps"]]
	for path in filepaths:
		filenames = os.listdir(path)
		for filename in filenames:
			string_start = filename.find(initial_string)
			if string_start != -1:
				replacement_filename = filename.replace(initial_string, replacement_string)
				print "RENAMING " + filename + " AS " + replacement_filename
				#time.sleep(.1)
				os.rename(path + "/"+ filename, path + "/"+replacement_filename)

if __name__ == '__main__':
	main()
