#!/usr/bin/python

#Script that deletes downloaded data based on their cotter version and subverion numbers

import os

def main():

	version = 0
	subversion = 0

	all_nodes = ["eor-00", "eor-02", "eor-03", "eor-04", "eor-05", "eor-06", "eor-07", "eor-08", "eor-09", "eor-10", "eor-11", "eor-12", "eor-13", "eor-14"]
	all_nodes = ["/nfs/" + nodename + "/r1/EoRuvfits/" for nodename in all_nodes]
	
	dir_deleted = False
	for node_path in all_nodes:
		if os.path.isdir(node_path):
			directory_contents = os.listdir(node_path)
			for dirname in directory_contents:
				if dirname.endswith("v" + str(version) + "_" + str(subversion)):
					print dirname + " FOUND IN DIRECTORY " + node_path + " ... DELETING ..."
					os.system("rm " + node_path + dirname + " -r")
					dir_deleted = True
	if not dir_deleted:
		print "NO FILES FOUND"
	
if __name__ == '__main__':
	main()
