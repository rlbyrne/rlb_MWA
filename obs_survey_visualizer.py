#!/usr/bin/python

import matplotlib.pyplot as plt

def main():

	obsfile_name = "/nfs/eor-00/h1/rbyrne/sidelobe_survey_obsinfo.txt"
	obsfile = open(obsfile_name, "r")
	obsinfo = [line.split("\n") for line in obsfile.readlines()]
	obsfile.close()
	obsinfo = [obs[0] for obs in obsinfo]
	obsinfo = list(set(obsinfo))
	obsinfo.sort()

	obsids = []
	LSTs = []
	RAs = []
	Decs = []
	Azimuths = []
	Elevations = []
	
	for info in obsinfo:
		info = info.split(", ")
		obsids.append(info[0])
		LSTs.append(info[1])
		RAs.append(info[2])
		Decs.append(info[3])
		Azimuths.append(info[4])
		Elevations.append(info[5])
		
	save_filepath = "/nfs/eor-00/h1/rbyrne/python_practice/testplot.png"
	plt.figure()
	plt.plot([1,2,3],[1,2,3], "o", markersize=20, mfc="none")
	plt.axis([0,5,0,5])
	plt.grid(True)
	plt.savefig(save_filepath)



if __name__ == '__main__':
	main()
