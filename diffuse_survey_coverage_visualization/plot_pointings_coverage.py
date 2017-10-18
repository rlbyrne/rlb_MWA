#!/usr/bin/python

# OBSOLETE CODE: see survey_plotter.py function plot_radec_pointings_coverage

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	ra_range = [-60,160]
	dec_range = [-70,20]
	resolution = 1
	radius = 12
	colorbar_max = 10

	ra_vals = [ra_range[0]+resolution*i for i in range(int((ra_range[1]-ra_range[0])/resolution))]
	dec_vals = [dec_range[0]+resolution*i for i in range(int((dec_range[1]-dec_range[0])/resolution))]
	number_obsids = [[0 for i in range(len(ra_vals))] for i in range(len(dec_vals))]

	obsfile_name = "/nfs/eor-00/h1/rbyrne/sidelobe_survey_obsinfo.txt"
	obsfile = open(obsfile_name, "r")
	obsinfo = [line.split("\n") for line in obsfile.readlines()]
	obsfile.close()
	obsinfo = [obs[0] for obs in obsinfo]
	obsinfo = list(set(obsinfo))
	obsinfo.sort()

	#obsids = []
	#LSTs = []
	RAs = []
	Decs = []
	Azimuths = []
	Elevations = []
	#AzEls = []

	for info in obsinfo:
		info = info.split(", ")
		#obsids.append(int(info[0]))
		#LSTs.append(float(info[1]))
		if float(info[2]) < 160:
			RAs.append(float(info[2]))
		else:
			RAs.append(float(info[2])-360.)
		Decs.append(float(info[3]))
		Azimuths.append(round(float(info[4])))
		Elevations.append(round(float(info[5])))

	AzEls = zip(Azimuths,Elevations)

	for i, RA_target in enumerate(ra_vals):
		for j, Dec_target in enumerate(dec_vals):
			(number_obsids[-(j+1)])[-(i+1)] = check_obsids(RAs, Decs, RA_target, Dec_target, AzEls, radius, colorbar_max)

	plt.figure(figsize=(9,3))
	plt.imshow(number_obsids, interpolation='none', extent=[ra_range[1]/360.*24.,ra_range[0]/360.*24.,dec_range[0],dec_range[1]], aspect=24/360., cmap=plt.get_cmap('Blues'))
	plt.xlabel('Right Ascension (hours)')
	plt.ylabel('Declination (degrees)')
	plt.tight_layout()
	cbar = plt.colorbar(extend='max')
	cbar.set_label('Number of Unique Pointings')
	plt.savefig("/nfs/eor-00/h1/rbyrne/radec_plots/pointings_coverage.png")


def check_obsids(RAs, Decs, RA_target, Dec_target, AzEls, radius, colorbar_max):

	pointings = []
	for i in range(len(RAs)):
		distance_2 = (Decs[i]-Dec_target)**2+(RAs[i]-RA_target)**2
		if distance_2 <= radius**2:
			pointings.append(AzEls[i])
	number_pointings = len(set(pointings))
	if number_pointings > colorbar_max:
		number_pointings = colorbar_max
	return number_pointings

if __name__ == '__main__':
	main()
