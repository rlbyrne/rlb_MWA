#!/usr/bin/python

import matplotlib.pyplot as plt
import sys
import numpy
import math

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
	AzEls = []
	
	for info in obsinfo:
		info = info.split(", ")
		obsids.append(int(info[0]))
		LSTs.append(float(info[1]))
		RAs.append(float(info[2]))
		Decs.append(float(info[3]))
		Azimuths.append(float(info[4]))
		Elevations.append(float(info[5]))
	
	Decs_round = [int(Dec/6.)*6 for Dec in Decs] #round the declinations to the nearest 6th to clump them in 7 Dec bands
	Decs_set = list(set(Decs_round))
	Decs_set.sort()
	
	Azimuths = [round(term) for term in Azimuths]
	Elevations = [round(term) for term in Elevations]
	

	for band_index, dec_band in enumerate(Decs_set):
		Decs_band = []
		RAs_band = []
		Azimuths_band =[]
		Elevations_band =[]
		Decs_other = []
		RAs_other = []
		Azimuths_other =[]
		Elevations_other =[]
	
		for index in range(len(obsids)):
			if Decs_round[index] == dec_band:
				Decs_band.append(Decs[index])
				RAs_band.append(RAs[index])
				Azimuths_band.append(Azimuths[index])
				Elevations_band.append(Elevations[index])
			else:
				Decs_other.append(Decs[index])
				RAs_other.append(RAs[index])
				Azimuths_other.append(Azimuths[index])
				Elevations_other.append(Elevations[index])
				
		AzEls_set_band = zip(Azimuths_band, Elevations_band)
		AzEls_set_band = list(set(AzEls_set_band))
		if len(AzEls_set_band) != 5:
			print "ERROR: obsids not sorted into 5 pointings"
			sys.exit(1)
		AzEls_set_band.sort(key=lambda tup: tup[0])
		Az_set_band = [term[0] for term in AzEls_set_band]
		El_set_band = [term[1] for term in AzEls_set_band]
		
		AzEls_set_other = zip(Azimuths_other, Elevations_other)
		AzEls_set_other = set(AzEls_set_other)
		Az_set_other = [term[0] for term in AzEls_set_other]
		El_set_other = [term[1] for term in AzEls_set_other]
		
		Decs_by_pointing = [[],[],[],[],[]]
		RAs_by_pointing = [[],[],[],[],[]]
		for pointing_index in range(len(Az_set_band)):
			for obs_index in range(len(Decs_band)):
				if Azimuths_band[obs_index] == Az_set_band[pointing_index] and Elevations_band[obs_index] == El_set_band[pointing_index]:
					(Decs_by_pointing[pointing_index]).append(Decs_band[obs_index])					
					(RAs_by_pointing[pointing_index]).append(RAs_band[obs_index])
		
		plot_points = 1000			
		RA_vals = [i/100.*160.-10. for i in range(plot_points)]
		how_many_obs = [0 for i in range(plot_points)]
		#how_many_pointings = [0 for i in range(plot_points)]
		for index, val in enumerate(RA_vals):
			in_pointing = [False, False, False, False, False]
			for RA in RAs_band:
			 	how_many_obs[index] = how_many_obs[index] + math.e**(-(val-RA)**2/(2*1.5**2))
			 		#for pointing in range(5):
			 			#if RA in RAs_by_pointing[pointing]:
			 			#	in_pointing[pointing] = True
			#how_many_pointings[index] = in_pointing.count(True)
		
		colors = ["red", "darkorange", "yellow", "green", "skyblue", "violet"]			
		plt.figure(figsize=(8,12))
		
		plt.subplot(4,1,1)
		plt.plot(Az_set_other, El_set_other, "o", markersize=10, mfc="black", alpha=1)
		for pointing_index in range(len(Az_set_band)):
			plt.plot(Az_set_band[pointing_index], El_set_band[pointing_index], "o", markersize=10, mfc=colors[pointing_index], alpha=1)
		plt.xlabel("Azimuth")
		plt.ylabel("Elevation")
		plt.axis([-10,350,60,91])
		plt.grid(True)
		
		plt.subplot(4,1,2)
		plt.plot(RAs_other, Decs_other, "o", markersize=15, mfc="black", alpha=0.05)
		for pointing_index in range(len(Az_set_band)):
			plt.plot(RAs_by_pointing[pointing_index], Decs_by_pointing[pointing_index], "o", markersize=15, mfc=colors[pointing_index], alpha=0.5)
		plt.xlabel("RA")
		plt.ylabel("Dec")
		plt.axis([-10,150,-65,15])
		plt.grid(True)
		
		plt.subplot(4,1,3)
		for pointing_index in range(len(Az_set_band)):
			plt.plot(RAs_by_pointing[pointing_index], [pointing_index+1 for i in range(len(RAs_by_pointing[pointing_index]))], \
				"o", markersize=15, mfc=colors[pointing_index], alpha=0.5)
		plt.xlabel("RA")
		plt.ylabel("Pointing")
		plt.axis([-10,150,0,6])
		plt.grid(True)
		
		plt.subplot(4,1,4)
		plt.plot(RA_vals, how_many_obs, "b-")
		#plt.plot(RA_vals, how_many_pointings, "r-")
		plt.xlabel("RA")
		plt.ylabel("Pointing")
		plt.axis([-10,150,0,15])
		
		
		plt.savefig("/nfs/eor-00/h1/rbyrne/radec_plots/pointings_plot" + str(band_index+1) + ".png")


if __name__ == '__main__':
	main()
