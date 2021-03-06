#!/usr/bin/python

import sys

# script that finds the obsids closest to a given RA/Dec
# also gives information about the pointings of those obsids
# Obsolete: see find_obs.py

def main():

	RA_target_h = 0
	RA_target_m = 0
	RA_target_s = 0
	Dec_target_deg = -27
	Dec_target_m = 0
	Dec_target_s = 0

	RA_target = (RA_target_h + RA_target_m/60. + RA_target_s/3600.)/24.*360.
	Dec_target = Dec_target_deg + Dec_target_m/60. + Dec_target_s/3600.

	print "RA: " + str(RA_target)
	print "Dec: " + str(Dec_target)

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

	pointings_codes = [0]*len(obsids)
	pointings_code_options = [1,2,2,3,3]
	for band_index, dec_band in enumerate(Decs_set):
		Azimuths_band = []
		Elevations_band = []
		for obs_index in range(len(obsids)):
			if Decs_round[obs_index] == dec_band:
				Azimuths_band.append(Azimuths[obs_index])
				Elevations_band.append(Elevations[obs_index])
		AzEls_band = zip(Azimuths_band, Elevations_band)
		AzEls_band_set = list(set(AzEls_band))
		Elevations_band_set = [value[1] for value in AzEls_band_set]
		Elevations_band_set_sort = list(reversed(sorted(Elevations_band_set)))
		for obs_index in range(len(obsids)):
			if Decs_round[obs_index] == dec_band:
				pointings_codes[obs_index] = pointings_code_options[Elevations_band_set_sort.index(Elevations[obs_index])]

	obsids_sorted = [[],[],[]]
	distance_2 = [[],[],[]]
	for i, obsid in enumerate(obsids):
		use_distance_2 = min([(Decs[i]-Dec_target)**2+(use_RAs-RA_target)**2 for use_RAs in [RAs[i]-360,RAs[i],RAs[i]+360]])
		obsids_sorted[pointings_codes[i]-1].append(obsid)
		distance_2[pointings_codes[i]-1].append(use_distance_2)

	for i, obs_set in enumerate(obsids_sorted):
		distance_2_set = distance_2[i]
		print "......."
		if i == 0:
			print "For the BEST pointings:"
		if i == 1:
			print "For the SECOND BEST pointings:"
		if i == 2:
			print "For the WORST pointings:"
		print "Closest obsid: " + str(obs_set[n_min(distance_2_set,1)]) + ", distance = " + str(distance_2_set[n_min(distance_2_set,1)]**0.5)
		print "Second closest obsid: " + str(obs_set[n_min(distance_2_set,2)]) + ", distance = " + str(distance_2_set[n_min(distance_2_set,2)]**0.5)
		print "Third closest obsid: " + str(obs_set[n_min(distance_2_set,3)]) + ", distance = " + str(distance_2_set[n_min(distance_2_set,3)]**0.5)
		print "Fourth closest obsid: " + str(obs_set[n_min(distance_2_set,4)]) + ", distance = " + str(distance_2_set[n_min(distance_2_set,4)]**0.5)
		print "Fifth closest obsid: " + str(obs_set[n_min(distance_2_set,5)]) + ", distance = " + str(distance_2_set[n_min(distance_2_set,5)]**0.5)


def n_min(data,n):
	#returns the nth minimum of a list data
	use_data = [element for element in data]
	for iteration in range(n):
		minimum = min(use_data)
		use_data.remove(minimum)
	i = data.index(minimum)
	return i

if __name__ == '__main__':
	main()
