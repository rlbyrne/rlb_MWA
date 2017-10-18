#!/usr/bin/python

import sys

# Script that finds the pointings of observations in the diffuse survey
# Obsolete: see surveyview.py module

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
	if len(Decs_set) != 7:
		print "ERROR: Dividing observations into declination bands failed"
		sys.exit(1)

	Azimuths = [round(term) for term in Azimuths]
	Elevations = [round(term) for term in Elevations]

	pointings_codes = ['']*len(obsids)
	dec_pointings_options = [3,2,1,0,-1,-2,-3]
	azimuth_pointings_options = [2,1,0,-1,-2]
	for band_index, dec_band in enumerate(Decs_set):
		Azimuths_band = []
		Elevations_band = []
		for obs_index in range(len(obsids)):
			if Decs_round[obs_index] == dec_band:
				Azimuths_band.append(Azimuths[obs_index])
				Elevations_band.append(Elevations[obs_index])
		AzEls_band = zip(Azimuths_band, Elevations_band)
		AzEls_band_set = list(set(AzEls_band))
		Azimuth_band_set = [value[0] for value in AzEls_band_set]
		Elevations_band_set = [value[1] for value in AzEls_band_set]
		Az_zero = Azimuth_band_set[Elevations_band_set.index(max(Elevations_band_set))]
		Azimuth_band_set_sort = list(sorted(Azimuth_band_set))
		Az_wrap = 0
		for Az in Azimuth_band_set_sort:
			if Az > Az_zero + 180:
				Az_wrap += 1
		if Az_wrap > 0:
			Azimuth_band_set_sort = Azimuth_band_set_sort[len(Azimuth_band_set_sort)-Az_wrap:] + Azimuth_band_set_sort[:len(Azimuth_band_set_sort)-Az_wrap]

		for obs_index in range(len(obsids)):
			if Decs_round[obs_index] == dec_band:
				pointings_codes[obs_index] = '({},{})'.format(azimuth_pointings_options[Azimuth_band_set_sort.index(Azimuths[obs_index])],dec_pointings_options[band_index])

	outfile_name = '/nfs/eor-00/h1/rbyrne/diffuse_survey_pointings.csv'
	print 'saving data to ' + outfile_name
	outfile = open(outfile_name, 'w')
	outfile.write('obsid,pointing\n')
	for i, obsid in enumerate(obsids):
		outfile.write(str(obsid)+','+pointings_codes[i]+'\n')
	outfile.close()

if __name__ == '__main__':
	main()
