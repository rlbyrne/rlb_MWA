#!/usr/bin/python

#script that finds the obsids closest to a given RA/Dec

def main():

	RA_target_h = 5
	RA_target_m = 35
	RA_target_s = 17.3
	Dec_target_deg = -5
	Dec_target_m = 23
	Dec_target_s = 28
	
	RA_target = (RA_target_h + RA_target_m/60. + RA_target_s/3600.)/24.*360.
	Dec_target = Dec_target_deg + Dec_target_m/60. + Dec_target_s/3600.

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
	
	distance_2 = []	
	for i, obsid in enumerate(obsids):
		distance_2.append((Decs[i]-Dec_target)**2+(RAs[i]-RA_target)**2)
		
	print "Closest obsid: " + str(obsids[distance_2.index(min(distance_2))]) + ", distance = " + str(min(distance_2)**(.5))
	obsids.remove(obsids[distance_2.index(min(distance_2))])
	distance_2.remove(min(distance_2))
	print "Second closest obsid: " + str(obsids[distance_2.index(min(distance_2))]) + ", distance = " + str(min(distance_2)**(.5))
	obsids.remove(obsids[distance_2.index(min(distance_2))])
	distance_2.remove(min(distance_2))
	print "Third closest obsid: " + str(obsids[distance_2.index(min(distance_2))]) + ", distance = " + str(min(distance_2)**(.5))
		
if __name__ == '__main__':
	main()
