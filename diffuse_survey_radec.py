#!/usr/bin/python

import ephem
from astropy.time import Time
import numpy as np
import math

def main():
	
	obsfile_name = "/nfs/eor-00/h1/rbyrne/sidelobe_survey_obsIDs.txt"
	obsfile = open(obsfile_name, "r")
	obsids = [line.split( ) for line in obsfile.readlines()]
	obsids = [obs[0] for obs in obsids]
	obsfile.close()

	#t = Time([int(obsid) for obsid in obsids], format="gps", scale="utc")
	t = Time([1130772904], format = "gps", scale = "utc")
	jdates = t.jd
	jdate = jdates[0]

	az = np.radians(198.435)
	el = np.radians(67.98)

	#approximate values, fix these:
	lon = np.radians(116.6)
	lat = np.radians(-26.75)
	alt = 378

	jdate0 = ephem.julian_date(0)

	observer = ephem.Observer()
	observer.lon = lon
	observer.lat = lat
	observer.elevation = alt
	observer.date = jdate - jdate0

	lst = observer.sidereal_time()
	lst_hrs = float(lst)*23.9344699/(2*math.pi)
	print lst_hrs
	print lst
	ra,dec = observer.radec_of(az, el)

	print "EPHEM: %f %f" % (np.degrees(ra),np.degrees(dec))



if __name__ == '__main__':
	main()
