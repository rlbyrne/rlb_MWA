#!/usr/bin/python

import matplotlib.pyplot as plt
import sys

def main():

	datafile_name = "/nfs/eor-00/h1/rbyrne/ram_log.txt"
	datafile = open(datafile_name, "r")
	raminfo = [line.split("\n") for line in datafile.readlines()]
	datafile.close()
	
	times = []
	loads = []
	ram_use = []
	next_day = 0
	for i,datapoint in enumerate(raminfo):
		datapoint = datapoint[0]
		if i%2 == 0:
			datapoint = datapoint.split(" ")
			time_stamp = datapoint[1]
			time_stamp = time_stamp.split(":")
			time_stamp = int(time_stamp[0])+int(time_stamp[1])/60.+int(time_stamp[2])/3600.
			if len(times) > 0 and time_stamp + 24*next_day < times[-1]:
				next_day += 1
			time_stamp += 24*next_day
			times.append(time_stamp)
			load_1min = datapoint[-3]
			load_1min = float(load_1min.replace(",",""))
			loads.append(load_1min)
		else:
			datapoint = datapoint.split()
			ram_use.append(int(datapoint[2]))

	time_0 = times[0]
	times = [t-time_0 for t in times]
	
	plt.figure()
	plt.plot(times,loads,"o",markersize=1)
	plt.xlabel("time")
	plt.ylabel("load (1 min. ave.)")
	plt.xlim([0,10])
	#plt.show()
	plt.savefig("/nfs/eor-00/h1/rbyrne/load_plot.png")

	plt.figure()
	plt.plot(times,ram_use,"-")
	plt.xlabel("time")
	plt.ylabel("mem use (MB)")
	plt.xlim([0,10])
	plt.savefig("/nfs/eor-00/h1/rbyrne/ram_use_plot.png")

if __name__ == '__main__':
	main()
