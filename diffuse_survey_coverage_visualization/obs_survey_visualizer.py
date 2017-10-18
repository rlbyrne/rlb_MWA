#!/usr/bin/python

# OBSOLETE CODE: SEE survey_plotter.py functions plot_azels,
# plot_radecs_colorcode_decs, and generate_radec_animation

import matplotlib.pyplot as plt
import sys

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
		AzEls.append(str(int(float(info[4]))) + ", " + str(int(float(info[5]))))

	Decs_round = [int(Dec/6.)*6 for Dec in Decs] #round the declinations to the nearest 6th to clump them in 7 Dec bands
	Decs_set = list(set(Decs_round))
	Decs_set.sort()

	use_colors = ["black","red","green","magenta","cyan","yellow","blue","black","red","green","magenta","cyan","yellow"]

	if False:
		plt.figure()
		for i in range(len(Decs_set)):
			use_Az = []
			use_El = []
			for index in range(len(RAs)):
				if Decs_round[index] == Decs_set[i]:
					use_Az.append(Azimuths[index])
					use_El.append(Elevations[index])
			use_AzEls = zip(use_Az, use_El)
			use_AzEls_set = list(set(use_AzEls))
			print use_AzEls_set
			use_Az_set = [term[0] for term in use_AzEls_set]
			use_El_set = [term[1] for term in use_AzEls_set]
			plt.plot(use_Az_set,use_El_set, "o", markersize=10, mfc=use_colors[i], alpha = 1)
		plt.xlabel("Azimuth")
		plt.ylabel("Elevation")
		plt.axis([-10,350,60,91])
		plt.grid(True)
		plt.savefig("/nfs/eor-00/h1/rbyrne/radec_plots/AzEls.png")


	for i, RA in enumerate(RAs):
		if RA > 250:
			RAs[i] = RAs[i] - 360

	if False:
		save_filepath = "/nfs/eor-00/h1/rbyrne/radec_plots/dec_colorcode_plot.png"
		plt.figure(figsize=(10,5))
		for i in range(len(Decs_set)):
			use_RAs = []
			use_Decs = []
			for index in range(len(RAs)):
				if Decs_round[index] == Decs_set[i]:
					use_RAs.append(RAs[index])
					use_Decs.append(Decs[index])
			plt.plot(use_RAs,use_Decs, "o", markersize=10, mfc=use_colors[i], alpha = 0.5)
		plt.xlabel("RA")
		plt.ylabel("Dec")
		plt.axis([-100,200,-65,15])
		plt.grid(True)
		plt.savefig(save_filepath)


	if True:
		for i, obs in enumerate(obsids):
		#if True:
			#obs = obsids[-1]
			#i = len(obsids)-1
			if i+1 < 10:
				filepath_num = "00" + str(i+1)
			else:
				if i+1 < 100:
					filepath_num = "0" + str(i+1)
				else:
					filepath_num = str(i+1)
			save_filepath = "/nfs/eor-00/h1/rbyrne/radec_plots/radec_plot" + filepath_num + ".png"
			#save_filepath = "/nfs/eor-00/h1/rbyrne/radec_plots/radec_last_plot.png"
			plt.figure(figsize=(17,5))
			plt.plot(RAs[0:i+1],Decs[0:i+1], "o", markersize=90, mfc="blue", alpha = 0.03)
			plt.plot(RAs[i],Decs[i], "o", markersize=90, mfc="none", mec="red")
			plt.plot(RAs[i],Decs[i], "x", markersize=10, mfc="red")
			plt.xticks(range(-100,200,10))
			plt.xlabel("RA")
			plt.ylabel("Dec")
			plt.axis("equal")
			plt.axis([-100,200,-65,15])
			plt.grid(which="both")
			plt.text(-80, 17, "Obsid: " + str(obs))
			plt.savefig(save_filepath)
			plt.close()
			#if i == 50:
			#	break



if __name__ == '__main__':
	main()
