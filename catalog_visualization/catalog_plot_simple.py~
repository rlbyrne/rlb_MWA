#!/usr/bin/python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	catalog_name = '1130781304_run1_catalog'
	#catalog_name = 'GLEAM'
	catalog_data_path = '/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'.txt'
	
	flux_cut = 0
	ra_plot_range = [-150,210]
	dec_plot_range = [-90,90]
	markersize_range = [.00001,50]
	flux_plot_max = 50

	datafile = open(catalog_data_path, "r")
	catalog_info = datafile.readlines()
	datafile.close()
	catalog_info = catalog_info[1:] #remove header

	ra = []
	dec = []
	marker_sizes = []
	
	for info in catalog_info:
		info = info.split(",")
		ra_val = float(info[0])
		dec_val = float(info[1])
		flux_val = float(info[2])
		if ra_val >= ra_plot_range[1]:
			ra_val -= 360.
		if ra_val >= ra_plot_range[0] and ra_val <= ra_plot_range[1] and dec_val >= dec_plot_range[0] and dec_val <= dec_plot_range[1] and flux_val >= flux_cut:
			ra.append(ra_val)
			dec.append(dec_val)
			if flux_val > flux_plot_max:
				flux_val = flux_plot_max
			marker_size_val = flux_val/flux_plot_max*(markersize_range[1]-markersize_range[0])+markersize_range[0]
			marker_sizes.append(marker_size_val)

	plt.figure()
	plt.scatter(ra,dec,marker='o',s=marker_sizes)
	#plt.axis('equal')
	plt.xlim(ra_plot_range[0],ra_plot_range[1])
	plt.ylim(dec_plot_range[0],dec_plot_range[1])
	plt.gca().set_aspect('equal',adjustable='box')
	plt.savefig('/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'_plot')
		
if __name__ == '__main__':
	main()
