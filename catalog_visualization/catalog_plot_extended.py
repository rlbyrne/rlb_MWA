#!/usr/bin/python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	#catalog_name = '1130781304_run1_catalog'
	catalog_name = 'GLEAM'
	catalog_data_path = '/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'.txt'
	
	flux_cut = 0
	ra_plot_range = [-150,210]
	dec_plot_range = [-90,90]
	markersize_range = [.00001,50]
	flux_plot_max = 200

	datafile = open(catalog_data_path, "r")
	catalog_info = datafile.readlines()
	datafile.close()
	catalog_info = catalog_info[1:] #remove header

	ra_sources = []
	dec_sources = []
	marker_sizes_sources = []
	ra_components = []
	dec_components = []
	marker_sizes_components = []
	
	for info in catalog_info:
		info = info.split(",")
		ra_val = float(info[0])
		dec_val = float(info[1])
		flux_val = float(info[2])
		extended = int(info[3])
		if ra_val >= ra_plot_range[1]:
			ra_val -= 360.
		if ra_val >= ra_plot_range[0] and ra_val <= ra_plot_range[1] and dec_val >= dec_plot_range[0] and dec_val <= dec_plot_range[1] and flux_val >= flux_cut:
			if extended == 0:
				ra_sources.append(ra_val)
				dec_sources.append(dec_val)
				if flux_val > flux_plot_max:
					flux_val = flux_plot_max
				marker_size_val = flux_val/flux_plot_max*(markersize_range[1]-markersize_range[0])+markersize_range[0]
				marker_sizes_sources.append(marker_size_val)
			else:
				num_components = int(info[4])
				for i in range(num_components):
					ra_components.append(float(info[5+3*i]))
					dec_components.append(float(info[6+3*i]))
					marker_sizes_components.append(float(info[7+3*i])/flux_plot_max*(markersize_range[1]-markersize_range[0])+markersize_range[0])
					
	plt.figure()
	plt.scatter(ra_sources,dec_sources,marker='o',s=marker_sizes_sources)
	plt.scatter(ra_components,dec_components,marker='o',s=marker_sizes_components,c='red')
	if False:
		a_team_ras = [83.6331,79.9572,139.524,201.365,252.784,187.706,299.868,350.858] #[Crab,Pic A,Hydra A,Cen A,Her A,Vir A,Cygnus A,Cas A]
		for i, ra in enumerate(a_team_ras):
			if ra > ra_plot_range[1]:
				a_team_ras[i] = ra-360.
		a_team_decs = [22.0145,-45.7788,-12.0956,-43.0192,4.9925,12.3911,40.7339,58.8]
		plt.scatter(a_team_ras,a_team_decs,marker='o',s=markersize_range[1],c='orange')
		a_team_names = ['Crab','Pic A','Hydra A','Cen A','Her A','Vir A','Cygnus A','Cas A']
		for i, name in enumerate(a_team_names):
			plt.annotate(name, (a_team_ras[i],a_team_decs[i]))
	plt.xlim(ra_plot_range[0],ra_plot_range[1])
	plt.ylim(dec_plot_range[0],dec_plot_range[1])
	plt.xlabel('RA (deg)')
	plt.ylabel('Dec (deg)')
	plt.gca().set_aspect('equal',adjustable='box')
	plt.savefig('/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'_extended_plot')
		
if __name__ == '__main__':
	main()
