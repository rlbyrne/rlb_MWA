#!/usr/bin/python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	#catalog_name = '1130781304_run1_catalog'
	catalog_name = 'GLEAM_plus_extended_fornax'
	catalog_data_path = '/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'.txt'
	
	saveplot_name = catalog_name
	
	flux_cut = 0
	ra_plot_range = [48,52]
	dec_plot_range = [-39,-35]
	#ra_plot_range = [137,141]
	#dec_plot_range = [-14,-10]
	markersize_range = [1,20]
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
					ra_use = float(info[5+3*i])
					if ra_use >= ra_plot_range[1]:
						ra_use -= 360.
					ra_components.append(ra_use)
					dec_components.append(float(info[6+3*i]))
					marker_sizes_components.append(float(info[7+3*i])/flux_plot_max*(markersize_range[1]-markersize_range[0])+markersize_range[0])
					
	plt.figure()
	plt.scatter(ra_sources,dec_sources,marker='o',s=marker_sizes_sources)	
	plt.scatter(ra_components,dec_components,marker='o',s=marker_sizes_components,c='red')
	plt.grid()
	if False:
		a_team_ras = [83.6331,79.9572,139.524,201.365,252.784,187.706,299.868,350.858]
		for i, ra in enumerate(a_team_ras):
			if ra > ra_plot_range[1]:
				a_team_ras[i] = ra-360.
		a_team_decs = [22.0145,-45.7788,-12.0956,-43.0192,4.9925,12.3911,40.7339,58.8]
		plt.scatter(a_team_ras,a_team_decs,marker='o',s=markersize_range[1],c='orange')
		a_team_names = ['Crab','Pic A','Hydra A','Cen A','Her A','Vir A','Cygnus A','Cas A']
		for i, name in enumerate(a_team_names):
			plt.annotate(name, (a_team_ras[i],a_team_decs[i]))
	if False:
		other_bright_ras = [50.6738,128.836,83.8221]
		for i, ra in enumerate(other_bright_ras):
			if ra > ra_plot_range[1]:
				other_bright_ras[i] = ra-360.
		other_bright_decs = [-37.2083,-45.1764,-5.3911]
		plt.scatter(other_bright_ras,other_bright_decs,marker='o',s=markersize_range[1],c='cyan')
		other_bright_names = ['Fornax A','Vela','Orion']
		for i, name in enumerate(other_bright_names):
			plt.annotate(name, (other_bright_ras[i],other_bright_decs[i]))
	if False:
		other_bright_ras = [201.437,96.7921,260.117,303.615,333.607,350.866,50.6738,252.793,139.524,11.8991,79.9541,187.706,57.8988,71.1571,329.275,353.609,359.768,22.6158]
		for i, ra in enumerate(other_bright_ras):
			if ra > ra_plot_range[1]:
				other_bright_ras[i] = ra-360.
		other_bright_decs = [-42.9608,-5.88472,-0.979722,23.5814,-17.0267,58.8117,-37.2083,4.99806,-12.0831,-25.2886,-45.7649,12.3786,-27.7431,-28.1653,-69.6900,-41.4233,-60.9164,-26.1656]
		plt.scatter(other_bright_ras,other_bright_decs,marker='o',s=markersize_range[1]+1,c='yellow')
		other_bright_names = ['CenA','3C161','3C353','3C409','3C444','CasA','FornaxA','HerA','HydA','NGC0253','PicA','VirA','PKS0349-27','PKS0442-28','PKS2153-69','PKS2331-41','PKS2356-61','PKSJ0130-2610']
		for i, name in enumerate(other_bright_names):
			plt.annotate(name, (other_bright_ras[i],other_bright_decs[i]))


	plt.xlim(ra_plot_range[0],ra_plot_range[1])
	plt.ylim(dec_plot_range[0],dec_plot_range[1])
	plt.xlabel('RA (deg)')
	plt.ylabel('Dec (deg)')
	plt.gca().set_aspect('equal',adjustable='box')
	print 'Saving plot to /nfs/eor-00/h1/rbyrne/catalog_data/'+saveplot_name+'_plot.png'
	plt.savefig('/nfs/eor-00/h1/rbyrne/catalog_data/'+saveplot_name+'_plot')
		
if __name__ == '__main__':
	main()
