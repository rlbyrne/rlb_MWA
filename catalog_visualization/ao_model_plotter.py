#!/usr/bin/python

#Script that plots the AO-compatible GLEAM extended source models.

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ao_model_reader as ao_read

def main():

	name = 'CenA-50comp'
	source_models = ao_read.read_model(name)
	
	ras = []
	decs = []
	fluxes = []
	for i, source in enumerate(source_models):
		comps = source.get('components')
		for j, comp in enumerate(comps):
			ras.append(comp.get('ra_deg'))
			decs.append(comp.get('dec_deg'))
			fluxes.append(comp.get('flux_I_Jy'))
	
	flux_cut = 50
	markersize_range = [1,100]
	flux_plot_max = max(fluxes)

	marker_sizes = []
	
	for flux_val in fluxes:
		if flux_val > flux_plot_max:
			flux_val = flux_plot_max
		marker_size_val = flux_val/flux_plot_max*(markersize_range[1]-markersize_range[0])+markersize_range[0]
		marker_sizes.append(marker_size_val)

	plt.figure()
	plt.scatter(ras,decs,marker='o',s=marker_sizes,alpha=1)
	plt.xlabel('RA (deg)')
	plt.ylabel('Dec (deg)')
	plt.gca().set_aspect('equal',adjustable='box')
	plt.grid()
	print 'Saving plot to /nfs/eor-00/h1/rbyrne/catalog_data/'+name+'_plot.png'
	plt.savefig('/nfs/eor-00/h1/rbyrne/catalog_data/'+name+'_plot')
	
	
if __name__ == '__main__':
	main()
