#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():

	xlocs = [i for i in range(50)]
	ylocs = [i for i in range(50)]
	vx = [[1 for i in range(len(xlocs))] for j in range(len(ylocs))]
	vy = [[1 for i in range(len(xlocs))] for j in range(len(ylocs))]
	norm_factor = 2.
	
	x_vals = [xlocs]*len(ylocs)
	x_vals = [i for sublist in x_vals for i in sublist] #flatten list
	y_vals = [[i]*len(xlocs) for i in ylocs]
	y_vals = [i for sublist in y_vals for i in sublist] #flatten list
	
	plt.figure()
	#plt.plot(x_vals, y_vals, '.', markersize=2)
	for i in range(len(xlocs)):
		for j in range(len(ylocs)):
			vx_val = (vx[j])[i]/norm_factor
			vy_val = (vy[j])[i]/norm_factor
			length = (vx_val**2+vy_val**2)**0.5
			plt.plot([xlocs[i]-vx_val,xlocs[i]+vx_val],[ylocs[j]-vy_val,ylocs[j]+vy_val],color='black',lw=.3)
	plt.gca().set_aspect('equal',adjustable='box')
	plt.savefig('/nfs/eor-00/h1/rbyrne/drapery_plotting/vector_plot.png')


if __name__ == '__main__':
	main()
