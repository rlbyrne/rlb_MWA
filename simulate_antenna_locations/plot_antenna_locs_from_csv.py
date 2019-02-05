#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_antenna_locs_from_csv():

    csv_path = '/Users/rubybyrne/array_simulation_331/split_hex_array_sim_331_antenna_locs.csv'
    output_path = '/Users/rubybyrne/array_simulation_331/antenna_layout_plots_small'

    csv_file = open(csv_path, 'r')
    csv_contents = csv_file.readlines()
    csv_file.close()
    antenna_xlocs = np.zeros(len(csv_contents)-1)
    antenna_ylocs = np.zeros(len(csv_contents)-1)
    for ind, line in enumerate(csv_contents[1:]):
        line_split = line.split(',')
        antenna_xlocs[ind] = float(line_split[1])
        antenna_ylocs[ind] = float(line_split[2])
    random_amp = np.random.normal(loc=1., scale=.01, size=len(csv_contents)-1)

    plt.figure()
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.rcParams.update({'font.size': 13})
    cmap = plt.cm.get_cmap('plasma')
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=3.)
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp, cmap=cmap, vmin=.96, vmax=1.04, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    #cbar = plt.colorbar(sc)
    #cbar.ax.set_ylabel('Antenna Gain Amplitude', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/split_hex_array_331_layout_plot.png'.format(output_path))
    plt.close()


if __name__ == '__main__':
    plot_antenna_locs_from_csv()
