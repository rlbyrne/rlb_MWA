#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_antenna_locs_from_csv():

    csv_path = '/Users/ruby/EoR/hex_array_sim_331_antenna_locs.csv'
    output_path = '/Users/ruby/EoR'

    csv_file = open(csv_path, 'r')
    csv_contents = csv_file.readlines()
    csv_file.close()
    antenna_xlocs = np.zeros(len(csv_contents)-1)
    antenna_ylocs = np.zeros(len(csv_contents)-1)
    for ind, line in enumerate(csv_contents[1:]):
        line_split = line.split(',')
        antenna_xlocs[ind] = float(line_split[1])
        antenna_ylocs[ind] = float(line_split[2])
    print(antenna_xlocs)
    random_amp = np.random.normal(loc=1., scale=.01, size=len(csv_contents)-1)

    plt.figure()
    cmap = plt.cm.get_cmap('plasma')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp, cmap=cmap, vmin=.96, vmax=1.04, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Amplitude', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_solution_example1.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('plasma')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp+.02, cmap=cmap, vmin=.96, vmax=1.04, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Amplitude', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_solution_example2.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('plasma')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp-.02, cmap=cmap, vmin=.96, vmax=1.04, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Amplitude', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_solution_example3.png'.format(output_path))
    plt.close()


def plot_phase_terms_from_csv():

    csv_path = '/Users/ruby/EoR/hex_array_sim_331_antenna_locs.csv'
    output_path = '/Users/ruby/EoR'

    csv_file = open(csv_path, 'r')
    csv_contents = csv_file.readlines()
    csv_file.close()
    antenna_xlocs = np.zeros(len(csv_contents)-1)
    antenna_ylocs = np.zeros(len(csv_contents)-1)
    for ind, line in enumerate(csv_contents[1:]):
        line_split = line.split(',')
        antenna_xlocs[ind] = float(line_split[1])
        antenna_ylocs[ind] = float(line_split[2])
    print(antenna_xlocs)
    random_amp = np.random.normal(loc=0., scale=4e-5, size=len(csv_contents)-1)

    plt.figure()
    cmap = plt.cm.get_cmap('viridis')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp, cmap=cmap, vmin=-1e-4, vmax=1e-4, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Phase', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_phase_solution_example1.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('viridis')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp+4e-5, cmap=cmap, vmin=-1e-4, vmax=1e-4, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Phase', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_phase_solution_example2.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('viridis')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp+5e-5, cmap=cmap, vmin=-1e-4, vmax=1e-4, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Phase', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_phase_solution_example2.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('viridis')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp+2e-5+1e-6*antenna_xlocs, cmap=cmap, vmin=-1e-4, vmax=1e-4, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Phase', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_phase_solution_example3.png'.format(output_path))
    plt.close()

    plt.figure()
    cmap = plt.cm.get_cmap('viridis')
    #sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1., color=cmap(norm(random_amp)), cmap=cmap)
    sc = plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=50., c=random_amp+.5e-6*antenna_xlocs-1e-6*antenna_ylocs, cmap=cmap, vmin=-1e-4, vmax=1e-4, edgecolors='black', linewidths=.5)
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(sc)
    cbar.ax.set_ylabel('Antenna Gain Phase', rotation=270, labelpad=15)  # label colorbar
    plt.savefig('{}/hex_array_cal_phase_solution_example4.png'.format(output_path))
    plt.close()




if __name__ == '__main__':
    plot_antenna_locs_from_csv()
