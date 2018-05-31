#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import math
from pyuvdata import UVData
from pyuvdata import uvutils


def create_hex_array(side_length, antenna_spacing, save_uvfits=True,
                     plot_array=True):

    save_uvfits=False

    # get the number of antennas in an array with this side length
    antennas = ((side_length-1)**2 + (side_length-1))*3 + 1

    # initialize antenna locations with the center antenna
    antenna_xlocs = []
    antenna_ylocs = []

    add_extra_antenna = 1
    for row in range(0, side_length):
        yloc = row*antenna_spacing*(3.**.5)/2
        if row % 2 == 0:
            new_xlocs = [
                i*antenna_spacing for i in range(-side_length+row/2+1, side_length-row/2 + add_extra_antenna)
                ]
            antenna_xlocs.extend(new_xlocs)
            antenna_ylocs.extend([yloc]*len(new_xlocs))
            if row != 0:
                antenna_xlocs.extend(new_xlocs)
                antenna_ylocs.extend([-yloc]*len(new_xlocs))
        else:
            new_xlocs = [
                (i+.5)*antenna_spacing for i in range(int(-side_length+row/2.+.5), int(side_length-row/2.-.5))
                ]
            antenna_xlocs.extend(new_xlocs*2)
            antenna_ylocs.extend([yloc]*len(new_xlocs))
            antenna_ylocs.extend([-yloc]*len(new_xlocs))
        add_extra_antenna = 0
    antennas = len(antenna_xlocs)

    if plot_array:
        plt.figure()
        plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1.)
        # plt.xlim(ra_plot_range[0],ra_plot_range[1])
        # plt.ylim(dec_plot_range[0],dec_plot_range[1])
        plt.xlabel('East/West Location (m)')
        plt.ylabel('North/South Location (m)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig('/Users/rubybyrne/array_simulation/antenna_hex.png')
        plt.close()

    # find the array's radial distribution
    radial_dist = [
        (antenna_xlocs[i]**2+antenna_ylocs[i]**2)**.5 for i in range(len(antenna_xlocs))
        ]
    radial_hist, bin_edges = np.histogram(
        radial_dist, bins=int(side_length*.75)
        )
    bin_centers = [
        (bin_edges[i]+bin_edges[i+1])/2 for i in range(len(bin_edges)-1)
        ]
    plt.figure()
    plt.scatter(bin_centers, radial_hist)
    plt.savefig('/Users/rubybyrne/array_simulation/antenna_hex_dist.png')
    plt.close()

    if save_uvfits:
        create_uvfits(antennas, antenna_xlocs, antenna_ylocs, '/Users/rubybyrne/array_simulation/hex_array_sim.uvfits')

    print antennas

    return antennas, radial_hist, bin_centers


def create_random_array():

    save_uvfits = True
    plot_array = True

    antenna_size = 5.
    antennas, radial_hist, bin_centers = create_hex_array(7, 15.,
                                                          save_uvfits=False)
    print antennas
    radial_vals = np.arange(bin_centers[0], bin_centers[-1], .1)
    radial_distribution = np.interp(radial_vals, bin_centers, radial_hist)
    norm_factor = sum(radial_distribution)
    radial_distribution = [i/norm_factor for i in radial_distribution]

    antenna_xlocs = []
    antenna_ylocs = []
    for ant in range(antennas):

        random = np.random.uniform()
        dist_sum = 0.
        for i, dist_val in enumerate(radial_distribution):
            dist_sum += dist_val
            if random < dist_sum:
                break
        radius = radial_vals[i]

        check_overlap = False
        while not check_overlap:
            azimuth = np.random.uniform()*2*math.pi
            xloc = radius*math.cos(azimuth)
            yloc = radius*math.sin(azimuth)

            # Check that the antenna does not overlap others
            check_overlap = True
            for i in range(len(antenna_xlocs)):
                if abs(antenna_xlocs[i]-xloc) < antenna_size and abs(antenna_ylocs[i]-yloc) < antenna_size:
                    check_overlap = False
            if check_overlap:
                antenna_xlocs.append(xloc)
                antenna_ylocs.append(yloc)

    if plot_array:
        plt.figure()
        plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1.)
        # plt.xlim(ra_plot_range[0],ra_plot_range[1])
        # plt.ylim(dec_plot_range[0],dec_plot_range[1])
        plt.xlabel('East/West Location (m)')
        plt.ylabel('North/South Location (m)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig('/Users/rubybyrne/array_simulation/antenna_random1.png')
        plt.close()

    if save_uvfits:
        create_uvfits(antennas, antenna_xlocs, antenna_ylocs, '/Users/rubybyrne/array_simulation/random1_array_sim.uvfits')


def create_uvfits(antennas, antenna_xlocs, antenna_ylocs, save_filename):

    antenna_locs_ENU = np.zeros((antennas, 3))
    antenna_locs_ENU[:, 0] = antenna_xlocs
    antenna_locs_ENU[:, 1] = antenna_ylocs

    filename = '/Users/Shared/uvfits/4.1/1061316296.uvfits'
    UV = UVData()
    UV.read_uvfits(filename)
    phase_center_ra = UV.phase_center_ra
    phase_center_dec = UV.phase_center_dec
    phase_center_epoch = UV.phase_center_epoch
    UV.unphase_to_drift()  # unphase data
    for ant1 in range(antenna_locs_ENU.shape[0]):
        for ant2 in range(antenna_locs_ENU.shape[0]):
            baseline_inds = np.intersect1d(
                np.where(UV.ant_1_array == ant1)[0],
                np.where(UV.ant_2_array == ant2)[0]
                )
            UV.uvw_array[baseline_inds, :] = antenna_locs_ENU[ant2, :] - antenna_locs_ENU[ant1, :]
    antenna_locs_ECEF = uvutils.ECEF_from_ENU(
        antenna_locs_ENU.T, *UV.telescope_location_lat_lon_alt
        ).T
    UV.antenna_positions = antenna_locs_ECEF - UV.telescope_location
    UV.phase(phase_center_ra, phase_center_dec, phase_center_epoch)
    UV.write_uvfits(save_filename, spoof_nonessential=True)


if __name__ == '__main__':
    create_random_array()
