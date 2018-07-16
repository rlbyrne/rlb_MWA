#!/usr/bin/python

import numpy as np
import csv
from pyuvdata import UVData
from pyuvdata import uvutils


def write_out_antenna_locs():

    uvfits_path = '/Users/rubybyrne/array_simulation'
    save_path = '/Users/rubybyrne/array_simulation'

    obsids = ['hex_array_sim_10m','hex_array_sim_15m','hex_array_sim_20m',
        'random1_array_sim_10m','random1_array_sim_15m','random1_array_sim_20m',
        'random2_array_sim_10m','random2_array_sim_15m','random2_array_sim_20m',
        'random3_array_sim_10m','random3_array_sim_15m','random3_array_sim_20m'
    ]

    for obsid in obsids:

        UV = UVData()
        print 'Reading {}/{}.uvfits'.format(uvfits_path, obsid)
        UV.read_uvfits('{}/{}.uvfits'.format(uvfits_path, obsid))
        antenna_locs_ENU = uvutils.ENU_from_ECEF(
                (UV.antenna_positions + UV.telescope_location).T,
                *UV.telescope_location_lat_lon_alt
            ).T
        antenna_nums = UV.antenna_numbers

        print 'Saving outputs to {}/{}_antenna_locs.csv'.format(save_path, obsid)
        outfile = open('{}/{}_antenna_locs.csv'.format(save_path, obsid), 'w')
        outfile_writer = csv.writer(outfile)
        outfile_writer.writerow(['Antenna Number', 'E-W Location (m)',
                                 'N-S Location (m)', 'Altitude (m)'])
        for i in range(len(antenna_nums)):
            outfile_writer.writerow(
                [antenna_nums[i]]+list(antenna_locs_ENU[i, :])
            )
        outfile.close()

if __name__ == '__main__':
    write_out_antenna_locs()
