#!/usr/bin/python

import numpy as np
import os
import sys
sys.path.append('/Users/rubybyrne/rlb_MWA/sky_imaging/')
import healpix_utils
import qa_utils
import scipy.io
import healpy as hp
import matplotlib.pyplot as plt


def get_radial_variance(
    fhd_run_path, average_maps_paths, weights_map_path, savepath,
    cube_names=['Residual_I']
):

    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    average_maps = []
    for map_ind in range(len(cube_names)):
        map = healpix_utils.load_map(average_maps_paths[map_ind])
        map.explicit_to_implicit_ordering()
        average_maps.append(map)
    weights = healpix_utils.load_map(weights_map_path)
    if weights.nside != average_maps[0].nside:
        weights.resample(average_maps[0].nside)
    weights.explicit_to_implicit_ordering()

    data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]

    bin_edges = np.arange(0., 15., .05)
    count = np.zeros(len(bin_edges)-1, dtype=long)
    dev_from_mean = np.zeros((len(bin_edges)-1, len(cube_names)), dtype=float)

    for obs_ind, obsid in enumerate(obs_list):
        print 'Loading observation {} of {}'.format(obs_ind, len(obs_list))
        obs_struct = scipy.io.readsav(
            '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
        )['obs']
        obs_vec = hp.pixelfunc.ang2vec(
            float(obs_struct['obsra']), float(obs_struct['obsdec']),
            lonlat=True
        )
        maps = []
        for cube_ind, cube in enumerate(cube_names):
            map = healpix_utils.load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
                fhd_run_path, obsid, cube
            ))
            map.resample(average_maps[cube_ind].nside)
            if map.nest != average_maps[cube_ind].nest:
                print 'Map nesting conventions do not match. Converting.'
                if average_maps[cube_ind].nside:
                    map.reorder_ring_to_nest()
                else:
                    map.reorder_nest_to_ring()
            if map.coords != average_maps[cube_ind].coords:
                print 'ERROR: Map coordinates do not match. Exiting.'
                sys.exit(1)
            maps.append(map)

        for pix in maps[0].pix_arr:  # Use the first cube for the pixel list
            pix_vec = hp.pix2vec(maps[0].nside, pix, nest=maps[0].nest)
            rad_dist = hp.rotator.angdist(pix_vec, obs_vec)*180./np.pi
            bin_index = 0
            while rad_dist > bin_edges[bin_index+1] and bin_index < len(bin_edges)-1:
                bin_index += 1
            if bin_index < len(bin_edges)-1 and weights.signal_arr[pix] > 1:
                count[bin_index] += 1
                for cube_ind in range(len(cube_names)):
                    dev_from_mean[bin_index, cube_ind] += (
                        maps[cube_ind].signal_arr[
                            np.where(maps[cube_ind].pix_arr == pix)
                        ][0] - average_maps[cube_ind].signal_arr[pix]
                    )**2
    average_dev = dev_from_mean/count[:, None]
    average_dev[np.where(count == 0), :] = 0.
    xvals = np.array([
        (bin_edges[ind+1]+bin_edges[ind])/2. for ind in range(len(bin_edges)-1)
    ])

    plt.figure()
    colors = ['black', 'red', 'blue']
    for cube_ind, cube in enumerate(cube_names):
        plt.scatter(
            xvals, average_dev[:, cube_ind], marker='o', s=1.,
            color=colors[cube_ind], label=cube
        )
    plt.xlabel('Distance from Observation Center')
    plt.ylabel('Average of the Deviation from the Mean Squared')
    plt.ylim(0,.0001)
    plt.legend(loc='upper left')
    print 'Saving plot to {}'.format(savepath)
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()

if __name__=='__main__':
    get_radial_variance(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_pol_leakage_correction_4pol_Mar2019',
        ['/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_ave.fits',
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_ave.fits',
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_ave.fits'],
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/Weights_60obs_ave.fits',
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/radial_variance.png',
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U']
    )
