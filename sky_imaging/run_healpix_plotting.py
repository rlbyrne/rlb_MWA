#!/usr/bin/python

import numpy as np
import healpy as hp
import sys
import os
import matplotlib
#matplotlib.use('Agg')  # use this if you don't have display access
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import healpix_utils
import scipy
from scipy.interpolate import griddata
import plot_healpix_map


def combine_maps_Feb26():

    combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        nside=256, cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_combined_60obs_closest.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_combined_60obs_closest.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_combined_60obs_closest.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_combined_60obs_closest.fits'
    )

    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    combined_maps, weights_map = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        obs_list, nside=256,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_combined_60obs_ave.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_combined_60obs_ave.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_combined_60obs_ave.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_combined_60obs_ave.fits'
    )
    weights_map.write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Weights_combined_60obs.fits'
    )


def combine_maps_Feb27():

    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    variance_maps, averaged_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        obs_list,
        saved_averaged_maps=[
            '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_combined_60obs_ave.fits',
            '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_combined_60obs_ave.fits',
            '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_combined_60obs_ave.fits',
            '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_combined_60obs_ave.fits'
        ],
        obs_weights=None, nside=256,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        apply_radial_weighting=False
    )

    variance_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_combined_60obs_variance.fits'
    )
    variance_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_combined_60obs_variance.fits'
    )
    variance_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_combined_60obs_variance.fits'
    )
    variance_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_combined_60obs_variance.fits'
    )


def combine_maps_Feb27_with_filtering():

    map = healpix_utils.load_map('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/1130773144_uniform_Residual_I_HEALPix.fits')
    map.resample(128)
    map.filter_map(lmin=None, lmax=50, filter_width=5)
    plot_healpix_map.plot_filled_pixels(
        map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/filteringtest.png'
    )

def combine_maps_Mar1():

    """combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        nside=256, cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_closest.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_closest.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_closest.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_closest.fits'
    )

    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    combined_maps, weights_map = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        obs_list, nside=256,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_ave.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_ave.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_ave.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_ave.fits'
    )"""

    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_closest.png')

    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_ave.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesI_filtered_60obs_ave.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_ave.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_filtered_60obs_ave.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_ave.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesU_filtered_60obs_ave.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_ave.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesV_filtered_60obs_ave.png')


def plot_maps_Mar4():

    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    print len(obs_list)

    for obsid in obs_list:
        obs_struct = scipy.io.readsav(
            '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/metadata/{}_obs.sav'.format(obsid)
        )['obs']
        if np.min([(float(obs_struct['obsra'])-val)**2 for val in [0, 360, -360]]) < 15.**2. and (float(obs_struct['obsdec'])+27)**2. < 15.**2.:
            for cube in ['I', 'Q', 'U', 'V']:
                map = healpix_utils.load_map('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/{}_uniform_Residual_{}_HEALPix.fits'.format(
                    obsid, cube
                ))
                map.resample(512)
                plot_healpix_map.plot_filled_pixels(
                    map,
                    '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/eor0_plots/{}_Stokes{}_EoR0.png'.format(obsid, cube),
                    ra_range=[-15, 15], dec_range=[-42, -12], colorbar_range=[-.05,.05]
                )


if __name__ == '__main__':

    plot_maps_Mar4()

    #map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Weights_combined_60obs.fits')
    #plot_healpix_map.plot_filled_pixels(
    #    map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Weights_combined_60obs.png',
    #    colorbar_label='Number of Observations'
    #)
