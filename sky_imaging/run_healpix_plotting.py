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


def plot_maps_Mar15():

    """combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_pol_leakage_correction_4pol_Mar2019',
        nside=256, cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_closest.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_closest.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_closest.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesV_60obs_closest.fits'
    )

    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_closest.png')
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesV_60obs_closest.fits')
    plot_healpix_map.plot_filled_pixels(map, '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesV_60obs_closest.png')
    """
    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_pol_leakage_correction_4pol_Mar2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    combined_maps, weights_map = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_pol_leakage_correction_4pol_Mar2019',
        obs_list, nside=256, apply_radial_weighting=True,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_ave.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_ave.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_ave.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesV_60obs_ave.fits'
    )
    weights_map.write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/Weights_60obs_ave.fits'
    )

    plot_healpix_map.plot_filled_pixels(combined_maps[0], '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesI_60obs_ave.png')
    plot_healpix_map.plot_filled_pixels(combined_maps[1], '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesQ_60obs_ave.png')
    plot_healpix_map.plot_filled_pixels(combined_maps[2], '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesU_60obs_ave.png')
    plot_healpix_map.plot_filled_pixels(combined_maps[4], '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/StokesV_60obs_ave.png')
    plot_healpix_map.plot_filled_pixels(weights_map, '/Users/rubybyrne/diffuse_survey_plotting_Mar2019/Weights_60obs_ave.png')


def plot_maps_Jul3():

    combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Jul2019', obs_list_file=None, nside=128, cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_combined_closest.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_combined_closest.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_combined_closest.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_combined_closest.fits'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[0], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_combined_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[1], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_combined_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[2], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_combined_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[3], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_combined_closest.png'
    )


def plot_maps_Jul4():

    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Jul2019',
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        apply_radial_weighting=True
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_ave.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_ave.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_ave.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_ave.fits'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[0], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_ave.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[1], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_ave.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[2], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_ave.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[3], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_ave.png'
    )


def plot_maps_Jul10():

    data_files = os.listdir('/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Jul2019/output_data/')
    data_files = [
        file for file in data_files
        if '_uniform_Residual_I_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]
    obs_list = [obs for obs in obs_list if obs not in ['1131713752','1131717112','1131719032','1131724672','1130776864','1131717232','1131731752','1131718912','1131726472','1131728032','1131728272',]]
    print len(obs_list)
    combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Jul2019',
        obs_list=obs_list,
        nside=256,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    combined_maps[0].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_closest.fits'
    )
    combined_maps[1].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_closest.fits'
    )
    combined_maps[2].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_closest.fits'
    )
    combined_maps[3].write_data_to_fits(
        '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_closest.fits'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[0], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesI_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[1], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesQ_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[2], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesU_closest.png'
    )
    plot_healpix_map.plot_filled_pixels(
        combined_maps[3], '/Users/rubybyrne/diffuse_survey_plotting_Jul2019/StokesV_closest.png'
    )


def plot_maps_Aug25():

    saveloc = '/Users/rubybyrne/diffuse_survey_plotting_Aug2019'
    combined_maps = healpix_utils.combine_maps_nearest_data(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
        nside=256,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V']
    )
    pol_names = ['I', 'Q', 'U', 'V']
    for map_ind in range(len(combined_maps)):
        combined_maps[map_ind].write_data_to_fits(
            '{}/Stokes{}_nearest_short_baselines.fits'.format(saveloc, pol_names[map_ind])
        )
    for map_ind in range(len(combined_maps)):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[map_ind], '{}/Stokes{}_nearest_short_baselines.png'.format(saveloc, pol_names[map_ind])
        )



def plot_maps_Dec3():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]
    print len(obs_list)
    sys.exit()
    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
        obs_list = obs_list,
        nside=512,
        cube_names=['Dirty_I', 'Dirty_Q', 'Dirty_U', 'Dirty_V'],
        apply_radial_weighting=True
    )
    #combined_maps = healpix_utils.combine_maps_nearest_data(
    #    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
    #    obs_list = obs_list,
    #    nside=512,
    #    cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
    #)
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Dec2019'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_dirty_averaged.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_dirty_averaged.png'.format(outdir, pol_name)
        )


def save_maps_Dec19():

    stokesI = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesI_averaged.fits'
    )
    stokesQ = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesQ_averaged.fits'
    )
    stokesU = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesU_averaged.fits'
    )
    stokesV = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesV_averaged.fits'
    )
    healpix_utils.write_data_to_standard_fits(
        [stokesI, stokesQ, stokesU, stokesV],
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/diffuse_averaged_maps.fits',
        history_str='53 observations processed with FHD and averaged'
    )

def plot_maps_Jan18():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]
    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Jan2020',
        obs_list = obs_list,
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True
    )
    #combined_maps = healpix_utils.combine_maps_nearest_data(
    #    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
    #    obs_list = obs_list,
    #    nside=512,
    #    cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
    #)
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Jan2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_averaged.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_averaged.png'.format(outdir, pol_name)
        )


def plot_maps_Jan20():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]
    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_no_pol_leakage_correction_Jan2020',
        obs_list = obs_list,
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True
    )
    #combined_maps = healpix_utils.combine_maps_nearest_data(
    #    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
    #    obs_list = obs_list,
    #    nside=512,
    #    cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
    #)
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Jan2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_no_pol_leakage_correction_averaged.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_no_pol_leakage_correction_averaged.png'.format(outdir, pol_name)
        )


def plot_maps_Feb03():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]
    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        obs_list = obs_list,
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True
    )
    #combined_maps = healpix_utils.combine_maps_nearest_data(
    #    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_baseline_cut_Aug2019',
    #    obs_list = obs_list,
    #    nside=512,
    #    cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
    #)
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Feb2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_averaged.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_averaged.png'.format(outdir, pol_name)
        )

def plot_maps_Mar20():

    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020',
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Mar2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_averaged.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_averaged.png'.format(outdir, pol_name)
        )

def plot_maps_Apr7():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624'
    ]

    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=False
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Apr2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_averaged_notaper.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_averaged_notaper.png'.format(outdir, pol_name)
        )


def test_rm_correction_Apr22():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        obs_lists = obs_list_1,
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=False,
        apply_rm_correction=True
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Apr2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_rm_corrected.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
        else:
            colorbar_range = [-2e3, 2e3]
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_rm_corrected.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )


def plot_variance_maps_May7():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624'
    ]

    averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=False,
        apply_rm_correction=True
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        averaged_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_average_map.fits'.format(outdir, pol_name)
        )
        variance_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_variance_map.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
            var_colorbar_range = [0, 1e4]
            snr_colorbar_range = [0, 2]
        else:
            colorbar_range = [-2e3, 2e3]
            var_colorbar_range = [0, 2e3]
            snr_colorbar_range = [0, 2]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )
        variance_maps[pol_ind].signal_arr = np.sqrt(variance_maps[pol_ind].signal_arr)
        plot_healpix_map.plot_filled_pixels(
            variance_maps[pol_ind],
            '{}/Stokes{}_stddev_map.png'.format(outdir, pol_name),
            colorbar_range=var_colorbar_range, colorbar_label='Standard Deviation (Jy/sr)'
        )
        plot_healpix_map.plot_filled_pixels(
            snr_maps[pol_ind],
            '{}/Stokes{}_snr_map.png'.format(outdir, pol_name),
            colorbar_range=snr_colorbar_range, colorbar_label='Signal Amplitude/Standard Dev.'
        )
    plot_healpix_map.plot_filled_pixels(
        weights_map,
        '{}/weights_map.png'.format(outdir)
    )


def plot_maps_May26():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        #'1131709792', stripes in V
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        #'1131740872', excess power
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        #'1131735952', # maybe excess power
        #'1131739552', excess power
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        #'1131713272', excess power
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    combined_maps, weight_maps = healpix_utils.average_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True,
        apply_rm_correction=True
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        combined_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_residual_averaged_more_obs.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
        else:
            colorbar_range = [-2e3, 2e3]
        plot_healpix_map.plot_filled_pixels(
            combined_maps[pol_ind],
            '{}/Stokes{}_residual_averaged_more_obs.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )


def plot_variance_maps_May26():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        #'1131709792', stripes in V
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        #'1131740872', excess power
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        #'1131735952', # maybe excess power
        #'1131739552', excess power
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        #'1131713272', excess power
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    print len(obs_list_1)+len(obs_list_2)

    averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True,
        apply_rm_correction=False
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        averaged_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_average_map_more_obs_no_rm_correction.fits'.format(outdir, pol_name)
        )
        variance_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_variance_map_more_obs_no_rm_correction.fits'.format(outdir, pol_name)
        )
    weights_map.write_data_to_fits(
        '{}/weights_map_more_obs.fits'.format(outdir)
    )
    nsamples_map.write_data_to_fits(
        '{}/nsamples_map_more_obs.fits'.format(outdir)
    )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
            var_colorbar_range = [0, 1e4]
            snr_colorbar_range = [0, 2]
        else:
            colorbar_range = [-2e3, 2e3]
            var_colorbar_range = [0, 2e3]
            snr_colorbar_range = [0, 2]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map_more_obs_no_rm_correction.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )
        variance_maps[pol_ind].signal_arr = np.sqrt(variance_maps[pol_ind].signal_arr)
        plot_healpix_map.plot_filled_pixels(
            variance_maps[pol_ind],
            '{}/Stokes{}_stddev_map_more_obs_no_rm_correction.png'.format(outdir, pol_name),
            colorbar_range=var_colorbar_range, colorbar_label='Standard Deviation (Jy/sr)'
        )
        plot_healpix_map.plot_filled_pixels(
            snr_maps[pol_ind],
            '{}/Stokes{}_snr_map_more_obs_no_rm_correction.png'.format(outdir, pol_name),
            colorbar_range=snr_colorbar_range, colorbar_label='Signal Amplitude/Standard Dev.'
        )
    plot_healpix_map.plot_filled_pixels(
        weights_map,
        '{}/weights_map_more_obs.png'.format(outdir), colorbar_label='Weights'
    )
    plot_healpix_map.plot_filled_pixels(
        nsamples_map,
        '{}/nsamples_map_more_obs.png'.format(outdir),
        colorbar_label='Number of Observations'
    )


def plot_maps_no_rm_correction_May27():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        #'1131709792', stripes in V
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        #'1131740872', excess power
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        #'1131735952', # maybe excess power
        #'1131739552', excess power
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        #'1131713272', excess power
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    print len(obs_list_1)+len(obs_list_2)

    averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True,
        apply_rm_correction=False
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        averaged_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_average_map_no_rm_correction_more_obs.fits'.format(outdir, pol_name)
        )

    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
            var_colorbar_range = [0, 1e4]
            snr_colorbar_range = [0, 2]
        else:
            colorbar_range = [-2e3, 2e3]
            var_colorbar_range = [0, 2e3]
            snr_colorbar_range = [0, 2]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map_no_rm_correction_more_obs.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )


def plot_projected_maps_Jun1():

    for pol in ['I', 'Q', 'U', 'V']:
        if pol == 'I':
            colorbar_range = [-2e4, 2e4]
        else:
            colorbar_range = [-5e3, 5e3]
        map = healpix_utils.load_map(
            '/Users/rubybyrne/diffuse_survey_plotting_May2020/Stokes{}_average_map_more_obs.fits'.format(pol)
        )
        plot_healpix_map.plot_projection(
            map,
            title='Stokes {}'.format(pol),
            save_filename='/Users/rubybyrne/diffuse_survey_plotting_May2020/Stokes{}_average_map_more_obs_proj.png'.format(pol),
            colorbar_range=colorbar_range
        )

def write_out_images_Jun3():

    maps = []
    for pol in ['I', 'Q', 'U', 'V']:
        map = healpix_utils.load_map(
            '/Users/rubybyrne/diffuse_survey_plotting_May2020/Stokes{}_average_map_more_obs.fits'.format(pol)
        )
        maps.append(map)
    healpix_utils.write_data_to_standard_fits(
        maps,
        '/Users/rubybyrne/diffuse_survey_plotting_May2020/polarized_diffuse_map.fits',
        history_str='produced by Ruby Byrne, U. Washington, June 2020'
    )


def plot_variance_maps_Aug4():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        #'1131709792', stripes in V
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        #'1131740872', excess power
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        #'1131735952', # maybe excess power
        #'1131739552', excess power
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        #'1131713272', excess power
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
        obs_lists = [obs_list_1, obs_list_2],
        nside=128,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='weighted',
        apply_radial_weighting=True,
        apply_rm_correction=True,
        rm_file='/Users/rubybyrne/rm_empirical_calculation/Jul2020_align_with_avg/diffuse_survey_rm_empirical_in_eor0_Aug2020.csv'
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        averaged_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(outdir, pol_name)
        )
        variance_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_variance_map_empirical_rm_in_eor0.fits'.format(outdir, pol_name)
        )
    weights_map.write_data_to_fits(
        '{}/weights_map.fits'.format(outdir)
    )
    nsamples_map.write_data_to_fits(
        '{}/nsamples_map.fits'.format(outdir)
    )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
            var_colorbar_range = [0, 1e4]
            snr_colorbar_range = [0, 2]
        else:
            colorbar_range = [-2e3, 2e3]
            var_colorbar_range = [0, 2e3]
            snr_colorbar_range = [0, 2]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range
        )
        variance_maps[pol_ind].signal_arr = np.sqrt(variance_maps[pol_ind].signal_arr)
        plot_healpix_map.plot_filled_pixels(
            variance_maps[pol_ind],
            '{}/Stokes{}_stddev_map_empirical_rm_in_eor0.png'.format(outdir, pol_name),
            colorbar_range=var_colorbar_range, colorbar_label='Standard Deviation (Jy/sr)'
        )
    plot_healpix_map.plot_filled_pixels(
        weights_map,
        '{}/weights_map.png'.format(outdir), colorbar_label='Weights'
    )
    plot_healpix_map.plot_filled_pixels(
        nsamples_map,
        '{}/nsamples_map.png'.format(outdir),
        colorbar_label='Number of Observations'
    )


def undo_rm_correction_Aug4():

    obsid = '1131454296'
    rm_file='/Users/rubybyrne/rm_empirical_calculation/Jul2020_align_with_avg/diffuse_survey_rm_empirical_in_eor0_Aug2020.csv'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    maps = [None]
    for pol_name in ['Q', 'U']:
        map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(outdir, pol_name)
        )
        maps.append(map)
    maps.append(None)
    print np.mean(np.sqrt(maps[1].signal_arr**2+maps[2].signal_arr**2))
    maps = healpix_utils.undo_rm_correction(
        obsid, maps, rm_file=rm_file
    )
    print np.mean(np.sqrt(maps[1].signal_arr**2+maps[2].signal_arr**2))
    maps[1].write_data_to_fits(
        '{}/StokesQ_average_map_{}_rm_undone.fits'.format(outdir, obsid)
    )
    maps[2].write_data_to_fits(
        '{}/StokesU_average_map_{}_rm_undone.fits'.format(outdir, obsid)
    )
    colorbar_range = [-2e3, 2e3]
    plot_healpix_map.plot_filled_pixels(
        maps[1],
        '{}/StokesQ_average_map_{}_rm_undone.png'.format(outdir, obsid),
        colorbar_range=colorbar_range
    )
    plot_healpix_map.plot_filled_pixels(
        maps[2],
        '{}/StokesU_average_map_{}_rm_undone.png'.format(outdir, obsid),
        colorbar_range=colorbar_range
    )


def overplot_points_Aug11():

    obs_list_1 = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        #'1131454176', high power and systematics in Stokes V
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        #'1131456216',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        #'1131475896',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        #'1131477816',
        '1131459696',
        '1130780464',
        '1131726352',
        #'1131715312',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        #'1131724672',
        '1131544424'
    ]

    obs_list_2 = ['1131542504',
        #'1131717112',
        '1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        #'1131728032',
        '1130787544',
        '1130776744',
        #'1131726472',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1130776624',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        #'1131728272',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        #'1131713752',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        #'1131720712',
        #'1131711832',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        #'1131711712',
        '1131472536',
        #'1131729832',
        '1130773024',
        #'1131720952',
        #'1131718912',
        #'1131719032',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        #'1131709792', stripes in V
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        #'1131740872', excess power
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        #'1131735952', # maybe excess power
        #'1131739552', excess power
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        #'1131713272', excess power
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    outpath = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'

    ras = np.zeros(len(obs_list_1 + obs_list_2))
    decs = np.zeros(len(obs_list_1 + obs_list_2))
    alts = np.zeros(len(obs_list_1 + obs_list_2))
    for obsind, obsid in enumerate(obs_list_1 + obs_list_2):
        if obsid in obs_list_1:
            path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020'
        elif obsid in obs_list_2:
            path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'
        obs_struct = scipy.io.readsav(
            '{}/metadata/{}_obs.sav'.format(path, obsid)
        )['obs']
        ras[obsind] = float(obs_struct['obsra'])
        decs[obsind] = float(obs_struct['obsdec'])
        alts[obsind] = float(obs_struct['obsalt'])
    ras[np.where(ras > 270.)] -= 360.

    maps = []
    for pol_name in ['I', 'Q', 'U']:
        map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(outpath, pol_name)
        )
        maps.append(map)

    # Plot RMs
    rm_file = '/Users/rubybyrne/rm_empirical_calculation/Jul2020_align_with_avg/diffuse_survey_rm_empirical_in_eor0_Aug2020.csv'
    rm_data = np.genfromtxt(
        rm_file, delimiter=',', dtype=None, names=True, encoding=None
    )
    rms = np.array([
        rm_data['RM'][np.where(rm_data['ObsID'] == int(obsid))][0] for obsid in obs_list_1+obs_list_2
    ])
    vmin = -4
    vmax = 0

    use_colormap = 'plasma'
    cm = plt.cm.get_cmap(use_colormap)
    plt.figure(figsize=[15,5])
    plt.scatter(
        ras, decs, c=rms, vmin=vmin, vmax=vmax, s=30, cmap=cm,
        edgecolor='black', linewidth=.5
    )
    plt.xlim([130,-40])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('RM (rad/m^2)', rotation=270, labelpad=15)
    plt.title('Rotation Measures')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.gca().set_aspect('equal')
    print 'Saving scatter plot to {}/scatterplot_rms.png'.format(outpath)
    plt.savefig('{}/scatterplot_rms.png'.format(outpath), format='png', dpi=300)
    plt.close()

    for pol_ind, pol_name in enumerate(['I', 'Q', 'U']):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
        else:
            colorbar_range = [-2e3, 2e3]
        plot_healpix_map.plot_filled_pixels(
            maps[pol_ind],
            '{}/Stokes{}_with_rms.png'.format(outpath, pol_name),
            colorbar_range=colorbar_range,
            overplot_points=True, point_ras=ras, point_decs=decs,
            point_values=rms,
            overplot_points_vmin=vmin, overplot_points_vmax=vmax,
            overplot_points_colormap=use_colormap
        )

    # Plot Altitudess
    vmin = np.min(alts)
    vmax = 90

    use_colormap = 'plasma'
    cm = plt.cm.get_cmap(use_colormap)
    plt.figure(figsize=[15,5])
    plt.scatter(
        ras, decs, c=alts, vmin=vmin, vmax=vmax, s=30, cmap=cm,
        edgecolor='black', linewidth=.5
    )
    plt.xlim([130,-40])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Altitudes (Deg)', rotation=270, labelpad=15)
    plt.title('Altitudes')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.gca().set_aspect('equal')
    print 'Saving scatter plot to {}/scatterplot_alts.png'.format(outpath)
    plt.savefig('{}/scatterplot_alts.png'.format(outpath), format='png', dpi=300)
    plt.close()

    for pol_ind, pol_name in enumerate(['I', 'Q', 'U']):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
        else:
            colorbar_range = [-2e3, 2e3]
        plot_healpix_map.plot_filled_pixels(
            maps[pol_ind],
            '{}/Stokes{}_with_alts.png'.format(outpath, pol_name),
            colorbar_range=colorbar_range,
            overplot_points=True, point_ras=ras, point_decs=decs,
            point_values=alts,
            overplot_points_vmin=vmin, overplot_points_vmax=vmax,
            overplot_points_colormap=use_colormap
        )

    # Plot Obsids
    obsids_int = [int(obsid) for obsid in obs_list_1 + obs_list_2]
    vmin = np.min(obsids_int)
    vmax = np.max(obsids_int)

    use_colormap = 'plasma'
    cm = plt.cm.get_cmap(use_colormap)
    plt.figure(figsize=[15,5])
    plt.scatter(
        ras, decs, c=obsids_int, vmin=vmin, vmax=vmax, s=30, cmap=cm,
        edgecolor='black', linewidth=.5
    )
    plt.xlim([130,-40])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Obsids', rotation=270, labelpad=15)
    plt.title('Obsids')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.gca().set_aspect('equal')
    print 'Saving scatter plot to {}/scatterplot_obsids.png'.format(outpath)
    plt.savefig('{}/scatterplot_obsids.png'.format(outpath), format='png', dpi=300)
    plt.close()

    for pol_ind, pol_name in enumerate(['I', 'Q', 'U']):
        if pol_name == 'I':
            colorbar_range = [-1e4, 1e4]
        else:
            colorbar_range = [-2e3, 2e3]
        plot_healpix_map.plot_filled_pixels(
            maps[pol_ind],
            '{}/Stokes{}_with_obsids.png'.format(outpath, pol_name),
            colorbar_range=colorbar_range,
            overplot_points=True, point_ras=ras, point_decs=decs,
            point_values=obsids_int,
            overplot_points_vmin=vmin, overplot_points_vmax=vmax,
            overplot_points_colormap=use_colormap
        )


def replot_maps_Sep23():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        if pol_ind == 0:
            colorbar_range = [-5e4, 5e4]
        else:
            colorbar_range = [-5e3, 5e3]
        map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        plot_healpix_map.plot_filled_pixels(
            map,
            '{}/Stokes{}_average_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, big=True
        )


def save_maps_Sep25():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020/polarized_diffuse_map.fits'
    pols = ['I', 'Q', 'U', 'V']
    maps = []
    for pol_ind, pol_name in enumerate(pols):
        new_map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        maps.append(new_map)
    healpix_utils.write_data_to_standard_fits(maps, save_filename)


def plot_maps_with_contours_Feb4():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        if pol_ind == 0:
            colorbar_range = [-5e4, 5e4]
        else:
            colorbar_range = [-5e3, 5e3]
        map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        plot_healpix_map.plot_filled_pixels(
            map,
            '{}/Stokes{}_average_map_with_annotations.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, overplot_mwa_beam_contours=True,
            overplot_hera_band=True, overplot_bright_sources=True
        )


def plot_variance_maps_Feb11():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        if pol_ind == 3:
            colorbar_range = [0, 5e5]
        else:
            colorbar_range = [0, 1e7]
        map = healpix_utils.load_map(
            '{}/Stokes{}_variance_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        plot_healpix_map.plot_filled_pixels(
            map,
            '{}/Stokes{}_variance_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, colorbar_label='Surface Brightness Variance (Jy/sr)${^2}$'
        )

    map = healpix_utils.load_map(
        '{}/nsamples_map.fits'.format(sourcedir)
    )
    plot_healpix_map.plot_filled_pixels(
        map,
        '{}/nsamples.png'.format(outdir),
        colorbar_range=[0, 25], colorbar_label='Number of Observations'
    )


def plot_stddev_maps_Mar8():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        if pol_ind == 3:
            colorbar_range = [0, 1e3]
        else:
            colorbar_range = [0, 5e3]
        map = healpix_utils.load_map(
            '{}/Stokes{}_variance_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        map.signal_arr = np.sqrt(map.signal_arr)
        plot_healpix_map.plot_filled_pixels(
            map,
            '{}/Stokes{}_stddev_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, colorbar_label='Surface Brightness Std. Dev. (Jy/sr)'
        )

def explore_map_Mar16():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    map = healpix_utils.load_map(
        '{}/StokesI_average_map_empirical_rm_in_eor0.fits'.format(sourcedir)
    )
    print(np.shape(map.signal_arr))
    print(map.nside)


def rm_correct_maps_Mar18():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/rm_corrected/fits_files'
    rm_file = '/Users/rubybyrne/mwa_2013_rm.csv'
    mapQ = healpix_utils.load_map(
        '{}/StokesQ_average_map_empirical_rm_in_eor0.fits'.format(sourcedir)
    )
    mapU = healpix_utils.load_map(
        '{}/StokesU_average_map_empirical_rm_in_eor0.fits'.format(sourcedir)
    )
    maps = ['', mapQ, mapU, '']

    rm_data = np.genfromtxt(
        rm_file, delimiter=',', dtype=None, names=True, encoding=None
    )
    obsids = rm_data['ObsID']
    for obsid in obsids:
        rotated_maps = healpix_utils.undo_rm_correction(obsid, maps, rm_file=rm_file)
        rotated_maps[1].write_data_to_fits(
            '{}/{}_RMcorrected_StokesQ.fits'.format(outdir, obsid)
        )
        rotated_maps[2].write_data_to_fits(
            '{}/{}_RMcorrected_StokesU.fits'.format(outdir, obsid)
        )


def replot_weights_Apr21():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020'
    map = healpix_utils.load_map(
        '{}/weights_map.fits'.format(sourcedir)
    )
    plot_healpix_map.plot_filled_pixels(
        map,
        '{}/weights_map.png'.format(outdir),
        colorbar_range=[.001, 10], big=False, colorbar_label='Weight'
    )


def make_fits_files_Apr22():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020/diffuse_map.fits'
    pols = ['I', 'Q', 'U', 'V']
    maps = []
    for pol_ind, pol_name in enumerate(pols):
        new_map = healpix_utils.load_map(
            '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        maps.append(new_map)
    healpix_utils.write_data_to_standard_fits(maps, save_filename)

    save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020/standard_deviation.fits'
    maps = []
    for pol_ind, pol_name in enumerate(pols):
        new_map = healpix_utils.load_map(
            '{}/Stokes{}_variance_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
        )
        new_map.signal_arr = np.sqrt(new_map.signal_arr)
        maps.append(new_map)
    healpix_utils.write_data_to_standard_fits(maps, save_filename, history_str='map standard deviations')


def combine_maps_with_variance_calculation_Jun14():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        '1131459696',
        '1130780464',
        '1131726352',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        '1131544424',
        '1131542504',
        #'1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        '1130787544',
        '1130776744',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        '1131472536',
        '1130773024',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map = healpix_utils.calculate_variance_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Jun2021',
        obs_lists = obs_list,
        nside=512,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        weighting='optimal',
        apply_radial_weighting=True,
        apply_rm_correction=True,
        rm_file='/Users/rubybyrne/diffuse_survey_rm_empirical_Jul2021.csv'
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Jun2021'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        variance_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_variance_map.fits'.format(outdir, pol_name)
        )
    weights_map.write_data_to_fits(
        '{}/weights_map.fits'.format(outdir)
    )
    nsamples_map.write_data_to_fits(
        '{}/nsamples_map.fits'.format(outdir)
    )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-5e4, 5e4]
        else:
            colorbar_range = [-5e3, 5e3]
        if pol_name == 'V':
            var_colorbar_range = [0, 1e3]
        else:
            var_colorbar_range = [0, 5e3]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, big=True
        )
        variance_maps[pol_ind].signal_arr = np.sqrt(variance_maps[pol_ind].signal_arr)
        plot_healpix_map.plot_filled_pixels(
            variance_maps[pol_ind],
            '{}/Stokes{}_stddev_map.png'.format(outdir, pol_name),
            colorbar_range=var_colorbar_range, colorbar_label='Standard Deviation (Jy/sr)', big=True
        )
    plot_healpix_map.plot_filled_pixels(
        weights_map,
        '{}/weights_map.png'.format(outdir), colorbar_label='Weights'
    )
    plot_healpix_map.plot_filled_pixels(
        nsamples_map,
        '{}/nsamples_map.png'.format(outdir),
        colorbar_label='Number of Observations'
    )


def combine_maps_Jun14():

    obs_list = [
        '1131551744',
        '1130783824',
        '1131562544',
        '1131709912',
        '1130776864',
        '1131461496',
        '1130782264',
        '1131715432',
        '1131733552',
        '1131542624',
        '1130773144',
        '1131461376',
        '1131557144',
        '1131454296',
        '1131731752',
        '1130778664',
        '1131470496',
        '1131559064',
        '1131717232',
        '1131463536',
        '1130773264',
        '1131463416',
        '1131717352',
        '1131713632',
        '1131478056',
        '1131468936',
        '1131468696',
        '1131535424',
        '1131463296',
        '1131465216',
        '1131710032',
        '1130776624',
        '1131456096',
        '1131540824',
        '1131711952',
        '1131459576',
        '1131477936',
        '1131733672',
        '1131564464',
        '1130787784',
        '1131461616',
        '1131558944',
        '1131470616',
        '1131549944',
        '1131553544',
        '1131459696',
        '1130780464',
        '1131726352',
        '1131470736',
        '1131548024',
        '1131710152',
        '1130785864',
        '1131544424',
        '1131542504',
        #'1131733432',
        '1131735232',
        '1131553664',
        '1131724432',
        '1131542744',
        '1131455976',
        '1131719152',
        '1131454416',
        '1130787544',
        '1130776744',
        '1130780224',
        '1131551624',
        '1131722632',
        '1131547904',
        '1131562664',
        '1131550064',
        '1131537104',
        '1131555224',
        '1131467136',
        '1131539024',
        '1131555344',
        '1131546104',
        '1131548144',
        '1131472416',
        '1131558824',
        '1131544304',
        '1130789584',
        '1131476136',
        '1130789344',
        '1131722872',
        '1130785744',
        '1131730072',
        '1131459816',
        '1131564584',
        '1131457776',
        '1131724552',
        '1130787664',
        '1130778424',
        '1131728152',
        '1131722752',
        '1131538904',
        '1131544544',
        '1130778544',
        '1131467016',
        '1131546344',
        '1130789464',
        '1131713512',
        '1131546224',
        '1131474336',
        '1130782144',
        '1131735472',
        '1130775064',
        '1130774824',
        '1131720832',
        '1130774944',
        '1131557264',
        '1130783944',
        '1131472296',
        '1131465096',
        '1131457896',
        '1131555464',
        '1131562424',
        '1131551864',
        '1131540704',
        '1130780344',
        '1131731632',
        '1131468816',
        '1131472536',
        '1130773024',
        '1131474096',
        '1131465336',
        '1131715552',
        '1131458016',
        '1131540944',
        '1131557024',
        '1131731872',
        '1131553424',
        '1131560864',
        '1130784064',
        '1131466896',
        '1130782024',
        '1131560624',
        '1131474216',
        '1131564344',
        '1131729952',
        '1131560744',
        '1130785624',
        '1131709432',
        '1131536624',
        '1131536384',
        '1131711112',
        '1131709192',
        '1131710992',
        '1131453456',
        '1131565304',
        '1131478776',
        '1131566504',
        '1131565184',
        '1131566624',
        '1131566744',
        '1131565064',
        '1131567944',
        '1131478656',
        '1131568544',
        '1131739432',
        '1130788504',
        '1130788264',
        '1131740752',
        '1131455736',
        '1131710392',
        '1131708952',
        '1131457176',
        '1131716512',
        '1131458976',
        '1131712192',
        '1131453936',
        '1131457536',
        '1131537704',
        '1131543584'
    ]

    averaged_maps, weights_maps = healpix_utils.average_healpix_maps(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Jun2021',
        obs_lists = obs_list,
        nside=512,
        cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
        apply_radial_weighting=True,
        weighting='optimal',
        apply_rm_correction=True,
        rm_file='/Users/rubybyrne/diffuse_survey_rm_empirical_Jul2021.csv'
    )
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Jun2021'
    pols = ['I', 'Q', 'U', 'V']
    for pol_ind, pol_name in enumerate(pols):
        averaged_maps[pol_ind].write_data_to_fits(
            '{}/Stokes{}_average_map.fits'.format(outdir, pol_name)
        )
    for pol_ind, pol_name in enumerate(pols):
        if pol_name == 'I':
            colorbar_range = [-5e4, 5e4]
        else:
            colorbar_range = [-5e3, 5e3]
        plot_healpix_map.plot_filled_pixels(
            averaged_maps[pol_ind],
            '{}/Stokes{}_average_map.png'.format(outdir, pol_name),
            colorbar_range=colorbar_range, big=True
        )


def rm_correct_maps_Jul6():

    sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Jun2021'
    outdir = '/Users/rubybyrne/diffuse_survey_plotting_Jun2021/rm_adjusted'
    rm_file = '/Users/rubybyrne/diffuse_survey_rm_empirical_Jul2021.csv'
    mapQ = healpix_utils.load_map(
        '{}/StokesQ_average_map.fits'.format(sourcedir)
    )
    mapU = healpix_utils.load_map(
        '{}/StokesU_average_map.fits'.format(sourcedir)
    )
    maps = ['', mapQ, mapU, '']

    rm_data = np.genfromtxt(
        rm_file, delimiter=',', dtype=None, names=True, encoding=None
    )
    obsid = '1131454296'
    rotated_maps = healpix_utils.undo_rm_correction(obsid, maps, rm_file=rm_file)
    rotated_maps[1].write_data_to_fits(
        '{}/{}_RMcorrected_StokesQ.fits'.format(outdir, obsid)
    )
    rotated_maps[2].write_data_to_fits(
        '{}/{}_RMcorrected_StokesU.fits'.format(outdir, obsid)
    )


if __name__ == '__main__':

    combine_maps_with_variance_calculation_Jun14()
