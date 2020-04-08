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


if __name__ == '__main__':

    plot_maps_Apr7()
