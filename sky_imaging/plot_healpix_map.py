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
from scipy.interpolate import griddata


def plot_filled_pixels(
    map, save_filename, title='', ra_range=[], dec_range=[],
    colorbar_range=[], log=False, colorbar_label='Flux Density (Jy/sr)',
    ra_cut=270.
):

    if map.coords == '':
        print 'WARNING: No map coordinate scheme supplied.'
        print 'Assuming equitorial coordinates.'
        coords = 'equitorial'
    else:
        coords = map.coords

    if len(ra_range) == 2 or len(dec_range) == 2:
        map.get_ra_dec(ra_cut=ra_cut)
        if len(ra_range) != 2:
            ra_range = [np.amin(map.ra_arr), np.amax(map.ra_arr)]
        else:
            ra_range = [np.amin(ra_range), np.amax(ra_range)]
        if len(dec_range) != 2:
            dec_range = [np.amin(map.dec_arr), np.amax(map.dec_arr)]
        use_indices = np.arange(len(map.signal_arr))[
            (map.ra_arr>ra_range[0]) & (map.ra_arr<ra_range[1])
            & (map.dec_arr>dec_range[0]) & (map.dec_arr<dec_range[1])
        ]
        use_map = healpix_utils.HealpixMap(
            map.signal_arr[use_indices], map.pix_arr[use_indices], map.nside,
            nest=map.nest, coords=coords
        )
    else:
        use_map = map

    use_map.get_pixel_corners(ra_cut=ra_cut)
    patches = []
    for ind in range(len(use_map.signal_arr)):
        polygon = Polygon(
            zip(use_map.pix_corner_ras_arr[ind],
                use_map.pix_corner_decs_arr[ind])
            )
        patches.append(polygon)
    colors = use_map.signal_arr

    # Establish axis ranges
    if len(ra_range) != 2:
        all_ras = [
            use_map.pix_corner_ras_arr[ind][poly_ind]
            for ind in range(len(use_map.pix_corner_ras_arr))
            for poly_ind in range(4)
        ]
        ra_range = [min(all_ras), max(all_ras)]
    if len(dec_range) != 2:
        all_decs = [
            use_map.pix_corner_decs_arr[ind][poly_ind]
            for ind in range(len(use_map.pix_corner_decs_arr))
            for poly_ind in range(4)
        ]
        dec_range = [min(all_decs), max(all_decs)]

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.05)
    collection.set_array(np.array(colors))  # set the data colors
    collection.set_edgecolor('face')  # make the face and edge colors match
    if log:  # set the color bar to a log scale
        collection.set_norm(LogNorm())
    if len(colorbar_range) == 2:  # set the colorbar min and max
        collection.set_clim(vmin=colorbar_range[0], vmax=colorbar_range[1])
    else:
        signal_mean = np.mean(colors)
        signal_std = np.std(colors)
        collection.set_clim(
            vmin=max([min(colors), signal_mean-5*signal_std]),
            vmax=min([max(colors), signal_mean+5*signal_std])
        )

    fig, ax = plt.subplots(figsize=(10, 8), dpi=500)
    ax.add_collection(collection)  # plot data

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([ra_range[1], ra_range[0], dec_range[0], dec_range[1]])
    plt.title(title)
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    # label colorbar
    cbar.ax.set_ylabel(colorbar_label, rotation=270, labelpad=15)
    print 'Saving plot to {}'.format(save_filename)
    plt.savefig(save_filename, format='png', dpi=300)


def plot_grid_interp(
    map, save_filename, resolution=.1, title='', ra_range=[], dec_range=[],
    colorbar_range=[None, None], log=False,
    colorbar_label='Flux Density (Jy/sr)'
):
    # resolution is in degrees

    map.get_ra_dec()
    for point in data:
        point.get_ra_dec(nside, nest, coords=coords)
    if len(ra_range) != 2:
        ra_range = [min(map.ra_arr), max(map.ra_arr)]
    if len(dec_range) != 2:
        dec_range = [min(map.dec_arr), max(map.dec_arr)]

    grid_dec, grid_ra = np.mgrid[
        dec_range[0]:dec_range[1]:resolution,
        ra_range[0]:ra_range[1]:resolution
        ]
    gridded_signal = griddata(
        (map.ra_arr, map.dec_arr), map.signal_arr, (grid_ra, grid_dec),
        method='linear'
        )

    fig, ax = plt.subplots(figsize=(21, 8), dpi=500)
    plt.imshow(
        gridded_signal, origin='lower', interpolation='none',
        extent=[ra_range[0], ra_range[1], dec_range[0], dec_range[1]],
        cmap='Greys_r', vmin=colorbar_range[0], vmax=colorbar_range[1]
        )
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([ra_range[1], ra_range[0], dec_range[0], dec_range[1]])
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar(extend='max')
    cbar.ax.set_ylabel(colorbar_label, rotation=270)  # label colorbar
    print 'Saving plot to {}'.format(save_filename)
    plt.savefig(save_filename, format='png', dpi=500)


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
    map.resample(8)
    map.filter_map(lmin=None, lmax=10, filter_width=2)
    plot_filled_pixels(
        map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/filteringtest.png', ra_range = [0,55], dec_range=[-50,0]
    )


if __name__ == '__main__':

    combine_maps_Feb27_with_filtering()

    #map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Weights_combined_60obs.fits')
    #plot_filled_pixels(
    #    map, '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Weights_combined_60obs.png',
    #    colorbar_label='Number of Observations'
    #)
