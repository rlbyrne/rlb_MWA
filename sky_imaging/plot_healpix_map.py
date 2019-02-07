#!/usr/bin/python

import numpy as np
import healpy as hp
import sys
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

    # Fast way to limit the plot region
    if len(ra_range) == 2 and len(dec_range) == 2 and coords == 'equitorial':
        tile_bounds_radec = [
                             [ra_range[0], dec_range[0]],
                             [ra_range[0], dec_range[1]],
                             [ra_range[1], dec_range[1]],
                             [ra_range[1], dec_range[0]]
                             ]
        tile_bounds_vec = np.array(
            [hp.pixelfunc.ang2vec(corner[0], corner[1], lonlat=True) for
                corner in tile_bounds_radec]
            )
        pixels_in_boundary = hp.query_polygon(
            map.nside, tile_bounds_vec, nest=map.nest
            )
        use_pixels = list(set(pixels_in_boundary) & set(map.pix_arr))
        use_signal = [
            map.signal_arr[(list(map.pix_arr)).index(pixel_val)]
            for pixel_val in use_pixels
            ]
        use_map = healpix_utils.HealpixMap(use_signal, use_pixels, map.nside, nest=map.nest,
                             coords=coords)
    # Slow way to limit the plot region
    elif len(ra_range) == 2 or len(dec_range) == 2:
        map.get_ra_dec(ra_cut=ra_cut)
        if len(ra_range) != 2:
            ra_range = [min(map.ra_arr), max(map.ra_arr)]
        if len(dec_range) != 2:
            dec_range = [min(map.dec_arr), max(map.dec_arr)]
        use_pixels = [
            map.pix_arr[ind] for ind in range(len(map.pix_arr))
            if ra_range[0] <= map.ra_arr[ind] <= ra_range[1]
            and dec_range[0] <= map.dec_arr[ind] <= dec_range[1]
            ]
        use_signal = [
            map.signal_arr[map.pix_arr.index(pixel_val)]
            for pixel_val in use_pixels
            ]
        use_map = HealpixMap(use_signal, use_pixels, map.nside, nest=map.nest,
                             coords=coords)
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
    # This part is buggy, fix it
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
        collection.set_clim(vmin=min(colors), vmax=max(colors))

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
    plt.savefig(save_filename, format='png', dpi=500)


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


if __name__ == '__main__':

    #map = healpix_utils.combine_maps_nearest_data(
    #    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Jan2019',
    #    nside=512, cube_name='Residual_{}'.format(cube_name)
    #)
    #map.write_data_to_fits(
    #    '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/Stokes{}_60obs_combined_nofilter.fits'.format(cube_name)
    #)
    map = healpix_utils.load_map('/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_60obs_combined_nofilter.fits')
    map.resample(128)
    plot_filled_pixels(
        map,
        '/Users/rubybyrne/diffuse_survey_plotting_Feb2019/StokesQ_60obs_combined_plot.png',
        ra_range = [-30, 135], dec_range=[-60, 5], colorbar_range=[-.02, .08]
    )
