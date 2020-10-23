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
    map, save_filename=None, title='', ra_range=[], dec_range=[],
    colorbar_range=[], log=False, colorbar_label='Flux Density (Jy/sr)',
    ra_cut=None,
    overplot_points=False, point_ras=None, point_decs=None, point_values=None,
    overplot_points_vmin=-np.pi, overplot_points_vmax=np.pi,
    overplot_points_colormap='seismic'
):

    if map.coords == '':
        print 'WARNING: No map coordinate scheme supplied.'
        print 'Assuming equitorial coordinates.'
        coords = 'equitorial'
    else:
        coords = map.coords

    # Set branch cut location
    if ra_cut is None:
        if len(ra_range) == 2:
            ra_cut = (ra_range[1]+ra_range[0])/2.+180.
        else:
            ra_cut = 270.

    if overplot_points:
        point_ras[np.where(point_ras > ra_cut)] -= 360.
        point_ras[np.where(point_ras < ra_cut-360.)] += 360.

    map.get_ra_dec(ra_cut=ra_cut)
    print np.min(map.ra_arr)
    print np.max(map.ra_arr)
    print np.min(map.dec_arr)
    print np.max(map.dec_arr)

    # Limit pixel calculation when axis ranges are set
    if len(ra_range) == 2 or len(dec_range) == 2:
        map.get_ra_dec(ra_cut=ra_cut)
        if len(ra_range) != 2:
            ra_range = [np.min(map.ra_arr), np.max(map.ra_arr)]
        if len(dec_range) != 2:
            dec_range = [np.min(map.dec_arr), np.max(map.dec_arr)]
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
        ra_range = [
            np.amin(use_map.pix_corner_ras_arr),
            np.amax(use_map.pix_corner_ras_arr)
        ]
    if len(dec_range) != 2:
        dec_range = [
            np.amin(use_map.pix_corner_decs_arr),
            np.amax(use_map.pix_corner_decs_arr)
        ]

    cm_use = 'Greys_r'
    #cm_use = 'viridis'
    collection = PatchCollection(patches, cmap=cm_use, lw=0.05)
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

    fig, ax = plt.subplots(figsize=(10, 4), dpi=500)
    ax.add_collection(collection)  # plot data

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    if overplot_points:
        cm_overplot_points = plt.cm.get_cmap(overplot_points_colormap)
        plt.scatter(
            point_ras, point_decs, c=point_values,
            vmin=overplot_points_vmin, vmax=overplot_points_vmax, s=30,
            cmap=cm_overplot_points, edgecolor='black', linewidth=.5
        )
    plt.axis([ra_range[1], ra_range[0], dec_range[0], dec_range[1]])
    plt.title(title)
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    # label colorbar
    cbar.ax.set_ylabel(colorbar_label, rotation=270, labelpad=15)
    if save_filename is not None:
        print 'Saving plot to {}'.format(save_filename)
        plt.savefig(save_filename, format='png', dpi=300)
        plt.close()
    else:
        plt.show()


def plot_grid_interp(
    map, save_filename=None, resolution=.1, title='', ra_range=[], dec_range=[],
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
    if save_filename is not None:
        print 'Saving plot to {}'.format(save_filename)
        plt.savefig(save_filename, format='png', dpi=500)
    else:
        plt.show()


def plot_projection(
    map, save_filename=None, title='', colorbar_range=[None, None]
):

    map.explicit_to_implicit_ordering()

    proj = hp.mollview(
        map=map.signal_arr, coord=map.coords_healpy_conv, nest=map.nest,
        title=title, min=colorbar_range[0], max=colorbar_range[1], cbar=True,
        cmap='Greys_r', return_projected_map=False, notext=True, unit='Jy/sr'
    )
    hp.graticule()
    plt.xlim([-1.7,.7])
    plt.ylim([-1.1, .25])
    plt.title = title
    if save_filename is None:
        plt.show()
    else:
        print 'Saving plot to {}'.format(save_filename)
        plt.savefig(save_filename, dpi=300)
        plt.close()


if __name__ == '__main__':

    map = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesI_averaged.fits'
    )
    plot_projection(
        map, title='Stokes I', colorbar_range=[-1.5,5],
        save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesI_averaged_alt_colorbar_mollweide.png'
    )
