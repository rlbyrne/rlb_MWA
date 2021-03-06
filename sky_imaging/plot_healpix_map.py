#!/usr/bin/python

import numpy as np
import healpy as hp
import sys
import os
import matplotlib
#matplotlib.use('Agg')  # use this if you don't have display access
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
from astropy.io import fits
import healpix_utils
from scipy.interpolate import griddata

def get_mwa_beam(
    beam_filepath='/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_subtract_GLEAM_and_diffuse_Jun2020/output_data',
    beam_obsid = '1061316296'
):

    beam_xx_filepath = '{}/{}_Beam_XX.fits'.format(beam_filepath, beam_obsid)
    beam_yy_filepath = '{}/{}_Beam_YY.fits'.format(beam_filepath, beam_obsid)

    contents = fits.open(beam_xx_filepath)
    use_hdu = 0
    beam_xx_vals = contents[use_hdu].data.T  # transpose so convention is [RA, Dec]
    header = contents[use_hdu].header

    cdelt1 = header['CD1_1']
    cdelt2 = header['CD2_2']

    ra_axis = np.array([
        header['crval1'] +
        cdelt1*(i-header['crpix1'])
        for i in range(header['naxis1'])
        ])
    dec_axis = np.array([
        header['crval2'] +
        cdelt2*(i-header['crpix2'])
        for i in range(header['naxis2'])
        ])

    contents = fits.open(beam_yy_filepath)
    beam_yy_vals = contents[use_hdu].data.T

    beam_ras, beam_decs = np.meshgrid(ra_axis, dec_axis)
    beam_amp = np.real((beam_xx_vals+beam_yy_vals)/2.)
    beam_peak = np.where(beam_amp == np.max(beam_amp))
    beam_ras -= (beam_ras[beam_peak])[0]
    beam_decs -= (beam_decs[beam_peak])[0]

    return beam_ras, beam_decs, beam_amp


def plot_filled_pixels(
    map, save_filename=None, title='', ra_range=[], dec_range=[],
    colorbar_range=[], log=False, colorbar_label='Surface Brightness (Jy/sr)',
    ra_cut=None,
    overplot_points=False, point_ras=None, point_decs=None, point_values=None,
    overplot_points_vmin=-np.pi, overplot_points_vmax=np.pi,
    overplot_points_colormap='seismic',
    overplot_mwa_beam_contours=False,
    mwa_beam_center_ras=[0, 60], mwa_beam_center_decs=[-27, -27],
    overplot_hera_band=False,
    overplot_bright_sources=False,
    big=False
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
            ra_cut = (ra_range[1]+ra_range[0])/2.+12.
        else:
            ra_cut = 18.

    if overplot_points:
        point_ras[np.where(point_ras > ra_cut)] -= 24.
        point_ras[np.where(point_ras < ra_cut-24.)] += 24.

    map.get_ra_dec(ra_cut=ra_cut*15.)

    # Limit pixel calculation when axis ranges are set
    if len(ra_range) == 2 or len(dec_range) == 2:
        if len(ra_range) != 2:
            ra_range = np.array([np.min(map.ra_arr), np.max(map.ra_arr)])/15.
        if len(dec_range) != 2:
            dec_range = [np.min(map.dec_arr), np.max(map.dec_arr)]
        use_indices = np.arange(len(map.signal_arr))[
            (map.ra_arr/15.>ra_range[0]) & (map.ra_arr/15.<ra_range[1])
            & (map.dec_arr>dec_range[0]) & (map.dec_arr<dec_range[1])
        ]
        use_map = healpix_utils.HealpixMap(
            map.signal_arr[use_indices], map.pix_arr[use_indices], map.nside,
            nest=map.nest, coords=coords
        )
    else:
        use_map = map

    use_map.get_pixel_corners(ra_cut=ra_cut*15.)
    patches = []
    for ind in range(len(use_map.signal_arr)):
        polygon = Polygon(
            zip(use_map.pix_corner_ras_arr[ind]/15.,
                use_map.pix_corner_decs_arr[ind])
            )
        patches.append(polygon)
    colors = use_map.signal_arr

    # Establish axis ranges
    if len(ra_range) != 2:
        ra_range = np.array([
            np.amin(use_map.pix_corner_ras_arr),
            np.amax(use_map.pix_corner_ras_arr)
        ])/15.
    if len(dec_range) != 2:
        dec_range = np.array([
            np.amin(use_map.pix_corner_decs_arr),
            np.amax(use_map.pix_corner_decs_arr)
        ])

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

    plt.rcParams.update({'font.size': 9})
    if big:
        fig, ax = plt.subplots(figsize=(10, 4), dpi=600)
    else:
        fig, ax = plt.subplots(figsize=(6, 0.6*4), dpi=600)
    ax.add_collection(collection)  # plot data

    plt.xlabel('RA (hours)')
    plt.ylabel('Dec (degrees)')
    #plt.axis('equal')
    ax.set_aspect(1./15.)
    ax.set_facecolor('gray')  # make plot background gray
    if overplot_points:
        cm_overplot_points = plt.cm.get_cmap(overplot_points_colormap)
        plt.scatter(
            point_ras, point_decs, c=point_values,
            vmin=overplot_points_vmin, vmax=overplot_points_vmax, s=30,
            cmap=cm_overplot_points, edgecolor='black', linewidth=.5
        )
    if overplot_mwa_beam_contours:
        beam_ras, beam_decs, beam_val = get_mwa_beam()
        for beam_ind, use_beam_center_ra in enumerate(mwa_beam_center_ras):
            use_beam_center_dec = mwa_beam_center_decs[beam_ind]
            plt.contour(
                (beam_ras+use_beam_center_ra)/15., beam_decs+use_beam_center_dec,
                beam_val, levels=[.5],
                colors='cyan', linestyles=['solid'], linewidths=0.7
            )
    if overplot_hera_band:
        hera_band_center = -30.
        hera_band_width = 11.
        plt.plot(
            ra_range, np.full(2, hera_band_center+hera_band_width/2),
            '--', color='cyan', linewidth=0.7
        )
        plt.plot(
            ra_range, np.full(2, hera_band_center-hera_band_width/2),
            '--', color='cyan', linewidth=0.7
        )
    if overplot_bright_sources:
        source_names = [
            'Pictor A', 'Fornax A'
        ]
        named_source_ras = np.array([
            79.9572, 50.6738
        ])/15.
        named_source_decs = np.array([
            -45.7788, -37.2083
        ])
        plt.plot(
            named_source_ras, named_source_decs, 'x',
            color='yellow', markersize=3
        )
        for source_ind, name in enumerate(source_names):
            plt.annotate(
                name, (named_source_ras[source_ind]-2/15., named_source_decs[source_ind]),
                fontsize=8., color='black',
                path_effects=[patheffects.withStroke(linewidth=0.5, foreground="white")]
            )
        sgp_ra = 0.857222
        sgp_dec = -27.1283
        plt.plot(
            [sgp_ra], [sgp_dec], '+', color='red', markersize=6
        )
        plt.annotate(
            'SGP', (sgp_ra, sgp_dec-8), fontsize=8,
            horizontalalignment='center', color='black',
            path_effects=[patheffects.withStroke(linewidth=0.5, foreground="white")]
        )

    plt.axis([ra_range[1], ra_range[0], dec_range[0], dec_range[1]])
    plt.title(title)
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    # label colorbar
    cbar.ax.set_ylabel(colorbar_label, rotation=270, labelpad=15)
    if save_filename is not None:
        print 'Saving plot to {}'.format(save_filename)
        plt.savefig(save_filename, format='png', dpi=600)
        plt.close()
    else:
        plt.show()


def plot_grid_interp(
    map, save_filename=None, resolution=.1, title='', ra_range=[], dec_range=[],
    colorbar_range=[None, None], log=False,
    colorbar_label='Surface Brightness (Jy/sr)'
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


#if __name__ == '__main__':

    #map = healpix_utils.load_map(
    #    '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesI_averaged.fits'
    #)
    #plot_projection(
    #    map, title='Stokes I', colorbar_range=[-1.5,5],
    #    save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Dec2019/StokesI_averaged_alt_colorbar_mollweide.png'
    #)
