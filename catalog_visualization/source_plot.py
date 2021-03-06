#!/usr/bin/python

# Script that plots catalogs

from astropy.io import fits
import numpy as np
import sys
import scipy.io
import astropy.io
import matplotlib
#matplotlib.use('Agg')  # use this if you don't have display access
import matplotlib.pyplot as plt


def plot_catalog(
    catalog_path, save_filename, flux_plot_max=None, title='',
    ra_range=None, dec_range=None, ra_cut_val=220, label_sources=False,
    label_gleam_supplement_sources=False
):

    catalog = scipy.io.readsav(catalog_path)['catalog']

    source_ras = []
    source_decs = []
    source_fluxes = []
    for source in catalog:
        ra = source['ra']
        if ra > ra_cut_val:
            ra -= 360.
        if ra < ra_cut_val-360.:
            ra += 360.
        source_ras.append(ra)
        source_decs.append(source['dec'])
        source_fluxes.append(source['flux']['I'][0])

    if flux_plot_max is None:
        flux_plot_max = max(source_fluxes)
    source_markersizes = []
    markersize_range = [.02, 3.]
    for flux in source_fluxes:
        if flux >= flux_plot_max:
            flux = flux_plot_max
        source_markersizes.append(
            flux/flux_plot_max*(markersize_range[1] - markersize_range[0])
            + markersize_range[0]
        )

    if ra_range is None:
        ra_min = min(source_ras)
        ra_max = max(source_ras)
        ra_range = [
            ra_min-(ra_max-ra_min)/10., ra_max+(ra_max-ra_min)/10.
        ]
    if dec_range is None:
        dec_min = min(source_decs)
        dec_max = max(source_decs)
        dec_range = [
            dec_min-(dec_max-dec_min)/10., dec_max+(dec_max-dec_min)/10.
        ]

    plt.figure()
    ax = plt.gca()
    plt.scatter(
        source_ras, source_decs, s=source_markersizes, facecolors='blue',
        edgecolors='none'
        )

    if label_sources:
        source_names = [
            'Crab', 'PicA', 'HydA', 'CenA',
            'HerA', 'VirA', 'CygA', 'CasA',
            '3C161', '3C353', '3C409',
            '3C444', 'CasA', 'ForA', 'HerA',
            'NGC0253', 'PicA', 'VirA',
            'PKS0349-27', 'PKS0442-28', 'PKS2153-69', 'PKS2331-41',
            'PKS2356-61', 'PKSJ0130-2610'
        ]
        named_source_ras = [
            83.6331, 79.9572, 139.524, 201.365,
            252.784, 187.706, 299.868, 350.858,
            96.7921, 260.117, 303.615,
            333.607, 350.866, 50.6738, 252.793,
            11.8991, 79.9541, 187.706,
            57.8988, 71.1571, 329.275, 353.609,
            359.768, 22.6158
        ]
        named_source_decs = [
            22.0145, -45.7788, -12.0956, -43.0192,
            4.9925, 12.3911, 40.7339, 58.8,
            -5.88472, -0.979722, 23.5814,
            -17.0267, 58.8117, -37.2083, 4.99806,
            -25.2886, -45.7649, 12.3786,
            -27.7431, -28.1653, -69.6900, -41.4233,
            -60.9164, -26.1656
        ]
        for ind, ra in enumerate(named_source_ras):
            if ra > ra_cut_val:
                named_source_ras[ind] = ra-360.
            if ra < ra_cut_val-360.:
                named_source_ras[ind] = ra+360.
        plt.scatter(
            named_source_ras, named_source_decs, marker='o',
            s=markersize_range[1], facecolors='none', edgecolors='red',
            linewidth=.1
        )
        for i, name in enumerate(source_names):
            plt.annotate(
                name, (named_source_ras[i], named_source_decs[i]),
                fontsize=3.
            )

    if label_gleam_supplement_sources:
        source_names = [
            'PicA', 'HydA', 'CenA', 'HerA', 'VirA', 'CasA', '3C161',
            '3C409', 'ForA'
        ]
        named_source_ras = [
            79.9572, 139.524, 201.365, 252.784, 187.706, 350.858, 96.7921, 303.615, 50.6738
        ]
        named_source_decs = [
            -45.7788, -12.0956, -43.0192, 4.9925, 12.3911, 58.8, -5.88472,
            23.5814, -37.2083
        ]
        for ind, ra in enumerate(named_source_ras):
            if ra > ra_cut_val:
                named_source_ras[ind] = ra-360.
            if ra < ra_cut_val-360.:
                named_source_ras[ind] = ra+360.
        plt.scatter(
            named_source_ras, named_source_decs, marker='o',
            s=markersize_range[1], facecolors='none', edgecolors='red',
            linewidth=.1
        )
        for i, name in enumerate(source_names):
            plt.annotate(
                name, (named_source_ras[i], named_source_decs[i]),
                fontsize=3.
            )

    plt.xlim(ra_range[1], ra_range[0])
    plt.ylim(dec_range[0], dec_range[1])
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title(title)
    ax.set_aspect('equal', adjustable='box')
    ax.set_facecolor('white')
    print 'Saving plot to {}'.format(save_filename)
    plt.savefig(save_filename, format='png', dpi=500)


def plot_extended_source_components(
    catalog_path, save_filename, source_ind=0, flux_plot_max=None, title='',
    ra_range=None, dec_range=None
):

    source = scipy.io.readsav(catalog_path)['catalog'][source_ind]
    source_ra = source['ra']
    source_dec = source['dec']
    components = source['extend']
    if len(components) == 0:
        print 'WARNING: Source is not extended.'

    comp_ras = []
    comp_decs = []
    comp_flux = []
    for comp in components:
        comp_ras.append(comp['ra'])
        comp_decs.append(comp['dec'])
        comp_flux.append(comp['flux']['I'][0])

    if flux_plot_max is None:
        flux_plot_max = max(comp_flux)
    comp_markersizes = []
    markersize_range = [.3, 3.]
    for flux in comp_flux:
        if flux >= flux_plot_max:
            flux = flux_plot_max
        comp_markersizes.append(
            flux/flux_plot_max*(markersize_range[1] - markersize_range[0])
            + markersize_range[0]
        )

    if ra_range is None:
        ra_min = min(comp_ras)
        ra_max = max(comp_ras)
        ra_range = [
            ra_min-(ra_max-ra_min)/10., ra_max+(ra_max-ra_min)/10.
        ]
    if dec_range is None:
        dec_min = min(comp_decs)
        dec_max = max(comp_decs)
        dec_range = [
            dec_min-(dec_max-dec_min)/10., dec_max+(dec_max-dec_min)/10.
        ]

    plt.figure()
    ax = plt.gca()
    plt.scatter(comp_ras, comp_decs, s=comp_markersizes, facecolors='white', edgecolors='none')
    plt.xlim(ra_range[1], ra_range[0])
    plt.ylim(dec_range[0], dec_range[1])
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title(title)
    ax.set_aspect('equal', adjustable='box')
    ax.set_facecolor('black')
    print 'Saving plot to {}'.format(save_filename)
    plt.savefig(save_filename, format='png', dpi=500)


def plot_gaussian_source_model(
    catalog_path, save_filename, source_ind=0, resolution=.01, title='',
    ra_range=None, dec_range=None, colorbar_range=[0., 1.],
    colorbar_label='Flux Density (Jy/sr)'
):
    # Note: this function does not support gaussian angles

    if ra_range is None:
        ra_range = [50, 51.25]
    if dec_range is None:
        dec_range = [-37.8, -36.7]

    grid_dec, grid_ra = np.mgrid[
        dec_range[0]:dec_range[1]:resolution,
        ra_range[0]:ra_range[1]:resolution
        ]
    plot_signal = np.zeros_like(grid_dec)

    source = scipy.io.readsav(catalog_path)['catalog'][source_ind]
    source_ra = source['ra']
    source_dec = source['dec']
    components = source['extend']
    if len(components) == 0:
        print 'WARNING: Source is not extended.'

    for comp in components:
        comp_ra = comp['ra']
        comp_dec = comp['dec']
        comp_flux = comp['flux']['I'][0]
        comp_size_x = comp['shape']['x'][0]/(7200.*np.sqrt(2*np.log(2.)))
        comp_size_y = comp['shape']['y'][0]/(7200.*np.sqrt(2*np.log(2.)))
        comp_size_angle = comp['shape']['angle'][0]

        if comp_size_x == 0:
            comp_size_x = resolution
        if comp_size_y == 0:
            comp_size_y = resolution

        for i in range(np.shape(grid_dec)[0]):
            for j in range(np.shape(grid_dec)[1]):
                pixel_val = (
                    comp_flux/(7200*np.pi*comp_size_x*comp_size_y)
                    * np.exp(-(grid_ra[i, j]-comp_ra)**2./(2*comp_size_x**2.))
                    * np.exp(-(grid_dec[i, j]-comp_dec)**2./(2*comp_size_y**2.))
                )
                plot_signal[i, j] += pixel_val

    fig, ax = plt.subplots(figsize=(21, 8), dpi=500)
    plt.imshow(
        plot_signal, origin='lower', interpolation='none',
        extent=[ra_range[0], ra_range[1], dec_range[0], dec_range[1]],
        cmap='Greys_r', vmin=colorbar_range[0], vmax=colorbar_range[1]
        )
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title(title)
    plt.axis('equal')
    ax.set_facecolor('black')  # make plot background gray
    plt.axis([ra_range[1], ra_range[0], dec_range[0], dec_range[1]])
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar(extend='max')
    cbar.ax.set_ylabel(colorbar_label, rotation=270, labelpad=15)  # label colorbar
    print 'Saving plot to {}'.format(save_filename)
    plt.savefig(save_filename, format='png', dpi=500)


if __name__ == '__main__':

    plot_catalog(
        '/Users/ruby/EoR/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav', '/Users/ruby/EoR/catalog_visualization/GLEAM_v2_plus_with_lables.png', title='', label_gleam_supplement_sources=True, dec_range=[-90, 90], flux_plot_max=100.
    )
