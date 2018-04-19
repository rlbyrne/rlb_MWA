#!/usr/bin/python

# Script that finds a source from a deconvolution catalog and finds its flux
# in the residual images. Used to get a back-of-the-envelope estimate of the
# source polarization fraction.
# Note: This code does not handle RA wrapping (i.e. only works when sources
# are far away from RA=0).


from astropy.io import fits
import numpy as np
import sys
import math
import scipy.io
import astropy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():

    fhd_run_path = '/home/rlbyrne/fhd_rlb_GLEAM+Fornax_cal_decon_4pol_Jan2018'
    obsid = '1130776744'
    min_patch_size = .1  # minimum box dimension in degrees

    decon_data = '{}/deconvolution/{}_fhd.sav'.format(fhd_run_path, obsid)
    sav = scipy.io.readsav(decon_data)
    source_array = sav['source_array']
    num_sources = len(source_array)

    # Find the brightest source
    source_fluxes = [0.]*num_sources
    source_ras = [0.]*num_sources
    source_decs = [0.]*num_sources
    for i, source in enumerate(source_array):
        source_fluxes[i] = source['flux']['I'][0]
        source_ras[i] = source['ra']
        source_decs[i] = source['dec']
    brightest_ind = source_fluxes.index(max(source_fluxes))
    brightest_loc = [source_ras[brightest_ind], source_decs[brightest_ind]]

    # Define a rectangle to integrate flux over
    if source_array[brightest_ind]['extend'] is not None:
        num_comps = len(source_array[brightest_ind]['extend'])
        comp_ras = [0.]*num_comps
        comp_decs = [0.]*num_comps
        for i in range(num_comps):
            comp_ras[i] = source_array[brightest_ind]['extend'][i]['ra']
            comp_decs[i] = source_array[brightest_ind]['extend'][i]['dec']
        ra_range = [
            min([min(comp_ras), brightest_loc[0]-min_patch_size/2.]),
            max([max(comp_ras), brightest_loc[0]+min_patch_size/2.])
            ]
        dec_range = [
            min([min(comp_decs), brightest_loc[1]-min_patch_size/2.]),
            max([max(comp_decs), brightest_loc[1]+min_patch_size/2.])
            ]
    else:
        ra_range = [brightest_loc[0]-min_patch_size/2.,
                    brightest_loc[0]+min_patch_size/2.]
        dec_range = [brightest_loc[1]-min_patch_size/2.,
                     brightest_loc[1]+min_patch_size/2.]

    # Make sure to grab all sources within the defined boundaries
    total_decon_stokes_I = 0.
    for i in range(len(source_array)):
        if ra_range[0] < source_ras[i] < ra_range[1] and dec_range[0] < source_decs[i] < dec_range[1]:
            total_decon_stokes_I += source_fluxes[i]

    print total_decon_stokes_I

    pol_mode = ['I', 'Q', 'U', 'V']
    res_flux = [0.]*4
    for pol_ind, pol in enumerate(pol_mode):
        res_file = '{}/output_data/{}_uniform_Residual_{}.fits'.format(fhd_run_path, obsid, pol)
        contents = astropy.io.fits.open(res_file)
        data = contents[0].data
        header = contents[0].header
        ra_axis = [
            header['crval1'] +
            header['cd1_1']*(i-header['crpix1'])
            for i in range(header['naxis1'])
            ]
        dec_axis = [
            header['crval2'] +
            header['cd2_2']*(i-header['crpix2'])
            for i in range(header['naxis2'])
            ]

        use_ra_inds = [i for i in range(len(ra_axis))
                       if ra_range[0] < ra_axis[i] < ra_range[1]
                       ]
        use_dec_inds = [i for i in range(len(dec_axis))
                        if dec_range[0] < dec_axis[i] < dec_range[1]
                        ]
        data_cut = data[
                        use_ra_inds[0]:use_ra_inds[-1]+1,
                        use_dec_inds[0]:use_dec_inds[-1]+1
                        ]
        res_flux[pol_ind] = np.sum(data_cut)

        if pol == 'Q':
            fig, ax = plt.subplots()
            plt.imshow(data_cut, origin='lower', interpolation='none', cmap='Greys_r', extent=[ra_axis[use_ra_inds[0]], ra_axis[use_ra_inds[-1]], dec_axis[use_dec_inds[0]], dec_axis[use_dec_inds[-1]]])
            plt.axis('equal')
            ax.set_facecolor('gray')  # make plot background gray
            plt.grid(which='both', zorder=10, lw=0.5)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar
            plt.savefig('q_res_plot.png')

    print ra_range
    print dec_range
    print res_flux






if __name__=='__main__':
    main()
