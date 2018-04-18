#!/usr/bin/python

# Script that finds a source from a deconvolution catalog and finds its flux
# in the residual images. Used to get a back-of-the-envelope estimate of the
# source polarization fraction.
# Note: This code does not handle RA wrapping (i.e. only works when sources
# are far away from RA=0).


from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import math
import matplotlib.pyplot as plt
import scipy.io


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
        num_comps = len(source[brightest_ind]['extend'])
        comp_ras = [0.]*num_comps
        comp_decs = [0.]*num_comps
        for i, comp in source[brightest_ind]['extend']:
            comp_ras[i] = comp['ra']
            comp_decs[i] = comp['dec']
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

    total_decon_stokes_I = 0.
    for i in range(len(source_array)):
        if ra_range[0] < source_ras[i] < ra_range[1] and dec_range[0] < source_decs[i] < dec_range[1]:
            total_decon_stokes_I += source_fluxes[i]

    print total_decon_stokes_I


if __name__=='__main__':
    main()
