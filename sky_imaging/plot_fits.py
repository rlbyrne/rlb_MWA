#!/usr/bin/python

from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import math
import matplotlib.pyplot as plt


def plot_file(data_filename):

    contents = fits.open(data_filename)
    data = contents[0].data

    fig, ax = plt.subplots()
    plt.imshow(data, origin='lower', interpolation='none', cmap='Greys_r',
               vmin=-.05, vmax=.05)
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar
    #plt.savefig(save_filename, format='png', dpi=500)
    plt.show()

def plot_difference(data_filename_1, data_filename_2):
    # Files must have exactly the same headers

    # Grab data
    contents_1 = fits.open(data_filename_1)
    data_1 = contents_1[0].data

    contents_2 = fits.open(data_filename_2)
    data_2 = contents_2[0].data

    # Get axes
    header = contents_1[0].header  # use the first file header for both
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

    data_diff = np.subtract(data_1, data_2)

    # Limit plot range
    ra_range = [30, 70]
    dec_range = [-60, -20]
    #ra_range = [min(ra_axis), max(ra_axis)]
    #dec_range = [min(dec_axis), max(dec_axis)]
    use_ra_inds = [i for i in range(len(ra_axis))
                   if ra_range[0] < ra_axis[i] < ra_range[1]
                   ]
    use_dec_inds = [i for i in range(len(dec_axis))
                    if dec_range[0] < dec_axis[i] < dec_range[1]
                    ]
    data_diff = data_diff[
                    use_dec_inds[0]:use_dec_inds[-1]+1,
                    use_ra_inds[0]:use_ra_inds[-1]+1
                    ]

    fig, ax = plt.subplots()
    #plt.imshow(data_diff, origin='lower', interpolation='none', cmap='Greys_r', vmin=-.00001, vmax=.00001)
    plt.imshow(data_diff, origin='lower', interpolation='none',
               cmap='Greys_r',
               extent=[ra_axis[use_ra_inds[0]], ra_axis[use_ra_inds[-1]],
                       dec_axis[use_dec_inds[0]], dec_axis[use_dec_inds[-1]]
                       ],
               vmin=-.01, vmax=.01, aspect='auto'
               )
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar
    plt.savefig('/Users/ruby/Desktop/4pol_decon_diff_Residual_V.png', format='png', dpi=1000)
    #plt.show()


if __name__=='__main__':
    #plot_file('/Users/ruby/EoR/full_pol_branch_ps_testing/workflow_3/sim_dirties/1061316296_uniform_Dirty_XX.fits')
    plot_difference('/Users/ruby/EoR/fhd_rlb_GLEAM+Fornax_cal_decon_4pol_Jan2018/output_data/1130776744_uniform_Residual_V.fits', '/Users/ruby/EoR/fhd_rlb_GLEAM+Fornax_cal_decon_4pol_fix_cross_phase_May2018/1130776744_uniform_Residual_V.fits')
