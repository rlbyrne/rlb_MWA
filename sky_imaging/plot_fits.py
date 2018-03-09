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

    contents_1 = fits.open(data_filename_1)
    data_1 = contents_1[0].data

    contents_2 = fits.open(data_filename_2)
    data_2 = contents_2[0].data

    data_diff = [[data_1[i][j]-data_2[i][j] for j in range(len(data_1[i]))] for i in range(len(data_1))]
    data_diff_flat = [data_diff[i][j] for j in range(len(data_diff[i])) for i in range(len(data_diff))]
    print min(data_diff_flat)
    print max(data_diff_flat)

    fig, ax = plt.subplots()
    plt.imshow(data_diff, origin='lower', interpolation='none', cmap='Greys_r', vmin=-.00001, vmax=.00001)
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar
    #plt.savefig('/Users/ruby/Desktop/dirty_2pol_XX_diff.png', format='png')
    plt.show()


if __name__=='__main__':
    #plot_file('/Users/ruby/EoR/full_pol_branch_ps_testing/workflow_3/sim_dirties/1061316296_uniform_Dirty_XX.fits')
    plot_difference('/Users/ruby/EoR/full_pol_branch_ps_testing/workflow_1/decon_dirties_no_recalc/1061316296_uniform_Dirty_XX.fits', '/Users/ruby/EoR/full_pol_branch_ps_testing/workflow_1/decon_dirties_no_recalc/stokes_I_sim_master_2pol_uniform_Dirty_XX.fits')
