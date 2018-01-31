#!/usr/bin/python

from mpl_toolkits.basemap import Basemap
import numpy as np
import healpy as hp
import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
import healpix_utils
from scipy.interpolate import griddata
sys.path.append('/Users/ruby/EoR/rlb_MWA/diffuse_survey_coverage_visualization')
import surveyview


def calculate_mean_and_var():

    obs_list_file = '/Users/ruby/EoR/obs_lists/diffuse_survey_good_pointings_successful.txt'
    data_loc = '/Users/ruby/EoR/Healpix_fits'
    nside = 256

    obs_list = open(obs_list_file, 'r').readlines()
    obs_list = [obs.strip() for obs in obs_list]  # strip newline characters
    obs_list = list(set(obs_list))  # remove duplicates

    all_data = []
    for i, obs in enumerate(obs_list):
        print 'Gathering data from obsid {} of {}'.format(i+1, len(obs_list))
        print 'Loading data'
        data, nside_old, nest = healpix_utils.load_map(
            '{}/{}_uniform_Residual_I_HEALPix.fits'.format(data_loc, obs))
        data = healpix_utils.healpix_downsample(
            data, nside_old, nside, nest)

        print 'Assembling data'
        # Convert to implicit indexing
        signal_vals_implicit = [hp.pixelfunc.UNSEEN]*12*nside**2
        for j in range(len(data)):
            signal_vals_implicit[data[j].pixelnum] = data[j].signal

        all_data.append(signal_vals_implicit)

    print 'Calculating data statistics'
    data_ave = []
    data_var = []
    data_npoints = []
    for i in range(12*nside**2):
        use_data = [data_set[i] for data_set in all_data if data_set[i] != hp.pixelfunc.UNSEEN]
        n_points = len(use_data)
        if n_points > 0:
            data_npoints.append(healpix_utils.HealpixPixel(i, n_points))
            data_ave.append(healpix_utils.HealpixPixel(i, np.mean(use_data)))
            data_var.append(healpix_utils.HealpixPixel(i, np.var(use_data)))

    print 'Saving data statistics'
    healpix_utils.write_data_to_fits(data_ave, nside, nest,
                                     '/Users/ruby/Desktop/data_ave.fits')
    healpix_utils.write_data_to_fits(data_var, nside, nest,
                                     '/Users/ruby/Desktop/data_var.fits')
    healpix_utils.write_data_to_fits(data_npoints, nside, nest,
                                     '/Users/ruby/Desktop/data_npoints.fits')


def hist2d_datastat_plot():

    obs_list_file = '/Users/ruby/EoR/obs_lists/diffuse_survey_good_pointings_successful.txt'
    data_loc = '/Users/ruby/EoR/Healpix_fits'

    obs_list = open(obs_list_file, 'r').readlines()
    obs_list = [obs.strip() for obs in obs_list]  # strip newline characters
    obs_list = list(set(obs_list))  # remove duplicates
    #obs_list = obs_list[:4]
    observations = surveyview.load_survey('/Users/ruby/EoR/sidelobe_survey_obsinfo.txt')
    observations = [obs for obs in observations if obs.obsid in obs_list]

    var_data, nside, nest = healpix_utils.load_map(
        '/Users/ruby/Desktop/data_var_nside256.fits')
    var_pixelvals = [data_point.pixelnum for data_point in var_data]

    #observations = observations[:2]

    data_diff = []
    data_dist = []
    for i, obs in enumerate(observations):
        print 'Gathering data from obsid {} of {}'.format(i+1, len(obs_list))
        print 'Loading data'
        data, nside_old, nest = healpix_utils.load_map(
            '{}/{}_uniform_Residual_I_HEALPix.fits'.format(data_loc, obs.obsid))
        data = healpix_utils.healpix_downsample(
            data, nside_old, nside, nest)

        print 'Calculating data statistics'
        for data_point in data:
            data_diff.append(
                data_point.signal - var_data[
                    var_pixelvals.index(data_point.pixelnum)].signal
                )
            data_point.get_ra_dec(nside, nest)
            data_dist.append(
                hp.rotator.angdist([obs.ra, obs.dec],
                                   [data_point.ra, data_point.dec],
                                   lonlat=True)[0]/np.pi*180.
                )

    hist, dist_edges, diff_edges = np.histogram2d(data_dist, data_diff, bins=100)
    hist = hist.T

    for j in range(len(hist[0])):
        total = sum([hist[i][j] for i in range(len(hist))])
        if total != 0:
            for i in range(len(hist)):
                hist[i][j] = hist[i][j]/total

    fig, ax = plt.subplots()
    plt.imshow(
        hist, origin='lower', interpolation='none',
        extent=[dist_edges[0], dist_edges[-1], diff_edges[0], diff_edges[-1]],
        vmin=0, vmax=0.1)
    ax.set_aspect((dist_edges[-1]-dist_edges[0])/(diff_edges[-1]-diff_edges[0]))
    plt.xlabel('Distance from Pointing Center (deg)')
    plt.ylabel('Difference from Mean (Jy/sr)')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts, Normalized', rotation=270)  # label colorbar
    plt.savefig('/Users/ruby/Desktop/data_var_2dhist_testing.png', format='png', dpi=500)

if __name__ == '__main__':
    hist2d_datastat_plot()
