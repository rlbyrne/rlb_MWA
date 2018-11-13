#!/usr/bin/python

# Obsolete script
# Replaced with healpix_utils.ppy and plot_healpix_map.py

import sys
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
sys.path.append('/Users/ruby/EoR/rlb_MWA/diffuse_survey_coverage_visualization')
import surveyview
import healpix_utils


def compare_observations():

    obs_list_file = '/Users/ruby/EoR/obs_lists/diffuse_survey_good_pointings_successful.txt'
    data_loc = '/Users/ruby/EoR/Healpix_fits'
    save_path = '/Users/ruby/Documents/2018 Fall Quarter/phys576b/jackknife_assignment'
    obs_list = open(obs_list_file, 'r').readlines()
    obs_list = [obs.strip() for obs in obs_list]  # strip newline characters
    obs_list = list(set(obs_list))  # remove duplicates
    observations = surveyview.load_survey(
        '/Users/ruby/EoR/sidelobe_survey_obsinfo.txt'
        )
    center_ra = 40.
    center_dec = -35.
    use_observations = []
    for obs in observations:
        if obs.obsid in obs_list:
            if (obs.ra-center_ra)**2.+ (obs.dec-center_dec)**2. < 10**2.:
                use_observations.append(obs)

    all_data = []
    for obs in use_observations:
        data, nside_old, nest = healpix_utils.load_map(
            '{}/{}_uniform_Residual_I_HEALPix.fits'.format(data_loc, obs.obsid)
            )
        #plot_filled_pixels(data, nside, nest,       '{}/{}_residual_I.png'.format(save_path, obs.obsid))
        nside = 256
        data = healpix_utils.healpix_downsample(data, nside_old, nside, nest)
        all_data.append(data)

    print 'Calculating data average...'
    ave_data, var_data, nsamples_data, ston_data = healpix_utils.average_healpix_maps(all_data)

    max_intersect_pixels = set([
        data_point.pixelnum for data_point in all_data[0]
        ])
    for data_set in all_data[1:]:
        max_intersect_pixels = max_intersect_pixels.intersection([
            data_point.pixelnum for data_point in data_set
            ])

    ave_data_sum = np.sum([data_point.signal for data_point in ave_data if data_point.pixelnum in max_intersect_pixels])
    all_data_norm = []
    for i, data_set in enumerate(all_data):
        data_sum = np.sum([data_point.signal for data_point in data_set if data_point.pixelnum in max_intersect_pixels])
        norm_factor = ave_data_sum/data_sum
        for point in data_set:
            point.signal = point.signal*norm_factor

    ave_data_norm, var_data_norm, nsamples_data, ston_data_norm = healpix_utils.average_healpix_maps(all_data)

    print 'Plotting data average...'
    plot_filled_pixels(ave_data_norm, nside, nest,       '{}/ave_norm_residual_I.png'.format(save_path))

    data_diff = healpix_utils.difference_healpix_maps(ave_data_norm, ave_data, nside, nest)

    plot_filled_pixels(data_diff, nside, nest,       '{}/ave_norm_minus_orig_residual_I.png'.format(save_path))
    #print 'Plotting data variance...'
    #plot_filled_pixels(var_data, nside, nest,       '{}/variance_residual_I.png'.format(save_path))
    #print 'Plotting data nsamples...'
    #plot_filled_pixels(nsamples_data, nside, nest,       '{}/nsamples_residual_I.png'.format(save_path))
    #plot_filled_pixels(ston_data, nside, nest,       '{}/ston_log_residual_I.png'.format(save_path))




def plot_filled_pixels(data, nside, nest, save_filename, title='', coords='equitorial'):

    #ra_min = -25
    #ra_max = 25
    #dec_min = -45
    #dec_max = -5
    ra_min = 30
    ra_max = 50
    dec_min = -45
    dec_max = -25

    # Only calculate pixel boundaries for pixels within the plot region
    # Only bother if data is in equitorial coordinates
    if coords == 'equatorial':
        tile_bounds_radec = [[ra_min, dec_min], [ra_min, dec_max],
                             [ra_max, dec_max], [ra_max, dec_min]]
        tile_bounds_vec = np.array(
            [hp.pixelfunc.ang2vec(corner[0], corner[1], lonlat=True) for
                corner in tile_bounds_radec]
            )
        use_pixels = hp.query_polygon(nside, tile_bounds_vec, nest=nest)
        data = [data_point for data_point in data if data_point.pixelnum in
                use_pixels]

    patches = []
    colors = []
    for point in data:

        if coords != 'equitorial':
            point.get_ra_dec(nside, nest, coords=coords)
            if not (ra_min < point.ra < ra_max and
                    dec_min < point.dec < dec_max):
                continue

        point.get_pixel_corners(nside, nest, coords=coords)
        polygon = Polygon(zip(point.pix_corner_ras, point.pix_corner_decs))
        patches.append(polygon)
        colors.append(point.signal)

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.05)
    collection.set_array(np.array(colors))  # set the data colors
    #collection.set_norm(LogNorm())
    #collection.set_clim(vmin=0, vmax=10)  # set the colorbar min and max
    collection.set_edgecolor('face')  # make the face and edge colors match

    fig, ax = plt.subplots(figsize=(10, 8), dpi=500)
    ax.add_collection(collection)  # plot data

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([ra_max, ra_min, dec_min, dec_max])
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270, labelpad=15)  # label colorbar
    #cbar.ax.set_ylabel('Number of Samples', rotation=270, labelpad=15)  # label colorbar

    plt.savefig(save_filename, format='png', dpi=500)


if __name__=='__main__':
    compare_observations()
