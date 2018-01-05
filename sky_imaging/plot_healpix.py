#!/usr/bin/python

# Based on code by Jared Canright, Univ of WA, and Kiana Henny, Whitman Col,
# 9/17


from mpl_toolkits.basemap import Basemap
import numpy as np
import healpy as hp
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import healpix_utils
from scipy.interpolate import griddata


def healpix_converter(data_filename):
    data_array = radcos_fio.data()
    # file pixels may be organized as ring or nest. Currently on nest
    ra = (np.pi * 2 / 360) * data_array[:, 0]
    dec = (np.pi * 2 / 360) * data_array[:, 1]
    flux = data_array[:, 2]
    ra[np.where(ra > np.pi)] -= 2 * np.pi
    pixel_refs = hp.pixelfunc.ang2pix(32, ra, dec, nest=False, lonlat=False)
    print pixel_refs

    return

    if nest_or_ring is 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    if nest_or_ring is 'nest':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)

    # plot of Galactic gas with coordinate projection
    min_ra = np.min(ra)
    print min_ra
    max_ra = np.max(ra)
    print max_ra
    mean_ra = np.mean(ra)
    min_dec = np.min(dec)
    print min_dec
    max_dec = np.max(dec)
    print max_dec
    mean_dec = np.mean(dec)

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # llc and urc are x & y ranges, and are specific to a location.
    # The latitude and longitude settings are part of basemap.
    m = Basemap(projection='hammer', llcrnrlon=-11, llcrnrlat=-15, urcrnrlon=13.5, urcrnrlat=-37, resolution='h', epsg=5520)
    # m = Basemap(projection='hammer', lon_0=mean_ra, lat_0=mean_dec, llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra, urcrnry=max_dec, resolution='h', epsg=5520)
    x, y = m(ra, dec)
    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    m.drawparallels(np.arange(-90., 120., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0., 420., 10.), labels=[0, 0, 0, 1])
    # creates a scatter plot of the selected data on a globe.
    m.scatter(x, y, 3, marker='o', linewidths=.1, c=data, cmap=plt.cm.coolwarm)
    m.colorbar()
    plt.show()


def plot_healpix_file(data_filename, save_filename):

    nside = 64
    data, nside_old, nest = healpix_utils.load_global_map(data_filename)
    print 'Downsampling...'
    data_downsampled = healpix_utils.healpix_downsample(data, nside_old, nside, nest)
    print 'Plotting...'
    #plot_filled_pixels(data_downsampled, nside, nest, save_filename, coords='galactic')
    plot_grid_interp(data_downsampled, nside, nest, save_filename, coords='galactic')


def overplot_haslam_contour():

    data, nside, nest = healpix_utils.load_map('/Users/ruby/EoR/mosaic_data.fits')

    fig, ax = plt.subplots(figsize=(24, 8), dpi=1000)

    # Define Healpix pixels to plot
    patches = []
    colors = []
    for point in data:
        point.get_pixel_corners(nside, nest)
        polygon = Polygon(zip(point.pix_corner_ras, point.pix_corner_decs))
        patches.append(polygon)
        colors.append(point.signal)

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.04)
    collection.set_array(np.array(colors))  # set the data colors
    collection.set_edgecolor('face')  # make the face and edge colors match
    ax.add_collection(collection)  # plot data

    # plot lines between tiles
    line_width = 2.0
    color = 'gray'
    order = 8
    plt.plot([130, -45], [-10, -10], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-20, -20], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-30, -30], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-40, -40], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-50, -50], lw=line_width, c=color, zorder=order)
    plt.plot([95, 95], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([85, 85], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([75, 75], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([65, 65], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([55, 55], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([45, 45], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([35, 35], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([25, 25], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([15, 15], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([5, 5], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([-5, -5], [-75, -20], lw=line_width, c=color, zorder=order)
    plt.plot([-15, -15], [-75, -20], lw=line_width, c=color, zorder=order)

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([110, -30, -50, 0])
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar

    data, nside, nest = healpix_utils.load_global_map('/Users/ruby/EoR/Healpix_fits/lambda_haslam408_dsds.fits')
    data = healpix_utils.healpix_downsample(data, nside, 512, nest)
    for point in data:
        point.get_ra_dec(nside, nest, coords='galactic')

    plt.tricontour(
        [point.ra for point in data],
        [point.dec for point in data],
        [point.signal for point in data],
        500, linewidths=0.5, colors='cyan')

    plt.savefig('/Users/ruby/EoR/haslam_overplot.png', format='png',
                dpi=1000)


def plot_healpix_tiling():

    data_type = 'Residual_I'
    normalization = 'uniform'
    data_dir = '/Users/ruby/EoR/Healpix_fits'

    tile_center_ras = [100, 100, 100, 100, 100,
                       90, 90, 90, 90, 90,
                       80, 80, 80, 80, 80,
                       70, 70, 70, 70, 70,
                       60, 60, 60, 60, 60,
                       50, 50, 50, 50, 50,
                       40, 40, 40, 40, 40,
                       30, 30, 30, 30, 30,
                       20, 20, 20, 20, 20,
                       10, 10, 10, 10, 10,
                       0, 0, 0, 0, 0,
                       -10, -10, -10, -10, -10,
                       -20, -20, -20, -20, -20]

    tile_center_decs = [-5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45,
                        -5, -15, -25, -35, -45]

    obsids = [1131478056, 1131564464, 1131477936, 1130787784, 1131477816,
              1131733672, 1131562544, 1131733552, 1130785864, 1131475896,
              1131559064, 1131558944, 1131731752, 1130782264, 1130783824,
              1131470736, 1131557144, 1131470616, 1130780464, 1131470496,
              1131468936, 1131553544, 1131726352, 1130778664, 1131468696,
              1131724672, 1131551744, 1131465216, 1130776864, 1130776624,
              1131463536, 1131549944, 1131463416, 1130773264, 1131463296,
              1131461616, 1131548024, 1131461496, 1130773144, 1131461376,
              1131717352, 1131544424, 1131717232, 1131459696, 1131459576,
              1131456216, 1131542624, 1131456096, 1131715432, 1131715312,
              1131711952, 1131540824, 1131454296, 1131713632, 1131454176,
              1131710152, 1131537224, 1131710032, 1131710032, 1131709912,
              1131535544, 1131535424, 1131535304, 1131710032, 1131709912]

    data = []
    nside = 512
    for i, obs in enumerate(obsids):
        print 'Gathering pixels from obsid {} of {}.'.format(i+1, len(obsids))
        obs_data, nside_old, nest = healpix_utils.load_map(
            '{}/{}_{}_{}_HEALPix.fits'.format(data_dir, obs, normalization,
                                              data_type)
            )
        obs_data = healpix_utils.healpix_downsample(obs_data, nside_old, nside, nest)
        tile_bounds_radec = [[tile_center_ras[i]-5, tile_center_decs[i]-5],
                             [tile_center_ras[i]-5, tile_center_decs[i]+5],
                             [tile_center_ras[i]+5, tile_center_decs[i]+5],
                             [tile_center_ras[i]+5, tile_center_decs[i]-5]]
        tile_bounds_vec = np.array(
            [hp.pixelfunc.ang2vec(corner[0], corner[1], lonlat=True) for
                corner in tile_bounds_radec]
            )
        use_pixels = hp.query_polygon(nside, tile_bounds_vec, nest=nest)
        data.extend(
            [data_point for data_point in obs_data if data_point.pixelnum in
                use_pixels]
            )

    #healpix_utils.write_data_to_fits(data, nside, nest, '/Users/ruby/EoR/mosaic_data.fits')

    # Define Healpix pixels to plot
    patches = []
    colors = []
    for point in data:
        point.get_pixel_corners(nside, nest)
        polygon = Polygon(zip(point.pix_corner_ras, point.pix_corner_decs))
        patches.append(polygon)
        colors.append(point.signal)

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.04)
    collection.set_array(np.array(colors))  # set the data colors
    collection.set_edgecolor('face')  # make the face and edge colors match

    fig, ax = plt.subplots(figsize=(24, 8), dpi=1000)
    ax.add_collection(collection)  # plot data

    # plot lines between tiles
    line_width = 2.0
    color = 'gray'
    order = 8
    plt.plot([130, -45], [-10, -10], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-20, -20], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-30, -30], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-40, -40], lw=line_width, c=color, zorder=order)
    plt.plot([130, -45], [-50, -50], lw=line_width, c=color, zorder=order)
    plt.plot([95, 95], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([85, 85], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([75, 75], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([65, 65], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([55, 55], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([45, 45], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([35, 35], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([25, 25], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([15, 15], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([5, 5], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([-5, -5], [-75, 20], lw=line_width, c=color, zorder=order)
    plt.plot([-15, -15], [-75, 20], lw=line_width, c=color, zorder=order)

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([110, -30, -50, 0])
    cbar = fig.colorbar(collection, ax=ax, extend='both')  # add colorbar
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar

    plt.savefig('/Users/ruby/EoR/Healpix_fits/mosaicplot.png', format='png',
                dpi=1000)


def plot_filled_pixels(data, nside, nest, save_filename, coords='equitorial'):

    # Collect Healpix pixels to plot
    patches = []
    colors = []
    for point in data:
        point.get_pixel_corners(nside, nest, coords=coords)
        polygon = Polygon(zip(point.pix_corner_ras, point.pix_corner_decs))
        patches.append(polygon)
        #if point.signal < 0.04:
        #    colors.append(point.signal)
        #else:
        #    colors.append(0.04)
        colors.append(point.signal)

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.04)
    collection.set_array(np.array(colors))  # set the data colors
    #collection.set_clim(vmin=-.02, vmax=.02)  # set the colorbar min and max
    collection.set_clim(vmin=0., vmax=60.)  # set the colorbar min and max
    collection.set_edgecolor('face')  # make the face and edge colors match

    fig, ax = plt.subplots(figsize=(21, 8), dpi=500)

    ax.add_collection(collection)  # plot data
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    #plt.axis([270, -90, -90, 90])
    plt.axis([105, -25, -50, 0])
    plt.grid(which='both', zorder=10, lw=0.5)
    #cbar = fig.colorbar(collection, ax=ax, extend='max')  # add colorbar
    cbar = fig.colorbar(collection, ax=ax, extend='max')
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar

    plt.savefig(save_filename, format='png', dpi=500)


def plot_grid_interp(data, nside, nest, save_filename, coords='equitorial'):

    stepsize = .1

    for point in data:
        point.get_ra_dec(nside, nest, coords=coords)
    ra_min = min([point.ra for point in data])
    ra_max = max([point.ra for point in data])
    dec_min = min([point.dec for point in data])
    dec_max = max([point.dec for point in data])

    grid_dec, grid_ra = np.mgrid[dec_min:dec_max:stepsize, ra_min:ra_max:stepsize]
    gridded_signal = griddata(([point.ra for point in data], [point.dec for point in data]), [point.signal for point in data], (grid_ra, grid_dec), method='linear')

    fig, ax = plt.subplots()
    plt.imshow(gridded_signal, origin='lower', interpolation='none', extent=[ra_min, ra_max, dec_min, dec_max], cmap='Greys_r', vmin=0., vmax=60.)
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.axis([ra_max, ra_min, dec_min, dec_max])
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar(extend='max')
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270)  # label colorbar
    plt.show()


if __name__ == '__main__':
    plot_healpix_file('/Users/ruby/EoR/Healpix_fits/lambda_haslam408_dsds.fits', '/Users/ruby/Desktop/haslam_nofilter.png')
