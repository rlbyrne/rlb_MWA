#!/usr/bin/python

# Based on code by Jared Canright, Univ of WA, and Kiana Henny, Whitman Col,
# 9/17


from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
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


def plot_filled_pixels(data_filename, save_filename):

    data = load_map(data_filename)
    fig, ax = plt.subplots()

    # Collect Healpix pixels to plot
    patches = []
    colors = []
    for point in data:
        point.get_pixel_corners()
        polygon = Polygon(zip(point.pix_corner_ras, point.pix_corner_decs))
        patches.append(polygon)
        colors.append(point.signal)

    collection = PatchCollection(patches, cmap='Greys_r', lw=0.04)
    collection.set_array(np.array(colors))  # set the data colors
    collection.set_edgecolor('face')  # make the face and edge colors match
    ax.add_collection(collection)  # plot data
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.axis('equal')
    ax.set_facecolor('black')  # make plot background black
    plt.axis([65, 40, -50, -30])
    plt.grid(which='both', zorder=10)
    fig.colorbar(collection, ax=ax)  # add colorbar

    plt.savefig(save_filename, format='png', dpi=2000)


def load_map(data_filename):

    contents = fits.open(data_filename)
    nside = contents[1].header['nside']
    ordering = contents[1].header['ordering']

    data = contents[1].data
    pixel_vals = data.field('PIXEL')
    signal_vals = data.field('SIGNAL')

    if len(pixel_vals) != len(signal_vals):
        print 'ERROR: Pixel index and data lengths do not match. Exiting.'
        sys.exit(1)

    pixel_data = []
    for i in range(len(pixel_vals)):
        data_point = HealpixPixel(nside, pixel_vals[i], ordering,
                                  signal_vals[i])
        data_point.get_ra_dec()
        pixel_data.append(data_point)

    return pixel_data


class HealpixPixel:

    def __init__(self, nside, pixelnum, ordering, signal):
        self.nside = int(nside)
        self.pixelnum = int(pixelnum)
        if ordering.lower() == 'ring':
            self.nest = False
        elif ordering.lower() == 'nested':
            self.nest = True
        else:
            print 'ERROR: Invalid ordering parameter.'
            print 'Ordering must be "ring" or "nested". Exiting.'
            sys.exit(1)
        self.signal = float(signal)

    def get_ra_dec(self):
        ra, dec = hp.pixelfunc.pix2ang(self.nside, self.pixelnum,
                                       nest=self.nest, lonlat=True)
        self.ra = ra
        self.dec = dec

    def get_pixel_corners(self):
        coords = hp.boundaries(self.nside, self.pixelnum, step=1,
                               nest=self.nest)
        ras, decs = hp.pixelfunc.vec2ang(np.transpose(coords), lonlat=True)
        self.pix_corner_ras = ras
        self.pix_corner_decs = decs


if __name__ == '__main__':
    plot_filled_pixels('/Users/ruby/EoR/1130776744_uniform_Residual_I_HEALPix.fits')
