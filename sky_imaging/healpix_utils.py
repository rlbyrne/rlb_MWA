#!/usr/bin/python

from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import math


class HealpixPixel:

    def __init__(self, pixelnum, signal):
        self.pixelnum = int(pixelnum)
        self.signal = float(signal)

    def get_ra_dec(self, nside, nest, coords='equitorial', ra_cut=270):
        if coords == 'galactic':
            theta_gal, phi_gal = hp.pixelfunc.pix2ang(nside, self.pixelnum,
                                                      nest=nest)
            rot = hp.rotator.Rotator(coord=['G', 'C'])
            theta_eq, phi_eq = rot(theta_gal, phi_gal)
            ra = phi_eq*180/math.pi
            dec = 90. - theta_eq*180/math.pi
        elif coords == 'equitorial':
            ra, dec = hp.pixelfunc.pix2ang(nside, self.pixelnum, nest=nest,
                                           lonlat=True)
        else:
            print 'ERROR: Coordinates must be galactic or equitorial.'
            sys.exit(1)

        if ra > ra_cut:
            self.ra = ra - 360.
        elif ra < ra_cut-360.:
            self.ra = ra + 360.
        else:
            self.ra = ra
        self.dec = dec

    def get_pixel_corners(self, nside, nest, coords='equitorial', ra_cut=270):
        corner_coords = hp.boundaries(nside, self.pixelnum, step=1,
                                      nest=nest)
        if coords == 'galactic':
            thetas_gal, phis_gal = hp.pixelfunc.vec2ang(
                np.transpose(corner_coords), lonlat=False)
            rot = hp.rotator.Rotator(coord=['G', 'C'])
            ras = [0]*len(thetas_gal)
            decs = [0]*len(thetas_gal)
            for corner in range(len(thetas_gal)):
                theta_eq, phi_eq = rot(thetas_gal[corner], phis_gal[corner])
                ras[corner] = phi_eq*180/math.pi
                decs[corner] = 90. - theta_eq*180/math.pi
        elif coords == 'equitorial':
            ras, decs = hp.pixelfunc.vec2ang(np.transpose(corner_coords),
                                             lonlat=True)
        else:
            print 'ERROR: Coordinates must be galactic or equitorial.'
            sys.exit(1)

        # Cluster pixel corners that are separated by the branch cut
        for i in range(1, len(ras)):
            if abs(ras[i]-ras[0]-360.) < abs(ras[i]-ras[0]):
                ras[i] -= 360.
            elif abs(ras[i]-ras[0]+360.) < abs(ras[i]-ras[0]):
                ras[i] += 360.

        if np.mean(ras) > ra_cut:
            ras = [ra-360. for ra in ras]
        elif np.mean(ras) < ra_cut-360.:
            ras = [ra+360. for ra in ras]

        self.pix_corner_ras = ras
        self.pix_corner_decs = decs


def load_map(data_filename):
    # Load a HEALPix map formatted with FHD conventions

    contents = fits.open(data_filename)
    nside = int(contents[1].header['nside'])
    ordering = contents[1].header['ordering']
    data = contents[1].data
    contents.close()

    signal_vals = data.field('SIGNAL')
    pixel_vals = data.field('PIXEL')

    if ordering.lower() == 'ring':
        nest = False
    elif ordering.lower() == 'nested':
        nest = True
    else:
        print 'ERROR: Invalid ordering parameter.'
        print 'Ordering must be "ring" or "nested". Exiting.'
        sys.exit(1)

    if len(pixel_vals) != len(signal_vals):
        print 'ERROR: Pixel index and data lengths do not match. Exiting.'
        sys.exit(1)

    pixel_data = []
    for i in range(len(pixel_vals)):
        data_point = HealpixPixel(pixel_vals[i], signal_vals[i])
        pixel_data.append(data_point)

    return pixel_data, nside, nest


def load_global_map(data_filename):
    # Load a HEALPix map formatted with Haslam-style conventions

    contents = fits.open(data_filename)
    nside = int(contents[1].header['nside'])
    ordering = contents[1].header['ordering']
    data = contents[1].data
    contents.close()

    if ordering.lower() == 'ring':
        nest = False
    elif ordering.lower() == 'nested':
        nest = True
    else:
        print 'ERROR: Invalid ordering parameter.'
        print 'Ordering must be "ring" or "nested". Exiting.'
        sys.exit(1)

    signal_vals = data.field('TEMPERATURE')

    pixel_data = []
    for i in range(len(signal_vals)):
        if signal_vals[i] != hp.pixelfunc.UNSEEN:  # Implicit indexing
            data_point = HealpixPixel(i, signal_vals[i])
            pixel_data.append(data_point)

    return pixel_data, nside, nest


def healpix_downsample(data, nside_in, nside_out, nest):

    if nest:
        ordering = 'nested'
    else:
        ordering = 'ring'

    # Convert data to implicit indexing
    signal_vals_implicit = [hp.pixelfunc.UNSEEN]*12*nside_in**2
    for i in range(len(data)):
        signal_vals_implicit[data[i].pixelnum] = data[i].signal

    downsampled_map = hp.pixelfunc.ud_grade(
        np.array(signal_vals_implicit), nside_out, pess=True,
        order_in=ordering)

    data_out = []
    for i in range(len(downsampled_map)):
        if downsampled_map[i] != hp.pixelfunc.UNSEEN:
            data_point = HealpixPixel(i, downsampled_map[i])
            data_out.append(data_point)

    return data_out


def write_data_to_fits(data, nside, nest, save_filename):

    signal_column = fits.Column(
        name='SIGNAL',
        array=np.array([data_point.signal for data_point in data]),
        format='1E'
        )
    pixelnum_column = fits.Column(
        name='PIXEL',
        array=np.array([data_point.pixelnum for data_point in data]),
        format='1J'
        )
    header = fits.Header()  # initialize header object
    header['nside'] = nside
    if nest:
        header['ordering'] = 'nested'
    else:
        header['ordering'] = 'ring'
    header['indxschm'] = 'explicit'

    hdu_0 = fits.PrimaryHDU()
    hdu_1 = fits.BinTableHDU.from_columns(
        [signal_column, pixelnum_column],
        header=header
        )
    hdu_list = fits.HDUList([hdu_0, hdu_1])
    hdu_list.writeto(save_filename)


def get_alm(data, nside, nest, lmax=None):

    # Convert data to implicit indexing
    signal_vals_implicit = [hp.pixelfunc.UNSEEN]*12*nside**2
    for i in range(len(data)):
        signal_vals_implicit[data[i].pixelnum] = data[i].signal

    # Convert to ring ordering if nested (map2alm doesn't support nested)
    if nest:
        signal_vals_implicit = hp.pixelfunc.reorder(
            signal_vals_implicit, n2r=True)

    alm = hp.sphtfunc.map2alm(signal_vals_implicit, lmax=lmax)
    alm = np.copy(alm)
    lmax = hp.sphtfunc.Alm.getlmax(len(alm))
    l, m = hp.sphtfunc.Alm.getlm(lmax)
    lm = np.zeros([len(alm), 2])
    lm[:, 0] = l
    lm[:, 1] = m

    return alm, lm


def filter_map(data, nside, nest, lmin=None, lmax=None, filter_width=0):
    # This function borrows from code by Miguel Morales

    if lmax is not None:
        alm_limit = int(math.ceil(lmax + filter_width/2.))
    else:
        alm_limit = None

    alm, lm = get_alm(data, nside, nest, lmax=alm_limit)

    # Define tukey windowing function
    window_fn = np.ones(int(max(lm[:, 0])+1))
    if lmin is not None:
        for i in range(int(math.ceil(lmin-filter_width/2.))):
            window_fn[i] = 0
        for i in range(int(math.ceil(lmin-filter_width/2.)),
                       int(math.floor(lmin+filter_width/2.))+1
                       ):
            window_fn[i] = (
                -math.cos((i-lmin+filter_width/2.)*math.pi/filter_width)+1
                )/2.
    if lmax is not None:
        for i in range(int(math.ceil(lmax-filter_width/2.)),
                       int(math.floor(lmax+filter_width/2.))+1
                       ):
            window_fn[i] = (
                math.cos((i-lmax+filter_width/2.)*math.pi/filter_width)+1
                )/2.
        for i in range(int(math.floor(lmax+filter_width/2.))+1,
                       len(window_fn)
                       ):
            window_fn[i] = 0

    # Apply windowing function
    for i in range(len(window_fn)):
        alm[np.where(lm[:, 0] == i)] = alm[
                                           np.where(lm[:, 0] == i)
                                           ]*window_fn[i]

    filtered_map = hp.sphtfunc.alm2map(alm, nside)

    # If input was nested ordering, convert to nested
    if nest:
        filtered_map = hp.pixelfunc.reorder(filtered_map, r2n=True)

    # Convert to explicit indexing
    output_data = []
    for i in range(len(filtered_map)):
        if filted_map[i] != hp.pixelfunc.UNSEEN:
            data_point = HealpixPixel(i, filted_map[i])
            output_data.append(data_point)

    return output_data
