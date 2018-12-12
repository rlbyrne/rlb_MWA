#!/usr/bin/python

from astropy.io import fits
import numpy as np
import healpy as hp
import sys
import math
import matplotlib.pyplot as plt


class ImageFromFits:

    def __init__(self, signal_arr, ra_axis=None, dec_axis=None,
                 ra_range=None, dec_range=None):
        n_ra_vals, n_dec_vals = np.shape(signal_arr)
        if ra_axis is not None:
            ra_axis = list(ra_axis)
            if len(ra_axis) != n_ra_vals:
                print 'ERROR: Number of axis elements does not match data axis. Exiting.'
                sys.exit(1)
        if dec_axis is not None:
            dec_axis = list(dec_axis)
            if len(dec_axis) != n_dec_vals:
                print 'ERROR: Number of axis elements does not match data axis. Exiting.'
                sys.exit(1)
        if ra_axis is None and ra_range is not None:
            if len(ra_range) != 2:
                print 'ERROR: parameters ra_range and dec_range must have the form [min, max]. Exiting.'
                sys.exit(1)
            ra_axis = np.linspace(ra_range[0], ra_range[1], n_ra_vals)
        if dec_axis is None and dec_range is not None:
            if len(dec_range) != 2:
                print 'ERROR: parameters ra_range and dec_range must have the form [min, max]. Exiting.'
                sys.exit(1)
            dec_axis = np.linspace(dec_range[0], dec_range[1], n_dec_vals)
        self.signal_arr = signal_arr
        self.ra_axis = ra_axis
        self.dec_axis = dec_axis

    def limit_data_range(self, ra_range=None, dec_range=None):
        if ra_range is not None:
            use_ra_inds = [i for i in range(len(self.ra_axis))
                           if ra_range[0] < self.ra_axis[i] < ra_range[1]
                           ]
        else:
            use_ra_inds = range(len(self.ra_axis))
        if dec_range is not None:
            use_dec_inds = [i for i in range(len(self.dec_axis))
                            if dec_range[0] < self.dec_axis[i] < dec_range[1]
                            ]
        else:
            use_dec_inds = range(len(self.dec_axis))
        self.signal_arr = self.signal_arr[
            use_dec_inds[0]:use_dec_inds[-1]+1,
            use_ra_inds[0]:use_ra_inds[-1]+1
        ]
        self.ra_axis = np.linspace(ra_range[0], ra_range[1], len(use_ra_inds))
        self.dec_axis = np.linspace(
            dec_range[0], dec_range[1], len(use_dec_inds)
        )


def load_image(data_filename):

    contents = fits.open(data_filename)
    use_hdu = 0
    data = contents[use_hdu].data
    header = contents[use_hdu].header

    if 'CD1_1' in header.keys() and 'CD2_2' in header.keys():  # FHD convention
        cdelt1 = header['CD1_1']
        cdelt2 = header['CD2_2']
        print 'WARNING: Ignoring curved sky effects.'
    elif 'CDELT1' in header.keys() and 'CDELT2' in header.keys():
        cdelt1 = header['CDELT1']
        cdelt2 = header['CDELT2']
    else:
        'ERROR: Header format not recognized.'
        sys.exit(1)

    ra_axis = [
        header['crval1'] +
        cdelt1*(i-header['crpix1'])
        for i in range(header['naxis1'])
        ]
    dec_axis = [
        header['crval2'] +
        cdelt2*(i-header['crpix2'])
        for i in range(header['naxis2'])
        ]

    fits_image = ImageFromFits(data, ra_axis=ra_axis, dec_axis=dec_axis)
    return fits_image


def difference_images(image1, image2):

    if image1.ra_axis == image2.ra_axis and image1.dec_axis == image2.dec_axis:
        data_diff = ImageFromFits(
            np.subtract(image1.signal_arr, image2.signal_arr),
            ra_axis=image1.ra_axis, dec_axis=image1.dec_axis
            )
    else:
        print 'WARNING: Image axes do not match. Interpolating image2 to image1 axes.'
        image2_signal_array_interp = griddata(
            (image2.ra_axis, image2.dec_axis),
            image2.signal_arr,
            (image1.ra_axis, image1.dec_axis),
            method='linear'
            )
        data_diff = ImageFromFits(
            np.subtract(image1.signal_arr, image2_signal_array_interp),
            ra_axis=image1.ra_axis, dec_axis=image1.dec_axis
            )
    return data_diff


def plot_fits_image(
    fits_image, save_filename='', title='', ra_range=None, dec_range=None,
    colorbar_range=[None, None], log=False,
    colorbar_label='Flux Density (Jy/sr)'
):

    if ra_range is not None or dec_range is not None:
        fits_image.limit_data_range(ra_range=ra_range, dec_range=dec_range)

    fig, ax = plt.subplots()
    plt.imshow(
        fits_image.signal_arr, origin='lower', interpolation='none',
        cmap='Greys_r',
        extent=[
            fits_image.ra_axis[0], fits_image.ra_axis[-1],
            fits_image.dec_axis[0], fits_image.dec_axis[-1]
            ],
        vmin=colorbar_range[0], vmax=colorbar_range[1], aspect='auto'
    )
    plt.axis('equal')
    ax.set_facecolor('gray')  # make plot background gray
    plt.grid(which='both', zorder=10, lw=0.5)
    cbar = plt.colorbar()
    # Label colorbar:
    cbar.ax.set_ylabel('Flux Density (Jy/sr)', rotation=270, labelpad=15)
    if save_filename == '':
        plt.show()
    else:
        plt.savefig(save_filename, format='png', dpi=500)


if __name__ == '__main__':

    #image_dirty = load_image('/Users/ruby/EoR/gaussian_model_debugging_Dec18/1130776864_uniform_Dirty_XX.fits')
    #image_model = load_image('/Users/ruby/EoR/gaussian_model_debugging_Dec18/1130776864_uniform_Model_XX.fits')
    #residual = difference_images(image_dirty, image_model)
    #plot_fits_image(residual, ra_range=[50, 52], dec_range=[-38, -36.5], save_filename='/Users/ruby/EoR/gaussian_model_debugging_Dec18/residual.png')
    image = load_image('/Users/ruby/Downloads/1130776864_uniform_Model_XX.fits')
    plot_fits_image(image, ra_range=[50, 52], dec_range=[-38, -36.5])
