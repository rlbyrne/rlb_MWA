#!/usr/bin/python

from astropy.io import fits
import numpy as np
import healpy as hp


contents = fits.open('/Users/ruby/Astro/diffuse_map.healfits')
nside = contents[0].header['nside']
ordering = contents[0].header['ordering']
signal_data = contents[0].data
freq = contents[0].header['crval2']  # Frequency in MHz
pixel_vals = contents[1].data['hpx_inds']
contents.close()

stokes_I = np.squeeze(signal_data[:, 0, 0])
stokes_Q = np.squeeze(signal_data[:, 0, 1])
stokes_U = np.squeeze(signal_data[:, 0, 2])
stokes_V = np.squeeze(signal_data[:, 0, 3])

coords = 'C'  # Map uses equitorial coordinates, you'll need this param for some healpy functions

if ordering.lower() == 'ring':
    nest = False
elif ordering.lower() == 'nested':
    nest = True

# Example of using healpy to calculate the pixel RA/Dec values
ra_arr, dec_arr = hp.pixelfunc.pix2ang(
    nside, pixel_vals, nest=nest, lonlat=True
)

# Example of converting Stokes I from explicit to implicit pixel ordering
signal_I_implicit = np.full(12*nside**2, hp.pixelfunc.UNSEEN)
signal_I_implicit[pixel_vals] = stokes_I
