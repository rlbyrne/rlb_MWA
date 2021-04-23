#!/usr/bin/python

import numpy as np
import healpy as hp
import sys
import os
import healpix_utils
import pyradiosky
from astropy.units import Quantity

sourcedir = '/Users/ruby/Downloads'
pols = ['I', 'Q', 'U', 'V']
maps = []
for pol_ind, pol_name in enumerate(pols):
    new_map = healpix_utils.load_map(
        '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
    )
    maps.append(new_map)

# Check that the pixel arrays are identical
for pol_ind in range(len(pols)-1):
    if np.max(np.abs(maps[pol_ind].pix_arr - maps[pol_ind+1].pix_arr)) != 0:
        print('ERROR: Mismatched pixel arrays.')
        sys.exit(1)

Ncomponents = np.shape(maps[0].signal_arr)[0]
stokes = np.zeros((4, 1, Ncomponents))
for pol_ind in range(len(pols)):
    stokes[pol_ind, :, :] = maps[pol_ind].signal_arr
stokes = Quantity(stokes, 'Jy/sr')
freq_array = Quantity([182000000], 'Hz')
if maps[0].nest:
    hpx_order = 'nested'
else:
    hpx_order = 'ring'

skymodel = pyradiosky.SkyModel(
    component_type='healpix',
    spectral_type='flat',
    stokes=stokes,
    freq_array=freq_array,
    hpx_inds=maps[0].pix_arr,
    hpx_order=hpx_order,
    nside=maps[0].nside
)
skymodel.write_skyh5('/Users/ruby/EoR/diffuse_map.skyh5')
