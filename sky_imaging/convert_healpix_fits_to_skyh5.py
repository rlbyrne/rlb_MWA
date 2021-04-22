#!/usr/bin/python

import numpy as np
import healpy as hp
import sys
import os
import healpix_utils

sourcedir = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
save_filename = '/Users/rubybyrne/diffuse_survey_plotting_Sept2020/polarized_diffuse_map.fits'
pols = ['I', 'Q', 'U', 'V']
maps = []
for pol_ind, pol_name in enumerate(pols):
    new_map = healpix_utils.load_map(
        '{}/Stokes{}_average_map_empirical_rm_in_eor0.fits'.format(sourcedir, pol_name)
    )
    maps.append(new_map)
healpix_utils.write_data_to_standard_fits(maps, save_filename)
