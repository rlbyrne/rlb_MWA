#!/usr/bin/python

from pyuvdata import UVData
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

use_uvfits = False

UV = UVData()
fhd_prefix = '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Oct2019/fhd_rlb_single_source_test_branch_onefreq_Dec2019'
fhd_obs = '1130773144'
if use_uvfits:
    UV.read_uvfits('{}/model_vis.uvfits'.format(fhd_prefix))
else:
    fhd_files = [fhd_prefix + '/metadata/' + fhd_obs + '_settings.txt',
             fhd_prefix + '/metadata/' + fhd_obs + '_params.sav',
             fhd_prefix + '/vis_data/' + fhd_obs + '_flags.sav',
             fhd_prefix + '/vis_data/' + fhd_obs + '_vis_model_XX.sav',
             fhd_prefix + '/vis_data/' + fhd_obs + '_vis_model_YY.sav',
             fhd_prefix + '/vis_data/' + fhd_obs + '_vis_model_XY.sav',
             fhd_prefix + '/vis_data/' + fhd_obs + '_vis_model_YX.sav']
    UV.read_fhd(fhd_files, use_model=True)
    use_freq_hz = 167555000.
    UV.select(frequencies=use_freq_hz)
    UV.write_uvfits('{}/model_vis.uvfits'.format(fhd_prefix), spoof_nonessential=True)
