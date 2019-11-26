#!/usr/bin/python

from pyuvdata import UVData
import scipy.io
import numpy as np
import sys

UV = UVData()
fhd_prefix = '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Oct2019/fhd_rlb_single_source_test_master_onefreq_Nov2019'
fhd_obs = '1130773144'
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
sys.exit()
start_phase_ra = UV.phase_center_ra
start_phase_dec = UV.phase_center_dec

# arrays are read in transposed
beamxx = scipy.io.readsav(
	'{}/beams_from_beam_image/single_source_sim_beam_xx.sav'.format(fhd_prefix)
)['beam_base'].T
beamyy = scipy.io.readsav(
	'{}/beams_from_beam_image/single_source_sim_beam_yy.sav'.format(fhd_prefix)
)['beam_base'].T
beamxy = scipy.io.readsav(
	'{}/beams_from_beam_image/single_source_sim_beam_xy.sav'.format(fhd_prefix)
)['beam_base'].T
beamyx = scipy.io.readsav(
	'{}/beams_from_beam_image/single_source_sim_beam_yx.sav'.format(fhd_prefix)
)['beam_base'].T

source_loc_x = 1230.39
source_loc_y = 1071.85
source_flux = 326.789
source_ra_deg = 50.4057
source_dec_deg = -37.1493

beam_val_xx = scipy.interpolate.interp2d(
	list(range(2048)), list(range(2048)), beamxx
)(source_loc_y, source_loc_x)[0]
apparent_flux = source_flux*beam_val_xx

print UV.data_array[:10,0,0,0]
start_data_array = np.copy(UV.data_array)
UV.unphase_to_drift()
UV.phase(np.deg2rad(source_ra_deg), np.deg2rad(source_dec_deg))
UV.data_array.fill(apparent_flux)
UV.unphase_to_drift()
UV.phase(start_phase_ra, start_phase_dec)

print UV.data_array[:10,0,0,0]

#UV.write_uvfits('/nfs/eor-00/h1/rbyrne/sim_visibilities/fornax_sim_fullpol_4pol.uvfits', spoof_nonessential=True)
