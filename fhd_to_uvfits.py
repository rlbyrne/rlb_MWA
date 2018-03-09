#!/usr/bin/python

from pyuvdata import UVData
UV=UVData() 
fhd_prefix = '/nfs/mwa-04/r1/EoRuvfits/DiffuseSurvey2015/fhd_rlb_fornax_sim_fullpol_4pol_Feb2018'
#fhd_obs = '1061316296'
fhd_obs = '1130776744'
fhd_files = [fhd_prefix + '/' + fhd_obs + '_settings.txt', 
	     fhd_prefix + '/metadata/' + fhd_obs + '_params.sav',
	     fhd_prefix + '/vis_data/' + fhd_obs + '_flags.sav',
	     fhd_prefix + '/vis_data/' + fhd_obs + '_vis_XX.sav', 
	     fhd_prefix + '/vis_data/' + fhd_obs + '_vis_YY.sav',
	     fhd_prefix + '/vis_data/' + fhd_obs + '_vis_XY.sav',
	     fhd_prefix + '/vis_data/' + fhd_obs + '_vis_YX.sav']
UV.read_fhd(fhd_files)
UV.write_uvfits('/nfs/eor-00/h1/rbyrne/sim_visibilities/fornax_sim_fullpol_4pol.uvfits', spoof_nonessential=True)

