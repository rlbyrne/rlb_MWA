#!/usr/bin/python

import numpy as np
import os
import sys
sys.path.append('/Users/rubybyrne/rlb_MWA/sky_imaging/')
import healpix_utils
import qa_utils
import scipy.io


def average_stokes_v(fhd_run_path, qa_params_filepath):

    data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
    data_files = [
        file for file in data_files
        if '_uniform_Residual_V_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]

    ave_power = []
    for obsid in obs_list:
        map = healpix_utils.load_map(
            '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019/output_data/{}_uniform_Residual_V_HEALPix.fits'.format(
                obsid
            )
        )
        map.implicit_to_explicit_ordering()
        ave_power.append(np.mean(map.signal_arr**2.)/len(map.signal_arr))

    qa_utils.write_params_to_csv(
        qa_params_filepath, 'average Stokes V power', obs_list, ave_power,
        overwrite=True
    )


def get_obs_metadata(fhd_run_path, qa_params_filepath):

    data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
    data_files = [
        file for file in data_files
        if '_uniform_Residual_V_HEALPix.fits' in file
    ]
    obs_list = [file[0:10] for file in data_files]

    print 'Processing {} observations'.format(len(obs_list))

    n_sources = []
    min_decon_fluxes = []
    mean_gains = []
    cross_phases = []
    for obsid in obs_list:
        print 'Processing obsid {}'.format(obsid)
        decon_catalog = scipy.io.readsav(
            '{}/deconvolution/{}_fhd.sav'.format(fhd_run_path, obsid)
            )['source_array']
        n_sources.append(len(decon_catalog))
        min_decon_fluxes.append(
            min([source['flux']['I'][0] for source in decon_catalog])
        )
        cal = scipy.io.readsav(
            '{}/calibration/{}_cal.sav'.format(fhd_run_path, obsid)
            )['cal']
        [mean_gain_x, mean_gain_y] = cal['mean_gain'][0]
        mean_gains.append(np.mean([mean_gain_x, mean_gain_y]))
        cross_phases.append(cal['cross_phase'][0])

    qa_utils.write_params_to_csv(
        qa_params_filepath, 'deconvolved sources', obs_list, n_sources
    )
    qa_utils.write_params_to_csv(
        qa_params_filepath, 'min decon flux', obs_list, min_decon_fluxes
    )
    qa_utils.write_params_to_csv(
        qa_params_filepath, 'mean gain', obs_list, mean_gains
    )
    qa_utils.write_params_to_csv(
        qa_params_filepath, 'cross phase', obs_list, cross_phases
    )

if __name__=='__main__':
    get_obs_metadata(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        '/Users/rubybyrne/diffuse_survey_qa/diffuse_qa_params.csv'
    )
