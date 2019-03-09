#!/usr/bin/python

import numpy as np
import os
import sys
sys.path.append('../sky_imaging/')
import healpix_utils
import qa_utils


def average_stokes_v(
    fhd_run_path, qa_params_filepath
):

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


if __name__=='__main__':
    average_stokes_v(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        '/Users/rubybyrne/diffuse_survey_qa/diffuse_qa_params.csv'
    )
