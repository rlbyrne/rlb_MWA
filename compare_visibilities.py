#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from pyuvdata import UVData


def load_visibilities():

    fhd_run_path = '/Volumes/Bilbo/rlb_fhd_outputs/polarimetry_bug_testing_Apr2019'
    run_codes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    run_codes = run_codes[:5]
    pols = ['XX', 'XY', 'YX', 'YY']
    bls = [(0,100),(10,20),(40,60)]
    visibilities = np.zeros(
        (len(run_codes), len(bls), len(pols)), dtype=complex
    )
    ras = [
        353.85-360., 356.85-360., 359.85-360., 2.85, 5.85,
        353.85-360., 356.85-360., 359.85-360., 2.85, 5.85
    ]
    decs = [
        -26.7836, -26.7836, -26.7836, -26.7836, -26.7836,
        -24.7836, -24.7836, -24.7836, -24.7836, -24.7836
    ]

    for ind, code in enumerate(run_codes):
        fhd_files = [
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_XX.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_XY.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_YX.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_YY.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_model_XX.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_model_XY.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_model_YX.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_vis_model_YY.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/metadata/1131713632_settings.txt'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/metadata/1131713632_params.sav'.format(fhd_run_path, code),
            '{}/fhd_rlb_pol_bug_test_source_sim{}_Mar2019/vis_data/1131713632_flags.sav'.format(fhd_run_path, code)
        ]
        uv = UVData()
        uv.read_fhd(fhd_files, use_model=True)
        uv.select(freq_chans=[192], times=uv.time_array[0], bls=bls)
        visibilities[ind, :, :] = np.squeeze(uv.data_array)

    fig, ax = plt.subplots()
    for bl in range(len(bls)):
        for pol in range(len(pols)):
            plt.plot(ras[:len(run_codes)], np.abs(visibilities[:, bl, pol], '-o'))
    plt.savefig(
        '/Users/rubybyrne/visibilities_plot.png', format='png', dpi=500
    )


if __name__=='__main__':
    load_visibilities()
