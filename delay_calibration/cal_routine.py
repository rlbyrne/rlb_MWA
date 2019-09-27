#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/ruby/EoR/pyuvdata')
from pyuvdata import UVData


"""coords = np.array(range(10))
input = np.exp(-coords**2.)
output = np.fft.fftshift(np.fft.fft(input))
output_coords = np.fft.fftshift(np.fft.fftfreq(len(coords)))

return_input = np.fft.ifft(np.fft.ifftshift(output))

print input
print input-return_input
plt.plot(output)
plt.plot(np.abs(return_input))
plt.show()
"""


def get_data(
    fhd_path='/Users/ruby/EoR/calibration_testing_Sept2019/fhd_rlb_single_source_test_single_beam_Sept2019',
    obsid='1130773144', pol='XX'
):

    filelist = ['{}/{}'.format(fhd_path, file) for file in [
        'vis_data/{}_vis_{}.sav'.format(obsid, pol),
        'vis_data/{}_vis_model_{}.sav'.format(obsid, pol),
        'vis_data/{}_flags.sav'.format(obsid),
        'metadata/{}_params.sav'.format(obsid),
        'metadata/{}_settings.txt'.format(obsid)
    ]]
    data = UVData()
    print 'Reading data...'
    data.read_fhd(filelist)
    model = UVData()
    print 'Reading model...'
    model.read_fhd(filelist, use_model=True)
    print 'Done.'

    # For testing, use one time only
    use_time = data.time_array[200000]
    data.select(times=use_time)
    model.select(times=use_time)

    data_times = np.unique(data.time_array)
    data_antennas = np.unique(np.append(data.ant_1_array, data.ant_2_array))

    # Data has shape (Nblts, Nspws, Nfreqs, Npols)
    # Reformat data array to have shape (Nants, Nants, Ntimes, Nfreqs)
    print 'Reformatting arrays...'
    data_vis_array = np.full(
        (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
        np.nan, dtype=complex
    )
    model_vis_array = np.full(
        (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
        np.nan, dtype=complex
    )
    for ind in range(data.Nblts):
        ant1 = np.where((data_antennas == data.ant_1_array[ind]))[0][0]
        ant2 = np.where((data_antennas == data.ant_2_array[ind]))[0][0]
        time = np.where((data_times == data.time_array[ind]))[0][0]
        data_vis_array[ant1, ant2, time, :] = data.data_array[ind, 0, :, 0]
        data_vis_array[ant2, ant1, time, :] = np.conj(
            data.data_array[ind, 0, :, 0]
        )
        model_vis_array[ant1, ant2, time, :] = model.data_array[ind, 0, :, 0]
        model_vis_array[ant2, ant1, time, :] = np.conj(
            model.data_array[ind, 0, :, 0]
        )
    print 'Done.'

    # Initialize weights:
    vis_weights_array = np.full(
        (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
        1., dtype=float
    )
    data_nans = np.where(np.isnan(data_vis_array))
    if len(data_nans[0]) > 0:
        vis_weights_array[data_nans[[0, 1, 2, 4]]] = 0.
    model_nans = np.where(np.isnan(model_vis_array))
    if len(model_nans[0]) > 0:
        vis_weights_array[model_nans[[0, 1, 2, 4]]] = 0.

    calibration_routine(data_vis_array, model_vis_array, vis_weights_array)


def calibration_routine(
    data_vis_array, model_vis_array, vis_weights_array, delay_calibrate=False,
    convergence_threshold=1, max_iter=1000, step_size=1e-3, use_autos=False,
    gains_init=None
):

    Nants, Ntimes, Nfreqs = data_vis_array.shape[1:4]

    # Initialize gains to 1
    if gains_init is None:
        gains = np.full((Nants, Nfreqs), 1., dtype=complex)
    else:
        gains = gains_init

    if not use_autos:  # autocalibration is currently not supported
        for ind in range(Nants):
            vis_weights_array[ind, ind, :] = 0.

    for mode in range(Nfreqs):
        gains_ant1_mat = np.repeat(np.repeat(
            gains[:, mode, np.newaxis], Nants, axis=1
            )[:, :, np.newaxis], Ntimes, axis=2
        )
        gains_ant2_mat = np.transpose(gains_ant1_mat, axes=[1, 0, 2])
        chi_squared = np.sum(vis_weights_array[:, :, :, mode] * np.abs(
            data_vis_array[:, :, :, mode]
            - gains_ant1_mat * np.conj(gains_ant2_mat)
            * model_vis_array[:, :, :, mode]
        )**2.)
        print "Initial X-squared for mode {}: {}".format(mode+1, chi_squared)

        # Initialize convergence parameter
        convergence = convergence_threshold+1
        iter = 0  # initialize iter parameter
        while iter <= max_iter and convergence >= convergence_threshold:

            # Calculate the gradient of the real part of the gains
            gain_grad_real = 4*np.sum(vis_weights_array[:, :, :, mode]*(
                np.real(gains_ant1_mat)
                * np.abs(gains_ant2_mat)**2.
                * np.abs(model_vis_array[:, :, :, mode])**2.
                - np.real(
                    np.conj(data_vis_array[:, :, :, mode])
                    * np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                )
                ), axis=(1, 2)
            )

            # Calculate the gradient of the imaginary part of the gains
            gain_grad_imag = 4*np.sum(vis_weights_array[:, :, :, mode]*(
                np.imag(gains_ant1_mat)
                * np.abs(gains_ant2_mat)**2.
                * np.abs(model_vis_array[:, :, :, mode])**2.
                + np.imag(
                    np.conj(data_vis_array[:, :, :, mode])
                    * np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                )
                ), axis=(1, 2)
            )

            convergence = np.sqrt(
                np.sum(gain_grad_real**2.) + np.sum(gain_grad_imag**2.)
            )/Nants

            # Update gains
            gains[:, mode] -= (
                (gain_grad_real + 1j * gain_grad_imag)
                / convergence * step_size
            )
            gains_ant1_mat = np.repeat(np.repeat(
                gains[:, mode, np.newaxis], Nants, axis=1
                )[:, :, np.newaxis], Ntimes, axis=2
            )
            gains_ant2_mat = np.transpose(gains_ant1_mat, axes=[1, 0, 2])
            new_chi_squared = np.sum(vis_weights_array[:, :, :, mode] * np.abs(
                data_vis_array[:, :, :, mode]
                - gains_ant1_mat * np.conj(gains_ant2_mat)
                * model_vis_array[:, :, :, mode]
            )**2.)
            chi_squared = new_chi_squared
            iter += 1
            print chi_squared
            print convergence
            print gains
        sys.exit()


if __name__ == '__main__':
    get_data()
