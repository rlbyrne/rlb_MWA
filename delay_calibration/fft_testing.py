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


def calibrate_test_data():

    metadata = UVData()
    metadata.read_uvfits('/Users/ruby/EoR/1061316296.uvfits', read_data=False)
    data = UVData()
    data.read_uvfits(
        '/Users/ruby/EoR/1061316296.uvfits', times=metadata.time_array[200000],
        polarizations=-5  # supports only one polarization: -5=XX
    )
    # Data has shape (Nblts, Nspws, Nfreqs, Npols)
    data_times = np.unique(data.time_array)
    data_antennas = np.unique(np.append(data.ant_1_array, data.ant_2_array))

    # Reformat data array to have shape (Nants, Nants, Ntimes, Nfreqs)
    data_vis_array = np.full(
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

    # Use dummy data for the model:
    model = data

    # Reformat model array to have shape (Nants, Nants, Ntimes, Nfreqs)
    model_vis_array = np.full(
        (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
        np.nan, dtype=complex
    )
    for ind in range(data.Nblts):
        ant1 = np.where((data_antennas == model.ant_1_array[ind]))[0][0]
        ant2 = np.where((data_antennas == model.ant_2_array[ind]))[0][0]
        time = np.where((data_times == model.time_array[ind]))[0][0]
        model_vis_array[ant1, ant2, time, :] = model.data_array[ind, 0, :, 0]
        model_vis_array[ant2, ant1, time, :] = np.conj(
            model.data_array[ind, 0, :, 0]
        )

    model_vis_array *= (1.2+0.001*1j)  # make model different from data
    model_vis_array[0, 1, :, :] = 22.
    model_vis_array[1, 0, :, :] = 22.

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
    convergence_threshold=1, max_iter=10000, gain_factor=1e-3, use_autos=False
):

    floating_pt_error_threshold = 1e-7
    Nants, Ntimes, Nfreqs = data_vis_array.shape[1:4]

    # Initialize gains to 1
    gains = np.full((Nants, Nfreqs), 1., dtype=complex)

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

        convergence = 1.  # initialize convergence parameter
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
            )

            new_gains_real = np.real(gains[:, mode])
            new_gains_imag = np.imag(gains[:, mode])
            for ant_i in range(Nants):
                if abs(gain_grad_real[ant_i]) > floating_pt_error_threshold:
                    new_gains_real[ant_i] -= (
                        gain_grad_real[ant_i] / convergence * gain_factor
                    )
                if abs(gain_grad_imag[ant_i]) > floating_pt_error_threshold:
                    new_gains_imag[ant_i] -= (
                        gain_grad_imag[ant_i] / convergence * gain_factor
                    )
            gains[:, mode] = new_gains_real + 1j * new_gains_imag
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
            print new_chi_squared
            chi_squared = new_chi_squared
            iter += 1
        print gains[:,0]
        sys.exit()


if __name__ == '__main__':
    """data_vis_array = np.full(
        (2, 2, 1, 1),
        1., dtype=complex
    )
    data_vis_array[0, 1, 0, 0] = 1+1j
    data_vis_array[1, 0, 0, 0] = 1-1j
    model_vis_array = np.full(
        (2, 2, 1, 1),
        1.2, dtype=complex
    )
    model_vis_array[0, 1, 0, 0] = 1.1+1j
    model_vis_array[1, 0, 0, 0] = 1.1-1j
    vis_weights_array = np.full(
        (2, 2, 1, 1),
        1., dtype=float
    )
    calibration_routine(
        data_vis_array, model_vis_array, vis_weights_array, delay_calibrate=False,
        max_iter=100, use_autos=False
    )"""
    calibrate_test_data()
