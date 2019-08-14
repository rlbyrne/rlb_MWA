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
"""
metadata = UVData()
metadata.read_uvfits('/Users/ruby/EoR/1061316296.uvfits', read_data=False)
data = UVData()
data.read_uvfits(
    '/Users/ruby/EoR/1061316296.uvfits', times=metadata.time_array[200000],
    polarizations=0  # supports only one polarization
)

# Reformat data array:
data_vis_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
    np.nan, dtype=complex
)
for bl_ind in range(data.Nbls):
    data_vis_array[
        data.ant_1_array[bl_ind], data.ant_2_array[bl_ind], :, :
    ] = data.data_array[bl_ind, 0, :, :]
    data_vis_array[
        data.ant_2_array[bl_ind], data.ant_1_array[bl_ind], :, :
    ] = np.conj(data.data_array[bl_ind, 0, :, :])
# Reformat model array:
model_vis_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
    np.nan, dtype=complex
)
for bl_ind in range(model.Nbls):
    model_vis_array[
        model.ant_1_array[bl_ind], model.ant_2_array[bl_ind], :, :
    ] = model.data_array[bl_ind, 0, :, :]
    model_vis_array[
        model.ant_2_array[bl_ind], model.ant_1_array[bl_ind], :, :
    ] = np.conj(model.data_array[bl_ind, 0, :, :])
# Initialize weights:
vis_weights_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs),
    1., dtype=float
)
vis_weights_array[np.where(np.isnan(data_vis_array))[[0, 1, 2, 4]]] = 0.
vis_weights_array[np.where(np.isnan(model_vis_array))[[0, 1, 2, 4]]] = 0.
"""

def calibration_routine(
    data_vis_array, model_vis_array, vis_weights_array, delay_calibrate=False,
    convergence_threshold=.2, max_iter=100, gain_factor=.1, use_autos=False
):

    Nants, Ntimes, Nfreqs = data_vis_array.shape[1:4]

    # Initialize gains to 1
    gains = np.full((Nants, Nfreqs), 1.)

    if not use_autos:
        for ind in range(Nants):
            vis_weights_array[ind, ind, :, :] = 0.

    for mode in range(Nfreqs):
        gains_ant1_mat = np.repeat(np.repeat(
            gains[:, mode, np.newaxis], Nants, axis=1
            )[:, :, np.newaxis], Ntimes, axis=2
        )
        gains_ant2_mat = np.transpose(gains_ant1_mat, axes=[1, 0, 2])
        chi_squared = np.sum(vis_weights_array * np.abs(
            data_vis_array[:, :, :, mode]
            - gains_ant1_mat * np.conj(gains_ant2_mat)
            * model_vis_array[:, :, :, mode]
        )**2.)
        print "Initial X-squared for mode {}: {}".format(mode+1, chi_squared)

        convergence = 1.  # initialize convergence parameter
        iter = 0  # initialize iter parameter
        while iter <= max_iter and convergence >= convergence_threshold:

            # Calculate the gain amplitude gradient:
            gain_amp_grad = 4*np.sum(
                vis_weights_array*np.abs(gains_ant1_mat)**2.
                * np.abs(model_vis_array[:, :, :, mode])**2.
                * np.abs(gains_ant2_mat),
                axis=(0, 2, 3)
            ) - 4*np.sum(
                vis_weights_array*np.real(
                    gains_ant1_mat * model_vis_array[:, :, :, mode]
                    * np.conj(data_vis_array[:, :, :, mode])
                    * np.conj(gains_ant2_mat)/np.abs(gains_ant2_mat)
                ), axis=(0, 2, 3)
            )
            print gain_amp_grad[0]
            # Calculate the gain phase gradient:
            gain_phase_grad = 4*np.sum(
                vis_weights_array*np.imag(data_vis_array[:, :, :, mode])
                * np.real(
                    gains_ant1_mat*np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                ), axis=(0, 2, 3)
            ) - 4*np.sum(
                vis_weights_array*np.real(data_vis_array[:, :, :, mode])
                * np.imag(
                    gains_ant1_mat*np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                ), axis=(0, 2, 3)
            ) - 4*np.sum(np.diagonal(
                vis_weights_array*np.imag(data_vis_array[:, :, :, mode])
                * np.real(
                    gains_ant1_mat*np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                ), axis1=0, axis2=1
            )) + 4*np.sum(np.diagonal(
                vis_weights_array*np.real(data_vis_array[:, :, :, mode])
                * np.imag(
                    gains_ant1_mat*np.conj(gains_ant2_mat)
                    * model_vis_array[:, :, :, mode]
                ), axis1=0, axis2=1
            ))
            new_gains_amp = (
                np.abs(gains[:, mode])-chi_squared/gain_amp_grad*gain_factor
            )
            new_gains_phase = (
                np.angle(gains[:, mode])
                - chi_squared/gain_phase_grad*gain_factor
            )
            gains[:, mode] = new_gains_amp * np.exp(
                np.complex(0, 1) * new_gains_phase
            )
            gains_ant1_mat = np.repeat(np.repeat(
                gains[:, mode, np.newaxis], Nants, axis=1
                )[:, :, np.newaxis], Ntimes, axis=2
            )
            gains_ant2_mat = np.transpose(gains_ant1_mat, axes=[1, 0, 2])
            break
            """
            new_chi_squared = np.sum(vis_weights_array * np.abs(
                data_vis_array[:, :, :, mode, :]
                - gains_ant1_mat * np.conj(gains_ant2_mat)
                * model_vis_array[:, :, :, mode, :]
            )**2.)
            convergence = (new_chi_squared-chi_squared)/chi_squared
            chi_squared = new_chi_squared
            iter += 1
            """


if __name__ == '__main__':
    data_vis_array = np.full(
        (5, 5, 1, 3),
        1., dtype=complex
    )
    model_vis_array = np.full(
        (5, 5, 1, 3),
        1.1, dtype=complex
    )
    vis_weights_array = np.full(
        (5, 5, 1, 3),
        1., dtype=float
    )
    calibration_routine(
        data_vis_array, model_vis_array, vis_weights_array, delay_calibrate=False,
        convergence_threshold=.2, max_iter=100, gain_factor=.1, use_autos=False
    )
