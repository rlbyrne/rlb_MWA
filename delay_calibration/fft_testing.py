#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from pyuvdata import UVData
import sys

coords = np.array(range(10))
input = np.exp(-coords**2.)
output = np.fft.fftshift(np.fft.fft(input))
output_coords = np.fft.fftshift(np.fft.fftfreq(len(coords)))

return_input = np.fft.ifft(np.fft.ifftshift(output))

print input
print input-return_input
plt.plot(output)
plt.plot(np.abs(return_input))
plt.show()

metadata = UVData()
metadata.read_uvfits('/Users/ruby/EoR/1061316296.uvfits', read_data=False)
data = UVData()
data.read_uvfits(
    '/Users/ruby/EoR/1061316296.uvfits', times=metadata.time_array[200000]
)

# Reformat data array:
data_vis_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs, data.Npols),
    np.nan, dtype=complex
)
for bl_ind in range(data.Nbls):
    data_vis_array[
        data.ant_1_array[bl_ind], data.ant_2_array[bl_ind], :, :, :
    ] = data.data_array[bl_ind, 0, :, :]
    data_vis_array[
        data.ant_2_array[bl_ind], data.ant_1_array[bl_ind], :, :, :
    ] = np.conj(data.data_array[bl_ind, 0, :, :])
# Reformat model array:
model_vis_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Nfreqs, data.Npols),
    np.nan, dtype=complex
)
for bl_ind in range(model.Nbls):
    model_vis_array[
        model.ant_1_array[bl_ind], model.ant_2_array[bl_ind], :, :, :
    ] = model.data_array[bl_ind, 0, :, :]
    model_vis_array[
        model.ant_2_array[bl_ind], model.ant_1_array[bl_ind], :, :, :
    ] = np.conj(model.data_array[bl_ind, 0, :, :])
# Initialize weights:
vis_weights_array = np.full(
    (data.Nants_data, data.Nants_data, data.Ntimes, data.Npols),
    1., dtype=float
)
vis_weights_array[np.where(np.isnan(data_vis_array))[[0, 1, 2, 4]]] = 0.
vis_weights_array[np.where(np.isnan(model_vis_array))[[0, 1, 2, 4]]] = 0.


def calibration_routine(
    data_vis_array, model_vis_array, delay_calibrate=False, convergence_threshold=.2, max_iter=100, use_autos=False
):

    # Initialize gains to 1
    gains = np.full((Nants, Nfreqs, Npols), 1.)
    iter = 0

    if not use_autos:
        for ind in range(Nants):
            vis_weights_array[ind, ind] = 0.

    if (data.ant_1_array != model.ant_1_array) or
    (data.ant_2_array != model.ant_2_array):
        print 'ERROR: Data and model antenna ordering does not match. Exiting.'
        sys.exit()

    while iter <= max_iter and convergence >= convergence_threshold:
        for mode in range(data.Nfreqs):
            gains_update = gains[:, :, mode]
            for ant in data.Nants_data:
                numerator = np.dot(model_update[])
                gains_update[ant] =
            data.data_array()
        iter += 1
