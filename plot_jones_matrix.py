#!/usr/bin/python

import scipy.io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def main():

    jones_path = '/Users/ruby/Astro/jones_matrix_plotting/zenith_jones.sav'
    jones_struct = scipy.io.readsav(jones_path)['jones']
    obs_path = '/Users/ruby/Astro/jones_matrix_plotting/1131454296_obs.sav'
    obs_struct = scipy.io.readsav(obs_path)['obs']
    degpix = obs_struct['degpix'][0]
    dim = np.shape(jones_struct[0,0])[0]
    print(dim)

    jones_p_amp = np.sqrt(
        np.abs(jones_struct[0,0])**2. + np.abs(jones_struct[0,1])**2.
    )
    jones_p_amp[np.where(jones_p_amp == 0.)] = np.nan
    print(jones_p_amp[1000,1000])
    jones_q_amp = np.sqrt(
        np.abs(jones_struct[1,0])**2. + np.abs(jones_struct[1,1])**2.
    )
    k_mat = np.zeros((2, 2, dim, dim))
    for ind in range(2):
        k_mat[0, ind, :, :] = np.abs((jones_struct[0, ind])/jones_p_amp)
        k_mat[1, ind, :, :] = np.abs((jones_struct[1, ind])/jones_q_amp)
    jones_p_amp /= np.nanmax(jones_p_amp)
    jones_q_amp /= np.nanmax(jones_q_amp)

    use_cmap = matplotlib.cm.get_cmap('viridis')
    use_cmap.set_bad(color='grey')
    plt.imshow(
        jones_p_amp, origin='lower', interpolation='none', vmin=0, vmax=1,
        cmap=use_cmap,
        extent = [-dim*degpix/2, dim*degpix/2, -dim*degpix/2, dim*degpix/2]
    )
    plt.show()


if __name__=='__main__':
    main()
