#!/usr/bin/python

import scipy.io
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def main():

    # Restore Jones matrix, as saved at the start of the rotate_jones_matrix
    # function. Jones matrix is in Az/ZA coordinates.
    jones_path = '/Users/ruby/Astro/jones_matrix_plotting/zenith_jones.sav'
    jones_struct = scipy.io.readsav(jones_path)['jones']
    jones_rot_path = '/Users/ruby/Astro/jones_matrix_plotting/zenith_jones_rot.sav'
    jones_rot_struct = scipy.io.readsav(jones_rot_path)['jones_rot']
    #obs_path = '/Users/ruby/Astro/jones_matrix_plotting/1131454296_obs.sav'
    #obs_struct = scipy.io.readsav(obs_path)['obs']
    #degpix = obs_struct['degpix'][0]
    deg_extent = 153.5  # Calculated by looking at the Dec. extent of the image
    dim = np.shape(jones_struct[0, 0])[0]

    output_path = '/Users/ruby/Astro/jones_matrix_plotting'

    jones_mat = np.zeros((2, 2, dim, dim), dtype=complex)
    for ind1 in range(2):
        for ind2 in range(2):
            jones_mat[ind1, ind2, :, :] = jones_struct[ind1, ind2]
    jones_mat[np.where(jones_mat == 0.)] = np.nan

    jones_rot_mat = np.zeros((2, 2, dim, dim), dtype=complex)
    for ind1 in range(2):
        for ind2 in range(2):
            jones_rot_mat[ind1, ind2, :, :] = jones_rot_struct[ind1, ind2]
    jones_rot_mat[np.where(jones_rot_mat == 0.)] = np.nan

    jones_amp = np.sqrt(
        np.real(jones_mat[:, 0, :, :])**2. + np.real(jones_mat[:, 1, :, :])**2.
    )
    k_mat = np.zeros((2, 2, dim, dim))
    for ind1 in range(2):
        for ind2 in range(2):
            k_mat[ind1, ind2, :, :] = np.real(
                jones_mat[ind1, ind2, :, :]/jones_amp[ind1, :, :]
            )
    for ind in range(2):
        jones_mat[ind, :, :, :] /= np.nanmax(jones_amp[ind, :, :])
        jones_rot_mat[ind, :, :, :] /= np.nanmax(jones_amp[ind, :, :])
        jones_amp[ind, :, :] /= np.nanmax(jones_amp[ind, :, :])

    # Plot elements of the Jones matrix
    use_cmap = matplotlib.cm.get_cmap('Spectral')
    use_cmap.set_bad(color='whitesmoke')
    plt.rcParams.update({'mathtext.default':  'regular'})
    for ind1 in range(2):
        for ind2 in range(2):
            fig, ax = plt.subplots()
            plt.imshow(
                np.real(jones_mat[ind1, ind2, :, :]), origin='lower',
                interpolation='none',
                vmin=-1, vmax=1,
                cmap=use_cmap,
                extent=[-dim/2, dim/2, -dim/2, dim/2]
            )
            plt.axis('off')  # Turn off axes
            ax = add_polar_axes(ax, dim, deg_extent)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel(
                'Beam Response, Real Part', rotation=270, labelpad=15
            )
            plt.title('$J^{{ZA}}$ [{}, {}]'.format(ind1, ind2))
            plt.savefig(
                '{}/jones_mat_{}{}.png'.format(output_path, ind1, ind2), dpi=300
            )

    # Plot elements of the rotated (RA/Dec) Jones matrix
    use_cmap = matplotlib.cm.get_cmap('Spectral')
    use_cmap.set_bad(color='whitesmoke')
    plt.rcParams.update({'mathtext.default':  'regular'})
    for ind1 in range(2):
        for ind2 in range(2):
            fig, ax = plt.subplots()
            plt.imshow(
                np.real(jones_rot_mat[ind1, ind2, :, :]), origin='lower',
                interpolation='none',
                vmin=-1, vmax=1,
                cmap=use_cmap,
                extent=[-dim/2, dim/2, -dim/2, dim/2]
            )
            plt.axis('off')  # Turn off axes
            ax = add_polar_axes(ax, dim, deg_extent)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel(
                'Beam Response, Real Part', rotation=270, labelpad=15
            )
            plt.title('$J^{{RD}}$ [{}, {}]'.format(ind1, ind2))
            plt.savefig(
                '{}/jones_rot_mat_{}{}.png'.format(output_path, ind1, ind2), dpi=300
            )

    # Plot beam amplitudes
    use_cmap = matplotlib.cm.get_cmap('viridis')
    use_cmap.set_bad(color='whitesmoke')
    for pol in range(2):
        if pol == 0:
            pol_name = 'P'
        else:
            pol_name = 'Q'
        fig, ax = plt.subplots()
        plt.imshow(
            jones_amp[pol, :, :], origin='lower', interpolation='none',
            vmin=0, vmax=1,
            cmap=use_cmap,
            extent=[-dim/2, dim/2, -dim/2, dim/2]
        )
        plt.axis('off')  # Turn off axes
        ax = add_polar_axes(ax, dim, deg_extent)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(
            'Beam Response Amplitude', rotation=270, labelpad=15
        )
        plt.title('Beam Amplitude, Polarization {}'.format(pol_name))
        plt.savefig('{}/beam_amp_{}.png'.format(output_path, pol_name), dpi=300)

    # Plot instrumental basis
    use_cmap = matplotlib.cm.get_cmap('Greys')
    colors = ['tab:blue', 'tab:red']
    line_length = 70.
    plot_vals = np.zeros((dim, dim))
    fig, ax = plt.subplots()
    # Plot zeroed data to set plot extent
    plt.imshow(
        plot_vals, origin='lower', interpolation='none',
        vmin=0, vmax=1,
        cmap=use_cmap,
        extent=[-dim/2, dim/2, -dim/2, dim/2]
    )
    plt.axis('off')  # Turn off axes
    x_pixel_vals = np.linspace(-dim/2., dim/2., dim)
    y_pixel_vals = np.linspace(-dim/2., dim/2., dim)
    use_x_pixels = np.rint(np.linspace(0, dim-1, 20)).astype(int)
    use_y_pixels = np.rint(np.linspace(0, dim-1, 20)).astype(int)
    for xpix in use_x_pixels:
        for ypix in use_y_pixels:
            if (np.isfinite(np.mean(jones_amp[:, xpix, ypix]))):
                par_ang = np.arctan2(
                    y_pixel_vals[ypix], x_pixel_vals[xpix]
                )
                for pol in range(2):
                    plot_angle = np.arctan2(
                        k_mat[pol, 1, xpix, ypix], -k_mat[pol, 0, xpix, ypix]
                    ) + par_ang
                    plt.plot([
                        x_pixel_vals[xpix] - line_length/2.*np.cos(plot_angle),
                        x_pixel_vals[xpix] + line_length/2.*np.cos(plot_angle)
                    ], [
                        y_pixel_vals[ypix] - line_length/2.*np.sin(plot_angle),
                        y_pixel_vals[ypix] + line_length/2.*np.sin(plot_angle)
                    ], color=colors[pol], linewidth=1.2)
    ax = add_polar_axes(ax, dim, deg_extent)

    plt.title('Instrumental Basis')
    plt.savefig('{}/instr_basis.png'.format(output_path), dpi=300)


def add_polar_axes(ax, dim, deg_extent):

    # Plot azimuth lines
    azimuth_lines = 8
    annotate_pad_from_edge = 50
    az_lines_angles = np.linspace(0, 360, azimuth_lines, endpoint=False)
    for line_ang in az_lines_angles:
        end_coords = [
            np.cos(np.radians(line_ang)),
            np.sin(np.radians(line_ang))
        ]
        end_coords = [
            end_coords[ind]*dim/(2*np.max(np.abs(end_coords)))
            for ind in range(len(end_coords))
        ]
        if line_ang < 180:
            plt.plot(
                [-end_coords[0], end_coords[0]],
                [-end_coords[1], end_coords[1]],
                color='grey', linewidth=0.5
            )
        annotation_loc = end_coords
        for ind in range(2):
            if annotation_loc[ind] > dim/2-annotate_pad_from_edge:
                annotation_loc[ind] = dim/2-annotate_pad_from_edge
            elif annotation_loc[ind] < -(dim/2-annotate_pad_from_edge):
                annotation_loc[ind] = -(dim/2-annotate_pad_from_edge)
        ax.annotate(
            '${}^\circ$'.format(round(line_ang)),
            (annotation_loc[0], annotation_loc[1]),
            horizontalalignment='center', verticalalignment='center'
        )

    # Plot ZA contours
    plot_za_values = [30., 45., 60., 75.]
    use_rad = dim/2./np.sin(np.radians(deg_extent/2.))
    for za in plot_za_values:
        circ_radius = use_rad*np.sin(np.radians(za))
        circle = plt.Circle(
            (0, 0), circ_radius,
            edgecolor='grey', facecolor='none', linewidth=0.5
        )
        ax.add_patch(circle)
        annotate_ang = 30
        annotate_pad = 70
        annotation_loc = (circ_radius+annotate_pad)*np.array([
            np.cos(np.radians(annotate_ang)), np.sin(np.radians(annotate_ang))
        ])
        ax.annotate(
            '${}^\circ$'.format(round(za)),
            (annotation_loc[0], annotation_loc[1]),
            horizontalalignment='center', verticalalignment='center'
        )

    return ax


if __name__ == '__main__':
    main()
