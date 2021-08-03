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
    obs_path = '/Users/ruby/Astro/jones_matrix_plotting/1131454296_obs.sav'
    obs_struct = scipy.io.readsav(obs_path)['obs']
    degpix = obs_struct['degpix'][0]
    deg_extent = 153.5  # Calculated by looking at the Dec. extent of the image
    dim = np.shape(jones_struct[0, 0])[0]

    jones_p_amp = np.sqrt(
        np.real(jones_struct[0, 0])**2. + np.real(jones_struct[0, 1])**2.
    )
    jones_p_amp[np.where(jones_p_amp == 0.)] = np.nan
    jones_q_amp = np.sqrt(
        np.real(jones_struct[1, 0])**2. + np.real(jones_struct[1, 1])**2.
    )
    jones_q_amp[np.where(jones_q_amp == 0.)] = np.nan
    k_mat = np.zeros((2, 2, dim, dim))
    for ind in range(2):
        k_mat[0, ind, :, :] = np.real((jones_struct[0, ind])/jones_p_amp)
        k_mat[1, ind, :, :] = np.real((jones_struct[1, ind])/jones_q_amp)
    jones_p_amp /= np.nanmax(jones_p_amp)
    jones_q_amp /= np.nanmax(jones_q_amp)

    angle_p = np.squeeze(np.arctan2(k_mat[0, 1, :, :], k_mat[0, 0, :, :]))
    angle_q = np.squeeze(np.arctan2(k_mat[1, 1, :, :], k_mat[1, 0, :, :]))

    plot_amp = False
    if plot_amp:
        plot_vals = jones_p_amp
        use_cmap = matplotlib.cm.get_cmap('viridis')
        use_cmap.set_bad(color='grey')
    else:
        plot_vals = np.zeros_like(jones_p_amp)
        use_cmap = matplotlib.cm.get_cmap('Greys')

    fig, ax = plt.subplots()
    plt.imshow(
        plot_vals, origin='lower', interpolation='none', vmin=0, vmax=1,
        cmap=use_cmap,
        extent=[-dim/2, dim/2, -dim/2, dim/2]
    )
    plt.axis('off')  # Turn off axes
    if plot_amp:
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(
            'Beam Response, Abs. Value', rotation=270, labelpad=15
        )

    if not plot_amp:
        x_pixel_vals = np.linspace(-dim/2., dim/2., dim)
        y_pixel_vals = np.linspace(-dim/2., dim/2., dim)
        use_x_pixels = np.rint(np.linspace(0, dim-1, 20)).astype(int)
        use_y_pixels = np.rint(np.linspace(0, dim-1, 20)).astype(int)

        for xpix in use_x_pixels:
            for ypix in use_y_pixels:
                if (
                    np.isfinite(jones_p_amp[xpix, ypix])
                    and np.isfinite(jones_q_amp[xpix, ypix])
                ):
                    par_ang = np.arctan2(
                        y_pixel_vals[ypix], x_pixel_vals[xpix]
                    )
                    for pol in range(2):
                        plot_angle = np.arctan2(k_mat[pol, 1, xpix, ypix], -k_mat[pol, 0, xpix, ypix])+par_ang

                        line_length = 50.
                        plt.plot(
                            [
                                x_pixel_vals[xpix] - line_length/2.*np.cos(plot_angle),
                                x_pixel_vals[xpix] + line_length/2.*np.cos(plot_angle)
                            ],
                            [
                                y_pixel_vals[ypix] - line_length/2.*np.sin(plot_angle),
                                y_pixel_vals[ypix] + line_length/2.*np.sin(plot_angle)
                            ],
                            color='black', linewidth=0.5
                        )

    # Plot azimuth lines
    azimuth_lines = 8
    annotate_pad_from_edge = 100
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
                color='black', linewidth=0.5
            )
        annotation_loc = end_coords
        for ind in range(2):
            if annotation_loc[ind] > dim/2-annotate_pad_from_edge:
                annotation_loc[ind] = dim/2-annotate_pad_from_edge
            elif annotation_loc[ind] < -(dim/2-annotate_pad_from_edge):
                annotation_loc[ind] = -(dim/2-annotate_pad_from_edge)
        plt.annotate(
            '${}^\circ$'.format(round(line_ang)),
            (annotation_loc[0], annotation_loc[1]),
            horizontalalignment='center', verticalalignment='center'
        )

    # Plot ZA contours
    plot_za_values = [30., 50., 70.]
    use_rad = dim/2./np.tan(np.radians(deg_extent/2.))
    for za in plot_za_values:
        circ_radius = use_rad*np.tan(np.radians(za))
        circle = plt.Circle(
            (0, 0), circ_radius,
            edgecolor='black', facecolor='none', linewidth=0.5
        )
        ax.add_patch(circle)
        annotate_ang = 30
        annotate_pad = 70
        annotation_loc = (circ_radius+annotate_pad)*np.array([
            np.cos(np.radians(annotate_ang)), np.sin(np.radians(annotate_ang))
        ])
        plt.annotate(
            '${}^\circ$'.format(round(za)),
            (annotation_loc[0], annotation_loc[1]),
            horizontalalignment='center', verticalalignment='center'
        )

    plt.show()


if __name__ == '__main__':
    main()
