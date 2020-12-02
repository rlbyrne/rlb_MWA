#!/usr/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import scipy
import scipy.optimize
import scipy.stats
sys.path.append('/Users/ruby/EoR/pyuvdata')
import pyuvdata


def calc_negloglikelihood_sky_cal(
    gains, data_visibilities, model_visibilities,
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev
):

    model_visibilities_expanded = np.matmul(a_mat, model_visibilities)
    gains_expanded = (
        np.matmul(gains_exp_mat_1, gains)
        * np.matmul(gains_exp_mat_2, np.conj(gains))
    )
    prob = np.sum(np.abs(
        data_visibilities - gains_expanded*model_visibilities_expanded
    )**2)

    return (prob/data_stddev**2.)


def calc_negloglikelihood_unified_cal(
    gains, fitted_visibilities, data_visibilities, model_visibilities,
    baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev, model_stddev
):

    vis_diff = fitted_visibilities-model_visibilities
    prior = np.abs(
        np.dot(np.matmul(np.conj(vis_diff), baseline_cov_inv), vis_diff)
    )

    term1 = calc_negloglikelihood_sky_cal(
        gains, data_visibilities, fitted_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )

    return (term1 + prior/model_stddev**2.)


def calc_negloglikelihood_relative_cal(
    gains, fitted_visibilities, data_visibilities,
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev
):

    return calc_negloglikelihood_sky_cal(
        gains, data_visibilities, fitted_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )


def calc_gains_grad(
    gains, fitted_visibilities, data_visibilities,
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev
):

    gains1_expanded = np.matmul(gains_exp_mat_1, gains)
    gains2_expanded = np.matmul(gains_exp_mat_2, gains)
    vis_expanded = np.matmul(a_mat, fitted_visibilities)

    gains_grad_term_1 = (
        np.abs(gains1_expanded*vis_expanded)**2.*gains2_expanded
        - np.conj(data_visibilities)*gains1_expanded*vis_expanded
    )
    gains_grad_term_2 = (
        np.abs(np.conj(gains2_expanded)*vis_expanded)**2.*gains1_expanded
        - data_visibilities*gains2_expanded*np.conj(vis_expanded)
    )
    gains_grad = (2./data_stddev**2.)*(
        np.matmul(gains_exp_mat_2.T, gains_grad_term_1)
        + np.matmul(gains_exp_mat_1.T, gains_grad_term_2)
    )

    return gains_grad


def calc_vis_grad_unified_cal(
    gains, fitted_visibilities, data_visibilities, model_visibilities,
    baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev, model_stddev
):

    vis_diff = fitted_visibilities-model_visibilities
    gains1_expanded = np.matmul(gains_exp_mat_1, gains)
    gains2_expanded = np.matmul(gains_exp_mat_2, gains)
    gains_expanded = gains1_expanded*np.conj(gains2_expanded)

    vis_grad_term_1 = (2./data_stddev**2.) * (
        np.matmul(np.abs(gains_expanded)**2, a_mat**2.)*fitted_visibilities
    )
    vis_grad_term_2 = (-2./data_stddev**2.) * (
        np.matmul(data_visibilities*np.conj(gains_expanded), a_mat)
    )
    vis_grad_term_3 = (
        (2./model_stddev**2.) * np.matmul(baseline_cov_inv, vis_diff)
    )
    vis_grad = vis_grad_term_1 + vis_grad_term_2 + vis_grad_term_3

    return vis_grad


def calc_vis_grad_relative_cal(
    gains, fitted_visibilities, data_visibilities,
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    data_stddev
):

    gains1_expanded = np.matmul(gains_exp_mat_1, gains)
    gains2_expanded = np.matmul(gains_exp_mat_2, gains)
    gains_expanded = gains1_expanded*np.conj(gains2_expanded)

    vis_grad_term_1 = (2./data_stddev**2.) * (
        np.matmul(np.abs(gains_expanded)**2, a_mat**2.)*fitted_visibilities
    )
    vis_grad_term_2 = (-2./data_stddev**2.) * (
        np.matmul(data_visibilities*np.conj(gains_expanded), a_mat)
    )
    vis_grad = vis_grad_term_1 + vis_grad_term_2

    return vis_grad


def cost_function_unified_cal(
    x,
    N_red_baselines,
    N_ants,
    baseline_cov_inv,
    model_visibilities,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev,
    model_stddev
):

    fitted_visibilities = (
        x[-2*N_red_baselines:-N_red_baselines] + 1j*x[-N_red_baselines:]
    )
    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    cost = calc_negloglikelihood_unified_cal(
        gains, fitted_visibilities, data_visibilities, model_visibilities,
        baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev, model_stddev
    )
    return cost


def jac_function_unified_cal(
    x,
    N_red_baselines,
    N_ants,
    baseline_cov_inv,
    model_visibilities,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev,
    model_stddev
):

    fitted_visibilities = (
        x[-2*N_red_baselines:-N_red_baselines] + 1j*x[-N_red_baselines:]
    )
    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    gains_grad = calc_gains_grad(
        gains, fitted_visibilities, data_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )
    vis_grad = calc_vis_grad_unified_cal(
        gains, fitted_visibilities, data_visibilities, model_visibilities,
        baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev, model_stddev
    )

    grads = np.zeros(N_ants*2+N_red_baselines*2)
    grads[:N_ants] = np.real(gains_grad)
    grads[N_ants:2*N_ants] = np.imag(gains_grad)
    grads[-2*N_red_baselines:-N_red_baselines] = np.real(vis_grad)
    grads[-N_red_baselines:] = np.imag(vis_grad)
    return grads


def cost_function_sky_cal(
    x,
    N_ants,
    model_visibilities,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev
):

    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    cost = calc_negloglikelihood_sky_cal(
        gains, data_visibilities, model_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )
    return cost


def jac_function_sky_cal(
    x,
    N_ants,
    model_visibilities,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev
):

    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    gains_grad = calc_gains_grad(
        gains, model_visibilities, data_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )

    grads = np.zeros(N_ants*2)
    grads[:N_ants] = np.real(gains_grad)
    grads[N_ants:2*N_ants] = np.imag(gains_grad)
    return grads


def cost_function_relative_cal(
    x,
    N_red_baselines,
    N_ants,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev
):

    fitted_visibilities = (
        x[-2*N_red_baselines:-N_red_baselines] + 1j*x[-N_red_baselines:]
    )
    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    cost = calc_negloglikelihood_relative_cal(
        gains, fitted_visibilities, data_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )
    return cost


def jac_function_relative_cal(
    x,
    N_red_baselines,
    N_ants,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev
):

    fitted_visibilities = (
        x[-2*N_red_baselines:-N_red_baselines] + 1j*x[-N_red_baselines:]
    )
    gains = x[:N_ants]+1j*x[N_ants:2*N_ants]

    gains_grad = calc_gains_grad(
        gains, fitted_visibilities, data_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )
    vis_grad = calc_vis_grad_relative_cal(
        gains, fitted_visibilities, data_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )

    grads = np.zeros(N_ants*2+N_red_baselines*2)
    grads[:N_ants] = np.real(gains_grad)
    grads[N_ants:2*N_ants] = np.imag(gains_grad)
    grads[-2*N_red_baselines:-N_red_baselines] = np.real(vis_grad)
    grads[-N_red_baselines:] = np.imag(vis_grad)
    return grads


def cost_function_absolute_cal(
    x,
    gains_fit,
    ant_locs,
    model_visibilities,
    a_mat,
    gains_exp_mat_1,
    gains_exp_mat_2,
    data_visibilities,
    data_stddev
):

    amp = x[0]
    phase_grad_1 = x[1]
    phase_grad_2 = x[2]

    phase_offsets = phase_grad_1*ant_locs[:, 0] + phase_grad_2*ant_locs[:, 1]
    gains = amp*(np.cos(phase_offsets) + 1j*np.sin(phase_offsets))*gains_fit

    cost = calc_negloglikelihood_sky_cal(
        gains, data_visibilities, model_visibilities,
        a_mat, gains_exp_mat_1, gains_exp_mat_2,
        data_stddev
    )
    return cost


def optimize_unified_cal_simple(
    data_visibilities, model_visibilities,
    baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
    N_red_baselines, N_ants, data_stddev, model_stddev,
    ant_locs,  # Not used
    gains_init=None, fitted_visibilities_init=None, quiet=True
):

    method = 'Powell'
    maxiter = 100000
    xtol = 1e-20
    ftol = 1e-20

    if gains_init is None:  # Initialize the gains to 1
        gains_init = np.full(N_ants, 1.+0.j)
    # Initialize the fitted visibilities to the model visibilities
    if fitted_visibilities_init is None:
        fitted_visibilities_init = model_visibilities
    # Expand the initialized values
    x0 = np.concatenate((
        np.real(gains_init), np.imag(gains_init),
        np.real(fitted_visibilities_init), np.imag(fitted_visibilities_init)
    ))

    # Minimize the cost function
    result = scipy.optimize.minimize(
        cost_function_unified_cal, x0,
        args=(
            N_red_baselines, N_ants, baseline_cov_inv,
            model_visibilities, a_mat, gains_exp_mat_1, gains_exp_mat_2,
            data_visibilities, data_stddev, model_stddev
        ),
        method=method, options={'xtol': xtol, 'ftol': ftol, 'maxiter': maxiter}
    )
    if not quiet:
        print(result.message)

    gains_fit = result.x[:N_ants]+1j*result.x[N_ants:2*N_ants]
    vis_fit = (
        result.x[-2*N_red_baselines:-N_red_baselines]
        + 1j*result.x[-N_red_baselines:]
    )
    # Ensure that the angle of the gains is mean-zero
    avg_angle = np.arctan2(
        np.mean(np.sin(np.angle(gains_fit))),
        np.mean(np.cos(np.angle(gains_fit)))
    )
    gains_fit *= np.cos(avg_angle) - 1j*np.sin(avg_angle)

    return gains_fit, vis_fit


def optimize_unified_cal(
    data_visibilities, model_visibilities,
    baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
    N_red_baselines, N_ants, data_stddev, model_stddev,
    ant_locs,  # Not used
    gains_init=None, fitted_visibilities_init=None, quiet=True
):

    method = 'CG'
    maxiter = 100000

    if gains_init is None:  # Initialize the gains to 1
        gains_init = np.full(N_ants, 1.+0.j)
    # Initialize the fitted visibilities to the model visibilities
    if fitted_visibilities_init is None:
        fitted_visibilities_init = model_visibilities
    # Expand the initialized values
    x0 = np.concatenate((
        np.real(gains_init), np.imag(gains_init),
        np.real(fitted_visibilities_init), np.imag(fitted_visibilities_init)
    ))

    # Minimize the cost function
    result = scipy.optimize.minimize(
        cost_function_unified_cal, x0, jac=jac_function_unified_cal,
        args=(
            N_red_baselines, N_ants, baseline_cov_inv,
            model_visibilities, a_mat, gains_exp_mat_1, gains_exp_mat_2,
            data_visibilities, data_stddev, model_stddev
        ),
        method=method, options={'disp': True, 'maxiter': maxiter}
    )
    if not quiet:
        print(result.message)

    gains_fit = result.x[:N_ants]+1j*result.x[N_ants:2*N_ants]
    vis_fit = (
        result.x[-2*N_red_baselines:-N_red_baselines]
        + 1j*result.x[-N_red_baselines:]
    )
    # Ensure that the angle of the gains is mean-zero
    avg_angle = np.arctan2(
        np.mean(np.sin(np.angle(gains_fit))),
        np.mean(np.cos(np.angle(gains_fit)))
    )
    gains_fit *= np.cos(avg_angle) - 1j*np.sin(avg_angle)

    return gains_fit, vis_fit


def optimize_sky_cal(
    data_visibilities, model_visibilities,
    baseline_cov_inv,  # Not used
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    N_red_baselines, N_ants, data_stddev,
    model_stddev,  # Not used
    ant_locs,  # Not used
    gains_init=None, fitted_visibilities_init=None, quiet=True
):

    method = 'CG'
    maxiter = 100000

    if gains_init is None:  # Initialize the gains to 1
        gains_init = np.full(N_ants, 1.+0.j)
    # Expand the initialized values
    x0 = np.concatenate((
        np.real(gains_init), np.imag(gains_init)
    ))

    # Minimize the cost function
    result = scipy.optimize.minimize(
        cost_function_sky_cal, x0, jac=jac_function_sky_cal,
        args=(
            N_ants, model_visibilities, a_mat,
            gains_exp_mat_1, gains_exp_mat_2,
            data_visibilities, data_stddev
        ),
        method=method, options={'disp': True, 'maxiter': maxiter}
    )
    if not quiet:
        print(result.message)

    gains_fit = result.x[:N_ants]+1j*result.x[N_ants:2*N_ants]
    # Ensure that the angle of the gains is mean-zero
    avg_angle = np.arctan2(
        np.mean(np.sin(np.angle(gains_fit))),
        np.mean(np.cos(np.angle(gains_fit)))
    )
    gains_fit *= np.cos(avg_angle) - 1j*np.sin(avg_angle)
    vis_fit = model_visibilities

    return gains_fit, vis_fit


def optimize_redundant_cal(
    data_visibilities, model_visibilities,
    baseline_cov_inv,  # Not used
    a_mat, gains_exp_mat_1, gains_exp_mat_2,
    N_red_baselines, N_ants, data_stddev,
    model_stddev,  # Not used
    ant_locs,
    gains_init=None, fitted_visibilities_init=None, quiet=True
):

    method = 'CG'
    maxiter = 100000

    if gains_init is None:  # Initialize the gains to 1
        gains_init = np.full(N_ants, 1.+0.j)
    # Initialize the fitted visibilities to the model visibilities
    if fitted_visibilities_init is None:
        fitted_visibilities_init = model_visibilities
    # Expand the initialized values
    x0 = np.concatenate((
        np.real(gains_init), np.imag(gains_init),
        np.real(fitted_visibilities_init), np.imag(fitted_visibilities_init)
    ))

    # Relative calibration
    result = scipy.optimize.minimize(
        cost_function_relative_cal, x0, jac=jac_function_relative_cal,
        args=(
            N_red_baselines, N_ants,
            a_mat, gains_exp_mat_1, gains_exp_mat_2,
            data_visibilities, data_stddev
        ),
        method=method, options={'disp': True, 'maxiter': maxiter}
    )
    if not quiet:
        print(result.message)

    gains_fit = result.x[:N_ants]+1j*result.x[N_ants:2*N_ants]
    vis_fit = (
        result.x[-2*N_red_baselines:-N_red_baselines]
        + 1j*result.x[-N_red_baselines:]
    )
    # Ensure that the angle of the gains is mean-zero
    avg_angle = np.arctan2(
        np.mean(np.sin(np.angle(gains_fit))),
        np.mean(np.cos(np.angle(gains_fit)))
    )
    gains_fit *= np.cos(avg_angle) - 1j*np.sin(avg_angle)

    # Absolute calibration
    x0 = np.array([1., 0., 0.])  # Amp, phase_grad_1, phase_grad_2
    method = 'Powell'
    result = scipy.optimize.minimize(
        cost_function_absolute_cal, x0,
        args=(
            gains_fit, ant_locs,
            model_visibilities,
            a_mat, gains_exp_mat_1, gains_exp_mat_2,
            data_visibilities, data_stddev
        ),
        method=method
    )
    if not quiet:
        print(result.message)
        print(result.x)
    abs_cal_params = result.x
    phase_offset = (
        abs_cal_params[1]*ant_locs[:, 0] + abs_cal_params[2]*ant_locs[:, 1]
    )
    shift_phase = np.cos(phase_offset) + 1j*np.sin(phase_offset)
    gains_fit *= abs_cal_params[0]*shift_phase
    vis_shift_phase = np.matmul(
        np.linalg.pinv(a_mat), (
            np.matmul(gains_exp_mat_1, shift_phase)
            * np.matmul(gains_exp_mat_2, shift_phase)
        )
    )
    vis_fit /= (abs_cal_params[0]**2. * vis_shift_phase)

    return gains_fit, vis_fit


def plot_model_errors(
    model_vis, data_vis, model_stddev, data_stddev,
    savepath=None
):

    plt.figure(figsize=[10, 8])
    plt.plot(
        np.real(model_vis - data_vis).flatten(),
        np.imag(model_vis - data_vis).flatten(),
        'x', color='black', markersize=5
    )
    plt.axis('square')
    plt.xlim(-1.15, 1.15)
    plt.ylim(-1.15, 1.15)
    plt.scatter([0], [0], marker='P', color='white', edgecolors='black', s=150)
    circle = plt.Circle((0, 0), model_stddev, fill=False, color='black')
    plt.gcf().gca().add_artist(circle)
    circle = plt.Circle(
        (0, 0), data_stddev, fill=False, color='blue', linestyle='dashed'
    )
    plt.gcf().gca().add_artist(circle)
    plt.xlabel('Model Visibility Error, Real Part (Jy)')
    plt.ylabel('Model Visibility Error, Imaginary Part (Jy)')
    plt.tight_layout()  # Ensure that axes labels don't get cut off
    if savepath is None:
        plt.show()
    else:
        plt.savefig(savepath, dpi=600)
    plt.close()


def produce_kde(data, xvals, yvals):

    data_real = np.real(data).flatten()
    data_imag = np.imag(data).flatten()
    rvs = np.append(data_real[:, np.newaxis], data_imag[:, np.newaxis], axis=1)

    kde = scipy.stats.kde.gaussian_kde(rvs.T)

    # Regular grid to evaluate KDE upon
    x, y = np.meshgrid(xvals, yvals)
    grid_coords = np.append(x.reshape(-1, 1), y.reshape(-1, 1), axis=1)

    kde_vals = kde(grid_coords.T)
    kde_vals = kde_vals.reshape(len(xvals), len(yvals))

    percent_vals = np.zeros_like(kde_vals)
    kde_total = np.sum(kde_vals)
    running_total = 0.
    # Sort values in reverse order
    for val in np.sort(kde_vals.flatten())[::-1]:
        percent_vals[np.where(kde_vals == val)] = running_total/kde_total
        running_total += val

    return kde_vals, percent_vals


def histogram_plot_2d(
    vals, gains=True, plot_range=None, nbins=50, colorbar_range=None,
    plot_contours=True, savepath=None, title=None
):

    matplotlib.rcParams.update({'font.size': 20})

    # Set defaults
    if gains:
        if plot_range is None:
            plot_range = .03
        if colorbar_range is None:
            colorbar_range = [0, .005]
        axis_label = 'Fit Gain'
    else:
        if plot_range is None:
            plot_range = .3
        if colorbar_range is None:
            colorbar_range = [0, .005]
        axis_label = 'Fit Vis. Error (Jy)'

    bins_imag = np.linspace(-plot_range, plot_range, num=nbins+1)
    if gains:
        bins_real = bins_imag+1
    else:
        bins_real = bins_imag

    for data_set_ind in range(np.shape(vals)[2]):

        plot_data = vals[:, :, data_set_ind]
        hist, x_edges, y_edges = np.histogram2d(
            np.real(plot_data).flatten(), np.imag(plot_data).flatten(),
            bins=[bins_real, bins_imag]
        )
        hist /= np.sum(hist)
        if plot_contours:
            kde, percent_plot = produce_kde(plot_data, bins_real, bins_imag)

        plt.figure(figsize=[7, 5.5])
        plt.imshow(
            hist.T, interpolation='none', origin='lower',
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
            vmin=colorbar_range[0], vmax=colorbar_range[1], aspect='equal',
            cmap='inferno'
        )
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Histogram Density', rotation=270, labelpad=20)
        if plot_contours:
            plt.contour(
                bins_real, bins_imag, percent_plot, levels=[.5, .9],
                colors='white', linestyles=['solid', 'dashed'], linewidths=1
            )
        if gains:
            plt.scatter(
                [1], [0], marker='P', color='white', edgecolors='black', s=150
            )
        else:
            plt.scatter(
                [0], [0], marker='P', color='white', edgecolors='black', s=150
            )
        plt.xlabel('{}, Real Part'.format(axis_label))
        #plt.xticks(rotation=45)
        plt.ylabel('{}, Imag. Part'.format(axis_label))
        if title is not None:
            plt.title(title)
        plt.tight_layout()  # Ensure that axes labels don't get cut off
        if savepath is None:
            plt.show()
        else:
            plt.savefig(savepath, dpi=600)
        plt.close()


def calibrate(
    data_path='/Users/ruby/EoR/compact_redundant_array_sim_May2020/square_grid_sim__results.uvh5',
    model_path='/Users/ruby/EoR/compact_redundant_array_sim_May2020/square_grid_100mjy_sim_results.uvh5',
    n_trials=100,
    data_stddev=.2,
    model_stddev_scaling=None,  # Can be list
    calibration_styles=['unified'],
    use_covariances=True,
    simple_optimization=False,
    create_model_errors_plot=False,
    model_errors_plot_savepath=None,
    perfect_model=False,
    model_injected_noise_stddev_jy=None,
    quiet=False
):

    if 'unified' not in calibration_styles:
        if use_covariances:
            print(
                'WARNING: Not using unified calibration. '
                'No covariances applied.'
            )
            use_covariances = False
        if simple_optimization:
            print(
                'WARNING: Simple optimization is implemented for unified '
                'calibration only.'
            )
        if model_stddev_scaling is not None:
            print(
                'WARNING: Simple model stddev scaling is implemented for '
                'unified calibration only.'
            )
            model_stddev_scaling = None

    if perfect_model:
        if model_path != data_path:
            print(
                'WARNING: perfect_model is set. Using file path {} for '
                'model visibilities'.format(data_path)
            )
        model_path = data_path

    uvw_match_tolerance = 1e-12

    # Load data from pyuvsim simulation:
    data_sim_compact = pyuvdata.UVData()
    data_sim_compact.read_uvh5(data_path)

    # Remove autos
    data_sim_compact.select(ant_str='cross')
    # Use only XX polarizations
    data_sim_compact.select(polarizations=[-5])

    # Convert baselines to have u>0
    data_sim_compact.conjugate_bls(
        convention='u>0', use_enu=False, uvw_tol=0.01
    )

    baseline_groups, vec_bin_centers, lengths, conjugates = (
        data_sim_compact.get_redundancies(
            tol=0.1, use_antpos=False, include_conjugates=True,
            include_autos=True, conjugate_bls=False
        )
    )

    # Define constants
    N_red_baselines = np.shape(baseline_groups)[0]
    N_ants = data_sim_compact.Nants_data

    # Reorder visibilities
    data_sim_vis_no_noise = np.zeros(N_red_baselines, dtype=np.complex_)
    for red_group in range(N_red_baselines):
        found_group = False
        for red_group_2 in range(N_red_baselines):
            if np.abs(np.sum(
                data_sim_compact.uvw_array[red_group]
                - vec_bin_centers[red_group_2]
            )) < uvw_match_tolerance:
                data_sim_vis_no_noise[red_group] = (
                    data_sim_compact.data_array[red_group_2, 0, 0, 0]
                )
                found_group = True
                break
        if not found_group:
            print('ERROR: Visibility not found.')

    # Make noiseless data
    data_sim_expanded = data_sim_compact.copy()
    data_sim_expanded.inflate_by_redundancy()

    # Define constant
    N_vis = data_sim_expanded.Nbls

    # Load data with missing sources from pyuvsim simulation:
    model_sim = pyuvdata.UVData()
    model_sim.read_uvh5(model_path)

    # Remove autos
    model_sim.select(ant_str='cross')
    # Use only XX polarizations
    model_sim.select(polarizations=[-5])

    # Convert baselines to have u>0
    model_sim.conjugate_bls(convention='u>0', use_enu=False, uvw_tol=0.01)

    model_sim_visibilities = np.zeros(N_red_baselines, dtype=np.complex_)
    for red_group in range(N_red_baselines):
        found_group = False
        for red_group_2 in range(N_red_baselines):
            if np.abs(np.sum(
                model_sim.uvw_array[red_group]-vec_bin_centers[red_group_2]
            )) < uvw_match_tolerance:
                model_sim_visibilities[red_group] = (
                    model_sim.data_array[red_group_2, 0, 0, 0]
                )
                found_group = True
                break
        if not found_group:
            print('ERROR: Visibility not found.')

    # Create the baseline covariance matrix
    baseline_cov_array = np.diag(np.full(N_red_baselines, 1.))
    if use_covariances:
        min_bl_length = 14.
        tolerance = .01
        for bl_1 in range(N_red_baselines):
            for bl_2 in [ind for ind in range(N_red_baselines) if ind != bl_1]:
                bl_separation_sq = (
                    (vec_bin_centers[bl_1, 0]-vec_bin_centers[bl_2, 0])**2
                    + (vec_bin_centers[bl_1, 1]-vec_bin_centers[bl_2, 1])**2
                )
                if (
                    (min_bl_length-tolerance)**2 <= bl_separation_sq
                    <= (min_bl_length+tolerance)**2
                ):
                    baseline_cov_array[bl_1, bl_2] = 0.1617
                elif (
                    2*(min_bl_length-tolerance)**2 <= bl_separation_sq
                    <= 2*(min_bl_length+tolerance)**2
                ):
                    baseline_cov_array[bl_1, bl_2] = 0.0176
        # Invert the matrix
        baseline_cov_inv = np.linalg.inv(baseline_cov_array)
    else:  # Identity matrix
        baseline_cov_inv = baseline_cov_array

    # Create the A matrix
    a_mat = np.zeros((N_vis, N_red_baselines))
    for vis_ind in range(N_vis):
        for red_group in range(N_red_baselines):
            if np.abs(np.sum(
                data_sim_expanded.uvw_array[vis_ind]-vec_bin_centers[red_group]
            )) < uvw_match_tolerance:
                a_mat[vis_ind, red_group] = 1
                break

    # Create gains expand matrices
    gains_exp_mat_1 = np.zeros((N_vis, N_ants), dtype=np.int)
    gains_exp_mat_2 = np.zeros((N_vis, N_ants), dtype=np.int)
    for baseline in range(N_vis):
        gains_exp_mat_1[baseline, data_sim_expanded.ant_1_array[baseline]] = 1
        gains_exp_mat_2[baseline, data_sim_expanded.ant_2_array[baseline]] = 1

    if 'redundant' in calibration_styles:  # Calculate antenna positions
        ant_locs = (
            data_sim_expanded.antenna_positions
            + data_sim_expanded.telescope_location
        )
        ant_locs = pyuvdata.utils.ENU_from_ECEF(
            ant_locs, *data_sim_expanded.telescope_location_lat_lon_alt
        )
        ant_locs = ant_locs[:, :2]  # Discard altitude information
    else:
        ant_locs = None

    if (
        model_injected_noise_stddev_jy != 0
        and model_injected_noise_stddev_jy is not None
    ):
        # Inject noise on the model
        model_sim_visibilities += (
            np.random.normal(
                0, model_injected_noise_stddev_jy, N_red_baselines
            ) + 1j*np.random.normal(
                0, model_injected_noise_stddev_jy, N_red_baselines
            )
        )

    # Calculate deviation between model and true data
    model_stddev_sim = np.sqrt(np.mean(np.abs(
        model_sim_visibilities-data_sim_vis_no_noise
    )**2)/2.)

    if create_model_errors_plot:
        plot_model_errors(
            model_sim_visibilities, data_sim_vis_no_noise,
            model_stddev_sim, data_stddev,
            savepath=model_errors_plot_savepath
        )

    if model_stddev_scaling is None:
        model_stddev_scaling = [1.]
    model_stddev_use_vals = [
        scaling*model_stddev_sim for scaling in model_stddev_scaling
    ]

    # Run optimization
    gain_vals = np.zeros((
        N_ants, n_trials, len(model_stddev_use_vals),
        len(calibration_styles)
    ), dtype=np.complex_)
    vis_diff_vals = np.zeros((
        N_red_baselines, n_trials, len(model_stddev_use_vals),
        len(calibration_styles)
    ), dtype=np.complex_)

    for trial in range(n_trials):

        # Generate noisy data
        data_visibilities = (
            data_sim_expanded.data_array[:, 0, 0, 0]
            + np.random.normal(0, data_stddev, N_vis)
            + 1j*np.random.normal(0, data_stddev, N_vis)
        )

        for style_ind, cal_style in enumerate(calibration_styles):

            if cal_style == 'sky':
                optimize_function_name = optimize_sky_cal
            elif cal_style == 'unified':
                if simple_optimization:
                    optimize_function_name = optimize_unified_cal_simple
                else:
                    optimize_function_name = optimize_unified_cal
            elif cal_style == 'redundant':
                optimize_function_name = optimize_redundant_cal
            else:
                print(
                    'ERROR: calibration_style options are "sky", "unified", '
                    'and "redundant". Exiting.'
                )
                sys.exit(1)

            for stddev_ind, model_stddev_use in enumerate(
                model_stddev_use_vals
            ):
                if not quiet:
                    print('***{} Calibration, Version {}, Trial {}***'.format(
                        cal_style, stddev_ind+1, trial+1
                    ))
                gains_fit, vis_fit = optimize_function_name(
                    data_visibilities, model_sim_visibilities,
                    baseline_cov_inv, a_mat, gains_exp_mat_1, gains_exp_mat_2,
                    N_red_baselines, N_ants, data_stddev, model_stddev_use,
                    ant_locs,
                    gains_init=None,
                    #fitted_visibilities_init=np.matmul(
                    #    np.linalg.pinv(a_mat), data_visibilities
                    #),
                    fitted_visibilities_init=model_sim_visibilities,
                    quiet=quiet
                )

                gain_vals[:, trial, stddev_ind, style_ind] = gains_fit
                vis_fit_diff = vis_fit-data_sim_vis_no_noise
                vis_diff_vals[
                    :, trial, stddev_ind, style_ind
                ] = vis_fit_diff

    return gain_vals, vis_diff_vals


if __name__ == '__main__':
    gain_vals, vis_diff_vals = calibrate(
        n_trials=100
    )
    histogram_plot_2d(
        gain_vals[:, :, :, 0], gains=True,
        savepath='/Users/ruby/EoR/unified_calibration_plotting_oct2020/test_abs_cal/gains_regular_unified_cal.png'
    )
