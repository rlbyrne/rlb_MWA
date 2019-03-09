#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os


def write_params_to_csv(
    filepath, param_name, obs_list, param_vals, overwrite=False
):

    # Get header
    file = open(filepath, 'r')
    header = file.readline().strip().split(',')
    file.close()
    params_old = np.genfromtxt(filepath, delimiter=',', skip_header=1)

    if param_name in header:
        params = params_old
    else:
        params = np.full(
            (np.shape(params_old)[0], np.shape(params_old)[1]+1), np.nan
        )
        params[
            0:np.shape(params_old)[0], 0:np.shape(params_old)[1]
        ] = params_old
        header.append(param_name)

    obs_list_old = list(params[:, header.index('obsid')])
    for obs_ind, obsid in enumerate(obs_list):
        obsid = int(obsid)
        if obsid in obs_list_old:
            if overwrite or np.isnan(params[
                obs_list_old.index(obsid), header.index(param_name)
            ]):
                params[
                    obs_list_old.index(obsid), header.index(param_name)
                ] = param_vals[obs_ind]
        else:
            params_new = np.full(
                (np.shape(params)[0]+1, np.shape(params)[1]), np.nan
            )
            params_new[0:np.shape(params)[0], 0:np.shape(params)[1]] = params
            params_new[-1, header.index('obsid')] = obsid
            params_new[-1, header.index(param_name)] = param_vals[obs_ind]
            params = params_new

    np.savetxt(
        filepath, params, delimiter=',', header=','.join(header), comments=''
    )


def qa_param_histogram(
    filepath, param_name, saveloc='', logx=False, logy=False
):

    # Get header
    file = open(filepath, 'r')
    header = file.readline().strip().split(',')
    file.close()

    params = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    use_params = params[:, header.index(param_name)]
    use_params = use_params[np.where(~np.isnan(use_params))[0]]

    # Histogram data
    if logx:
        hist, bin_edges = np.histogram(np.log(use_params), bins='auto')
        bin_edges = np.exp(bin_edges)
    else:
        hist, bin_edges = np.histogram(use_params, bins='auto')
    xvals = [bin_edges[0]]
    yvals = [0.]
    for i in range(len(bin_edges)-1):
        xvals.extend([bin_edges[i], bin_edges[i+1]])
        yvals.extend([hist[i], hist[i]])
    xvals.append(bin_edges[-1])
    yvals.append(0.)

    plt.figure()
    if logx:
        plt.xscale('log')
    plt.plot(xvals, yvals, linestyle='-', color='blue', linewidth=0.5,
             label='GLEAM flux densities')
    if logy:
        plt.yscale('log')
    else:  # Fill gets weird with a log y-axis
        plt.fill_between(xvals, 0, yvals, facecolor='blue', alpha=0.1)
    plt.xlabel(param_name)
    plt.ylabel('Histogram Count')
    if saveloc == '':
        plt.show()
    else:
        print 'Saving plot to {}'.format(saveloc)
        plt.savefig(saveloc, dpi=300)
        plt.close()


def qa_param_hist2d(
    filepath, x_param_name, y_param_name, saveloc='', log_colorbar=False,
    colorbar_range=None
):

    # Get header
    file = open(filepath, 'r')
    header = file.readline().strip().split(',')
    file.close()

    params = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    use_params_x = params[:, header.index(x_param_name)]
    use_params_y = params[:, header.index(y_param_name)]
    use_indices = np.where(~np.isnan(use_params_x + use_params_y))[0]
    use_params_x = use_params_x[use_indices]
    use_params_y = use_params_y[use_indices]

    # Histogram data
    hist, x_edges, y_edges = np.histogram2d(
        use_params_x, use_params_y, bins=10
    )
    if colorbar_range is None:
        colorbar_range = [np.min(hist), np.max(hist)]

    plt.figure()
    plt.imshow(
        hist, interpolation='none', origin='low',
        extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
        vmin=colorbar_range[0], vmax=colorbar_range[1], aspect='auto'
    )
    plt.xlabel(x_param_name)
    plt.ylabel(y_param_name)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Histogram Count', rotation=270, labelpad=15)
    if saveloc == '':
        plt.show()
    else:
        print 'Saving plot to {}'.format(saveloc)
        plt.savefig(saveloc, dpi=300)
        plt.close()


def qa_param_hist_set(filepath, param_names, saveloc='', colorbar_range=None):

    file = open(filepath, 'r')
    header = file.readline().strip().split(',')
    file.close()

    params = np.genfromtxt(filepath, delimiter=',', skip_header=1)

    fig = plt.figure(figsize=(6, 6))
    plt.rcParams.update({'font.size': 5})
    grid = plt.GridSpec(
        2*len(param_names)-1, 2*len(param_names)-1, wspace=0.1, hspace=0.1
    )

    for ind1, param1 in enumerate(param_names[:-1]):
        for ind2, param2 in enumerate(param_names[ind1+1:]):
            use_params_x = params[:, header.index(param1)]
            use_params_y = params[:, header.index(param2)]
            xrange = [np.min(use_params_x), np.max(use_params_x)]
            yrange = [np.min(use_params_y), np.max(use_params_y)]
            print grid[2*ind1+1:2*ind1+3, 2*(ind2+ind1):2*(ind2+ind1)+2]
            subfig = fig.add_subplot(
                grid[2*(ind2+ind1):2*(ind2+ind1)+2, 2*ind1+1:2*ind1+3]
            )
            subfig.plot(
                use_params_x, use_params_y, 'ok', markersize=2, alpha=0.2
            )
            subfig.axes.get_xaxis().set_visible(False)
            subfig.axes.get_yaxis().set_visible(False)

            if False:
                hist, x_edges, y_edges = np.histogram2d(
                    use_params_x, use_params_y, bins=100
                )
                x_centers = [
                    (x_edges[i]+x_edges[i+1])/2 for i in range(len(x_edges)-1)
                ]
                y_centers = [
                    (y_edges[i]+y_edges[i+1])/2 for i in range(len(y_edges)-1)
                ]
                subfig.contour(
                    x_centers, y_centers, hist.T,
                    levels=[np.max(hist)/10000., np.max(hist)/100.],
                    colors=['skyblue', 'dodgerblue']
                )

            if ind2 == 0:
                x_hist = fig.add_subplot(grid[-1, 2*ind1+1:2*ind1+3])
                x_hist.hist(
                    use_params_x,
                    orientation='vertical', color='gray', density=True
                )
                x_hist.invert_yaxis()
                x_hist.get_shared_x_axes().join(x_hist, subfig)
                x_hist.set_xlabel(param1)
                x_hist.axes.get_yaxis().set_visible(False)
            else:
                subfig.get_shared_x_axes().join(x_hist, subfig)

            if ind1 == 0:
                y_hist = fig.add_subplot(
                    grid[2*(ind2+ind1):2*(ind2+ind1)+2, 0], sharey=subfig
                )
                y_hist.hist(
                    use_params_y,
                    orientation='horizontal', color='gray', density=True
                )
                y_hist.invert_xaxis()
                y_hist.get_shared_y_axes().join(y_hist, subfig)
                y_hist.set_ylabel(param2)
                y_hist.axes.get_xaxis().set_visible(False)
            else:
                subfig.get_shared_y_axes().join(y_hist, subfig)

    if saveloc == '':
        plt.show()
    else:
        print 'Saving plot to {}'.format(saveloc)
        plt.savefig(saveloc, dpi=300)
        plt.close()


if __name__ == '__main__':

    qa_param_hist_set('/Users/ruby/EoR/diffuse_total_occ.csv', ['obsid', 'RFI occupancy', 'obsid', 'RFI occupancy'], colorbar_range=None, saveloc='/Users/ruby/Desktop/test_params_plot.png')
