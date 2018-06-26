#!/usr/bin/python

import csv
import sys
sys.path.append('/Users/ruby/EoR/rlb_MWA/diffuse_survey_coverage_visualization')
import surveyview
import re
import matplotlib.pyplot as plt
import numpy as np


def plot_leakage_fits():

    observations = surveyview.load_survey('/Users/ruby/EoR/sidelobe_survey_obsinfo.txt')
    observations = surveyview.get_pointings(observations)

    pol_leakage_params_file = open('/Users/ruby/EoR/pol_leakage_testing/plotting/pol_leakage_fit_params.csv', 'r')
    read_csv = csv.reader(pol_leakage_params_file)
    pol_leakage_params = [line for line in read_csv]
    pol_leakage_params_file.close()
    pol_leakage_params = pol_leakage_params[1:]  # remove header
    obsids = [obs_params[0] for obs_params in pol_leakage_params]
    observations = [obs for obs in observations if obs.obsid in obsids]

    xvals = list(range(512, 1536, 8))
    yvals = list(range(512, 1536, 8))

    for obs in observations:
        params = (pol_leakage_params[obsids.index(obs.obsid)])[1:]
        params = [float(param) for param in params]
        leakage_surf_q = np.zeros((len(xvals), len(yvals)))
        leakage_surf_u = np.zeros((len(xvals), len(yvals)))
        for i, xval in enumerate(xvals):
            for j, yval in enumerate(yvals):
                leakage_surf_q[i, j] = (
                    params[0]*xval**2 + params[1]*yval**2
                    + params[2]*xval*yval + params[3]*xval + params[4]*yval
                    + params[5]
                    )
                leakage_surf_u[i, j] = (
                    params[6]*xval**2 + params[7]*yval**2
                    + params[8]*xval*yval + params[9]*xval + params[10]*yval
                    + params[11]
                    )
        pointing = obs.pointing
        pointing = re.split('\(|\)|, ||', pointing)
        pointing = list(filter(None, pointing))

        fig, ax = plt.subplots()
        plt.imshow(leakage_surf_q, origin='lower', interpolation='none',
                   #cmap='Greys_r',
                   extent=[xvals[0], xvals[-1], yvals[0], yvals[-1]],
                   vmin=-.3, vmax=.3, aspect='auto'
                   )
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
        plt.title('{} Q Leakage Fit (Pointing {}, {})'.format(obs.obsid, pointing[0], pointing[1]))
        plt.axis('equal')
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Fractional Leakage', rotation=270)  # label colorbar
        plt.savefig('/Users/ruby/EoR/pol_leakage_testing/plotting/pointing_{}_{}/{}_q_leakage.png'.format(pointing[0], pointing[1], obs.obsid), format='png', dpi=500)
        plt.close()

        fig, ax = plt.subplots()
        plt.imshow(leakage_surf_u, origin='lower', interpolation='none',
                   #cmap='Greys_r',
                   extent=[xvals[0], xvals[-1], yvals[0], yvals[-1]],
                   vmin=-.3, vmax=.3, aspect='auto'
                   )
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
        plt.title('{} U Leakage Fit (Pointing {}, {})'.format(obs.obsid, pointing[0], pointing[1]))
        plt.axis('equal')
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Fractional Leakage', rotation=270)  # label colorbar
        plt.savefig('/Users/ruby/EoR/pol_leakage_testing/plotting/pointing_{}_{}/{}_u_leakage.png'.format(pointing[0], pointing[1], obs.obsid), format='png', dpi=500)
        plt.close()


if __name__ == '__main__':
    plot_leakage_fits()
