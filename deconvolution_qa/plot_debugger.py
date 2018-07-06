#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np



def test_flux_scatter():

    matched_ref_catalog_fluxes = [1,1,1,1,10,100,1000,500,6,999]
    matched_decon_catalog_fluxes = [i*.9 for i in matched_ref_catalog_fluxes]
    fit_param = sum(
        [(
          matched_decon_catalog_fluxes[i]/matched_ref_catalog_fluxes[i]
          )**2 for i in range(len(matched_ref_catalog_fluxes))]
    )/sum(
        [(
          matched_decon_catalog_fluxes[i]/matched_ref_catalog_fluxes[i]
          ) for i in range(len(matched_ref_catalog_fluxes))]
    )
    goodness_of_fit = sum(
        [(
          matched_ref_catalog_fluxes[i]-matched_decon_catalog_fluxes[i]/fit_param
          )**2/(matched_ref_catalog_fluxes[i]**2) for i in range(len(matched_ref_catalog_fluxes))]
    )/len(matched_ref_catalog_fluxes)

    print fit_param
    print goodness_of_fit

    plt.figure()
    # plot the 1-to-1 line
    plt.plot([min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
    [min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
        linestyle=':', linewidth=0.5, color='blue', label='1-to-1 line')
    # plot the power law fit
    plt.plot([min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
    [fit_param*min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), fit_param*max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
        linestyle=':', linewidth=0.5, color='red', label='power law fit: A={:.3f}, X^2={:.2f}'.format(fit_param, goodness_of_fit))
    plt.scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes, marker='o', s=1., color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('GLEAM flux density (Jy)')
    plt.ylabel('Deconvolved flux density (Jy)')
    plt.title('Deconvolved Flux Agreement With GLEAM')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc='upper left')
    plt.show()
    plt.close()


def test_flux_hist():

    matched_ref_catalog_fluxes = [1,1,1,1,10,100,1000,500,6,999]
    matched_decon_catalog_fluxes = [i*.9 for i in matched_ref_catalog_fluxes]

    ref_hist, bin_edges = np.histogram([np.log(val) for val in matched_ref_catalog_fluxes], bins='auto')
    decon_hist, bin_edges = np.histogram([np.log(val) for val in matched_decon_catalog_fluxes], bins=bin_edges)
    print len(ref_hist)
    print len(decon_hist)
    bin_edges = [np.exp(val) for val in bin_edges]
    xvals = []
    ref_yvals = []
    decon_yvals = []
    for i in range(len(bin_edges)-1):
        xvals.extend([bin_edges[i], bin_edges[i+1]])
        ref_yvals.extend([ref_hist[i], ref_hist[i]])
        decon_yvals.extend([decon_hist[i], decon_hist[i]])

    plt.figure()
    #plt.plot(ref_xvals, ref_yvals, linestyle='-')
    plt.xscale('log')
    #plt.fill_between(xvals, 0, ref_yvals, facecolor='blue', alpha=0.3)
    plt.fill_between(xvals, 0, decon_yvals, facecolor='red', alpha=0.3)
    #plt.yscale('log')
    #plt.xlabel('GLEAM flux density (Jy)')
    #plt.ylabel('Deconvolved flux density (Jy)')
    #plt.title('Deconvolved Flux Agreement With GLEAM')
    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.legend(loc='upper left')
    plt.show()
    plt.close()

if __name__=='__main__':
  test_flux_hist()
