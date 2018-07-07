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
    plt.title('Deconvolved Flux Density Agreement With GLEAM')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc='upper left')
    plt.show()
    plt.close()


def test_flux_hist():

    matched_ref_catalog_fluxes = [1,1,1,1,10,100,1000,500,6,999]
    matched_decon_catalog_fluxes = [1,1,1,1,10,100,100,100,6,999]

    ref_hist, bin_edges = np.histogram([np.log(val) for val in matched_ref_catalog_fluxes], bins='auto')
    decon_hist, bin_edges = np.histogram([np.log(val) for val in matched_decon_catalog_fluxes], bins=bin_edges)
    print len(ref_hist)
    print len(decon_hist)
    bin_edges = [np.exp(val) for val in bin_edges]
    xvals = [bin_edges[0]]
    ref_yvals = [0.]
    decon_yvals = [0.]
    for i in range(len(bin_edges)-1):
        xvals.extend([bin_edges[i], bin_edges[i+1]])
        ref_yvals.extend([ref_hist[i], ref_hist[i]])
        decon_yvals.extend([decon_hist[i], decon_hist[i]])
    xvals.append(bin_edges[-1])
    ref_yvals.append(0.)
    decon_yvals.append(0.)

    plt.figure()
    plt.xscale('log')
    plt.plot(xvals, ref_yvals, linestyle='-', color='blue', linewidth=0.5, label='GLEAM flux densities')
    plt.plot(xvals, decon_yvals, linestyle='-', color='red', linewidth=0.5, label='Deconvolved flux densities')
    plt.fill_between(xvals, 0, ref_yvals, facecolor='blue', alpha=0.1)
    plt.fill_between(xvals, 0, decon_yvals, facecolor='red', alpha=0.1)
    #plt.yscale('log')
    plt.xlabel('Flux Densities (Jy)')
    plt.ylabel('Histogram Count')
    plt.title('Source Flux Density Histogram')
    plt.legend(loc='upper left')
    plt.show()
    plt.close()


def plot_pos_offsets_scatter():

    match_radius = 2
    matched_ref_catalog_ras = [1,2,3,4]
    matched_ref_catalog_decs = [2,2,2,2]
    matched_decon_catalog_ras = [1,1,4,5]
    matched_decon_catalog_decs = [1,3,1,3]

    ra_offsets = [(matched_decon_catalog_ras[i]-matched_ref_catalog_ras[i])*60. for i in range(len(matched_ref_catalog_ras))]
    dec_offsets = [(matched_decon_catalog_decs[i]-matched_ref_catalog_decs[i])*60. for i in range(len(matched_ref_catalog_decs))]

    plt.figure()
    plt.grid(True, zorder=0)
    plt.scatter(ra_offsets, dec_offsets, marker='o', s=1., color='black', zorder=10)
    plt.xlabel('RA offset (arcmin)')
    plt.ylabel('Dec offset (arcmin)')
    plt.title('Source Positional Offsets, Deconvolved - GLEAM')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([-match_radius*60., match_radius*60., -match_radius*60., match_radius*60.])
    plt.show()
    plt.close()


def plot_pos_offsets_vectors():

    match_radius = 2
    search_radius = 12.
    arrow_scale_factor = 2.
    matched_ref_catalog_ras = [1,2,3,4]
    matched_ref_catalog_decs = [2,2,2,2]
    matched_decon_catalog_ras = [1,1,4,5]
    matched_decon_catalog_decs = [1,3,1,3]

    ra_offsets = [(
        matched_decon_catalog_ras[i]-matched_ref_catalog_ras[i]
    )*60. for i in range(len(matched_ref_catalog_ras))]
    dec_offsets = [(
        matched_decon_catalog_decs[i]-matched_ref_catalog_decs[i]
    )*60. for i in range(len(matched_ref_catalog_decs))]


    ave_ra = np.mean(matched_ref_catalog_ras)
    ave_dec = np.mean(matched_ref_catalog_decs)
    plt.figure()
    for i in range(len(ra_offsets)):
        if matched_ref_catalog_ras[i]-ave_ra > search_radius:
            use_ra = matched_ref_catalog_ras[i]-360.
        elif matched_ref_catalog_ras[i]-ave_ra < -search_radius:
            use_ra = matched_ref_catalog_ras[i]+360.
        else:
            use_ra = matched_ref_catalog_ras[i]
        plt.arrow(use_ra,
                  matched_ref_catalog_decs[i],
                  ra_offsets[i]*arrow_scale_factor,
                  dec_offsets[i]*arrow_scale_factor,
                  length_includes_head=True, width=.00001, head_width = .1,
                  fc='black')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title('Source Positional Offsets')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([ave_ra-search_radius-2., ave_ra+search_radius+2.,
              ave_dec-search_radius-2., ave_dec+search_radius+2.])
    plt.show()
    plt.close()



if __name__=='__main__':
  plot_pos_offsets_vectors()
