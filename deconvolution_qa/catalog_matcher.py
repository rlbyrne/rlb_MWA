#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np


def match_catalogs():

    fhd_run_path = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018'
    obsid = '1130773024'
    ref_catalog_path = '/Users/rubybyrne/FHD/catalog_data/GLEAM_plus_rlb2017.sav'
    output_path = '/Users/rubybyrne/diffuse_survey_qa_plots'
    search_radius = 12.  # match sources out to 12 degrees from pointing center
    match_radius = .05  # a match is within 3 arcmin of the source

    obs_sav = scipy.io.readsav(
        '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
        )
    obs_ra = float(obs_sav['obs']['obsra'][0])
    obs_dec = float(obs_sav['obs']['obsdec'][0])

    decon_catalog = scipy.io.readsav(
        '{}/deconvolution/{}_fhd.sav'.format(fhd_run_path, obsid)
        )['source_array']
    # grab sources near the pointing center
    decon_catalog_limited = []
    for source in decon_catalog:
        if min((float(source['ra'])-obs_ra)**2.,
               (float(source['ra'])-obs_ra+360.)**2.,
               (float(source['ra'])-obs_ra-360.)**2.
               ) + (float(source['dec'])-obs_dec)**2. < search_radius**2.:
            decon_catalog_limited.append(source)
    # sort sources by flux
    decon_catalog_limited.sort(key=lambda x: float(x['flux']['I']), reverse=True)
    decon_catalog_limited = decon_catalog_limited[0:100]  # for debugging

    ref_catalog = scipy.io.readsav(ref_catalog_path)['catalog']
    # grab sources near the pointing center
    ref_catalog_limited = []
    for source in ref_catalog:
        if min((float(source['ra'])-obs_ra)**2.,
               (float(source['ra'])-obs_ra+360.)**2.,
               (float(source['ra'])-obs_ra-360.)**2.
               ) + (float(source['dec'])-obs_dec)**2. < search_radius**2.:
            ref_catalog_limited.append(source)
    # sort sources by flux
    ref_catalog_limited.sort(key=lambda x: float(x['flux']['I']), reverse=True)
    ref_catalog_limited = ref_catalog_limited[0:100]  # for debugging

    plot_flux_hist([source['flux']['I'] for source in ref_catalog_limited],
                   [source['flux']['I'] for source in decon_catalog_limited],
                   '{}/{}_flux_hist.png'.format(output_path, obsid))

    ref_catalog_ras = [float(source['ra']) for source in ref_catalog_limited]
    ref_catalog_decs = [float(source['dec']) for source in ref_catalog_limited]

    # match sources
    match_index_list = [-1.]*len(decon_catalog_limited)
    ref_catalog_indices = list(range(len(ref_catalog_limited)))
    for i, decon_source in enumerate(decon_catalog_limited):
        ra = float(decon_source['ra'])
        dec = float(decon_source['dec'])
        for ref_index in ref_catalog_indices:
            if ((ref_catalog_ras[ref_index]-ra)**2.
                + (ref_catalog_decs[ref_index]-dec)**2.
                < match_radius**2.
            ):
                match_index_list[i] = ref_index
                ref_catalog_indices.remove(ref_index)
                break

    match_index_list = np.asarray(match_index_list)
    unmatched_indices = np.where(match_index_list == -1.)[0]
    decon_catalog_match_ratio = 1.-float(len(unmatched_indices))/float(len(match_index_list))
    ref_catalog_match_ratio = float(
        len(ref_catalog_limited)-len(match_index_list)
        )/float(len(ref_catalog_limited))

    matched_decon_catalog_fluxes = [source['flux']['I'] for i, source in enumerate(decon_catalog_limited) if match_index_list[i] != -1.]
    matched_ref_catalog_fluxes = [ref_catalog_limited[int(ind)]['flux']['I'] for ind in match_index_list if int(ind) != -1]

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

    plot_flux_scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes,
                      fit_param, goodness_of_fit, decon_catalog_match_ratio,
                      '{}/{}_flux_compare.png'.format(output_path, obsid))

    # Calculate and plot positional offsets
    matched_decon_catalog_ras = [source['ra'] for i, source in enumerate(decon_catalog_limited) if match_index_list[i] != -1.]
    matched_decon_catalog_decs = [source['dec'] for i, source in enumerate(decon_catalog_limited) if match_index_list[i] != -1.]
    matched_ref_catalog_ras = [ref_catalog_limited[int(ind)]['ra'] for ind in match_index_list if int(ind) != -1]
    matched_ref_catalog_decs = [ref_catalog_limited[int(ind)]['dec'] for ind in match_index_list if int(ind) != -1]

    ra_offsets = [(
        matched_decon_catalog_ras[i]-matched_ref_catalog_ras[i]
    )*60. for i in range(len(matched_ref_catalog_ras))]
    dec_offsets = [(
        matched_decon_catalog_decs[i]-matched_ref_catalog_decs[i]
    )*60. for i in range(len(matched_ref_catalog_decs))]

    plot_pos_offsets_scatter(ra_offsets, dec_offsets,
        '{}/{}_pos_offsets_scatter.png'.format(output_path, obsid)
    )
    plot_pos_offsets_vectors(matched_ref_catalog_ras, matched_ref_catalog_decs,
        ra_offsets, dec_offsets, match_radius, search_radius,
        '{}/{}_pos_offsets_vectors.png'.format(output_path, obsid)
    )


def plot_flux_scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes,
                      fit_param, goodness_of_fit, decon_catalog_match_ratio,
                      saveloc):
    plt.figure()
    # plot the 1-to-1 line
    plt.plot([min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
    [min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
        linestyle=':', linewidth=0.5, color='blue', label='1-to-1 line')
    # plot the power law fit
    plt.plot([min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
    [fit_param*min([min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)]), fit_param*max([max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)])],
        linestyle=':', linewidth=0.5, color='red', label='power law fit: A={:.3f}, X^2={:.2f}'.format(float(fit_param), float(goodness_of_fit)))
    plt.scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes,
                marker='o', s=1., color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('GLEAM flux density (Jy)')
    plt.ylabel('Deconvolved flux density (Jy)')
    plt.title('Deconvolved Flux Density Agreement With GLEAM: '
              'Match Ratio {:.3f}'.format(decon_catalog_match_ratio))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc='upper left')
    plt.savefig(saveloc)
    plt.close()


def plot_flux_hist(ref_catalog_fluxes, decon_catalog_fluxes, saveloc):

    # get bin edges
    null, bin_edges = np.histogram(
        [np.log(val) for val in ref_catalog_fluxes.extend(decon_cal_fluxes)],
        bins='auto'
    )
    # histogram data
    ref_hist, bin_edges = np.histogram(
        [np.log(val) for val in ref_catalog_fluxes], bins=bin_edges
    )
    decon_hist, bin_edges = np.histogram(
        [np.log(val) for val in decon_catalog_fluxes], bins=bin_edges
    )
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
    plt.plot(xvals, ref_yvals, linestyle='-', color='blue', linewidth=0.5,
             label='GLEAM flux densities')
    plt.plot(xvals, decon_yvals, linestyle='-', color='red', linewidth=0.5,
             label='Deconvolved flux densities')
    plt.fill_between(xvals, 0, ref_yvals, facecolor='blue', alpha=0.1)
    plt.fill_between(xvals, 0, decon_yvals, facecolor='red', alpha=0.1)
    #plt.yscale('log')
    plt.xlabel('Flux Densities (Jy)')
    plt.ylabel('Histogram Count')
    plt.title('Source Flux Density Histogram')
    plt.legend(loc='upper right')
    plt.savefig(saveloc)
    plt.close()


def plot_pos_offsets_scatter(ra_offsets, dec_offsets, saveloc):
    plt.figure()
    plt.grid(True, zorder=0)
    plt.scatter(ra_offsets, dec_offsets, marker='o', s=1., color='black', zorder=10)
    plt.xlabel('RA offset (arcmin)')
    plt.ylabel('Dec offset (arcmin)')
    plt.title('Source Positional Offsets, Deconvolved - GLEAM')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([-match_radius*60., match_radius*60., -match_radius*60., match_radius*60.])
    plt.savefig(saveloc)
    plt.close()


def plot_pos_offsets_vectors(matched_ref_catalog_ras, matched_ref_catalog_decs,
                             ra_offsets, dec_offsets,
                             match_radius, search_radius,
                             saveloc):
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
    plt.savefig(saveloc)
    plt.close()


if __name__=='__main__':
    match_catalogs()
