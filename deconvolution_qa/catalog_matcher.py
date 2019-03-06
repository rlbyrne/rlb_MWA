#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os


def match_catalogs_wrapper(
    fhd_run_path, output_path,
    ref_catalog_path='/Users/rubybyrne/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav'
):

    # Find all obs files in the run path directory
    obs_file_list = os.listdir('{}/metadata/'.format(fhd_run_path))
    obs_file_list = [
        filename.replace(
            '{}/metadata/'.format(fhd_run_path), ''
        )[0:10] for filename in obs_file_list if filename.endswith('_obs.sav')
    ]
    obs_file_list = list(set(obs_file_list))

    # Find all deconvolution output files in the run path directory
    decon_catalog_file_list = os.listdir(
        '{}/deconvolution/'.format(fhd_run_path)
    )
    decon_catalog_file_list = [
        filename.replace(
            '{}/deconvolution/'.format(fhd_run_path), ''
        )[0:10] for filename in decon_catalog_file_list
        if filename.endswith('_fhd.sav')
    ]
    decon_catalog_file_list = list(set(decon_catalog_file_list))

    # Run all obsids that have both their obs files and deconvolution outputs
    obsid_list = [
        obs for obs in obs_file_list if obs in decon_catalog_file_list
    ]

    print 'Restoring GLEAM catalog'
    ref_catalog = scipy.io.readsav(ref_catalog_path)['catalog']

    for i, obsid in enumerate(obsid_list):
        print 'Creating deconvolution QA plots for {}, obsid {}/{}'.format(
            obsid, i+1, len(obsid_list)
        )
        match_catalogs(fhd_run_path, output_path, ref_catalog, obsid)


def match_catalogs(fhd_run_path, output_path, ref_catalog, obsid):

    search_radius = 12.  # match sources out to 12 degrees from pointing center
    match_radius = .04  # a match is within this many degrees of the source
    flux_percent_tolerance = .25  # a match is within this flux tolerance

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
    decon_catalog_limited.sort(
        key=lambda x: float(x['flux']['I']), reverse=True
    )

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
    #ref_catalog_limited = ref_catalog_limited[0:100]  # for debugging

    plot_flux_hist([source['flux']['I'] for source in ref_catalog_limited],
                   [source['flux']['I'] for source in decon_catalog_limited],
                   '{}/{}_flux_hist.png'.format(output_path, obsid))

    ref_catalog_ras = [float(source['ra']) for source in ref_catalog_limited]
    ref_catalog_decs = [float(source['dec']) for source in ref_catalog_limited]
    ref_catalog_fluxes = [float(source['flux']['I'])
                          for source in ref_catalog_limited]

    # match sources
    ref_match_index_list = []
    decon_match_index_list = []
    ref_catalog_indices = list(range(len(ref_catalog_limited)))
    for i, decon_source in enumerate(decon_catalog_limited):
        ra = float(decon_source['ra'])
        dec = float(decon_source['dec'])
        flux = float(decon_source['flux']['I'])
        for ref_index in ref_catalog_indices:
            if ((ref_catalog_ras[ref_index]-ra)**2.
                + (ref_catalog_decs[ref_index]-dec)**2.
                < match_radius**2.
                and abs(ref_catalog_fluxes[ref_index]-flux) < flux_percent_tolerance*flux
            ):
                ref_match_index_list.append(ref_index)
                decon_match_index_list.append(i)
                ref_catalog_indices.remove(ref_index)
                break
    ref_match_index_list = np.array(ref_match_index_list)
    decon_match_index_list = np.array(decon_match_index_list)

    n_unmatched = len(decon_catalog_limited) - len(ref_match_index_list)
    decon_catalog_match_ratio = 1.-float(n_unmatched)/float(len(ref_match_index_list))
    ref_catalog_match_ratio = float(
        len(ref_catalog_limited)-len(ref_match_index_list)
        )/float(len(ref_catalog_limited))

    matched_decon_catalog_fluxes = np.array(decon_catalog_limited[
        decon_match_index_list
    ]['flux']['I'])
    matched_ref_catalog_fluxes = np.array(ref_catalog_fluxes[ref_match_index_list])

    fit_param = np.sum(
        (matched_decon_catalog_fluxes/matched_ref_catalog_fluxes)**2
    )/np.sum(matched_decon_catalog_fluxes/matched_ref_catalog_fluxes)
    goodness_of_fit = np.sum((
            matched_ref_catalog_fluxes-(matched_decon_catalog_fluxes/fit_param)
    )**2/matched_ref_catalog_fluxes**2)/len(matched_ref_catalog_fluxes)

    plot_flux_scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes,
                      fit_param, goodness_of_fit, decon_catalog_match_ratio,
                      '{}/{}_flux_compare.png'.format(output_path, obsid))

    # Calculate and plot positional offsets
    matched_decon_catalog_ras = [
        source['ra'] for i, source in enumerate(decon_catalog_limited)
        if match_index_list[i] != -1.
    ]
    matched_decon_catalog_decs = [
        source['dec'] for i, source in enumerate(decon_catalog_limited)
        if match_index_list[i] != -1.
    ]
    matched_ref_catalog_ras = [
        ref_catalog_limited[int(ind)]['ra'] for ind in match_index_list
        if int(ind) != -1
    ]
    matched_ref_catalog_decs = [
        ref_catalog_limited[int(ind)]['dec'] for ind in match_index_list
        if int(ind) != -1
    ]

    ra_offsets = [(
        matched_decon_catalog_ras[i]-matched_ref_catalog_ras[i]
    )*60. for i in range(len(matched_ref_catalog_ras))]
    dec_offsets = [(
        matched_decon_catalog_decs[i]-matched_ref_catalog_decs[i]
    )*60. for i in range(len(matched_ref_catalog_decs))]

    plot_pos_offsets_scatter(ra_offsets, dec_offsets, match_radius,
        '{}/{}_pos_offsets_scatter.png'.format(output_path, obsid)
    )
    plot_pos_offsets_vectors(matched_ref_catalog_ras, matched_ref_catalog_decs,
        ra_offsets, dec_offsets, obs_ra, obs_dec, match_radius, search_radius,
        '{}/{}_pos_offsets_vectors.png'.format(output_path, obsid)
    )


def plot_flux_scatter(
    matched_ref_catalog_fluxes, matched_decon_catalog_fluxes, fit_param,
    goodness_of_fit, decon_catalog_match_ratio, saveloc
):

    plt.figure()
    # plot the 1-to-1 line
    plt.plot(
        [np.min([
            np.min(matched_ref_catalog_fluxes),
            np.min(matched_decon_catalog_fluxes)
        ]), np.max([
            np.max(matched_ref_catalog_fluxes),
            np.max(matched_decon_catalog_fluxes)
        ])],
        [np.min([
            np.min(matched_ref_catalog_fluxes),
            np.min(matched_decon_catalog_fluxes)
        ]), np.max([
            np.max(matched_ref_catalog_fluxes),
            np.max(matched_decon_catalog_fluxes)
        ])],
        linestyle=':', linewidth=0.5, color='blue', label='1-to-1 line'
    )
    # plot the power law fit
    plt.plot(
        [np.min([
            np.min(matched_ref_catalog_fluxes),
            np.min(matched_decon_catalog_fluxes)
        ]), np.max([
            np.max(matched_ref_catalog_fluxes),
            np.max(matched_decon_catalog_fluxes)
        ])],
        [
            fit_param*np.min([np.min(matched_ref_catalog_fluxes),
            np.min(matched_decon_catalog_fluxes)]),
            fit_param*np.max([np.max(matched_ref_catalog_fluxes),
            np.max(matched_decon_catalog_fluxes)])
        ],
        linestyle=':',
        linewidth=0.5,
        color='red',
        label='power law fit: A={:.3f}, X^2={:.2f}'.format(
            float(fit_param), float(goodness_of_fit)
        )
    )
    print np.shape(matched_ref_catalog_fluxes)
    print np.shape(matched_decon_catalog_fluxes)
    plt.scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes,
                marker='o', s=1., color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('GLEAM flux density (Jy)')
    plt.ylabel('Deconvolved flux density (Jy)')
    plt.title('Deconvolved Flux Density Agreement With GLEAM: '
              'Match Ratio {:.3f}'.format(decon_catalog_match_ratio))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([
        min([
            min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)
        ]),
        max([
            max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)
        ]),
        min([
            min(matched_ref_catalog_fluxes), min(matched_decon_catalog_fluxes)
        ]),
        max([
            max(matched_ref_catalog_fluxes), max(matched_decon_catalog_fluxes)
        ])
    ])
    plt.legend(loc='upper left')
    print 'Saving plot to {}'.format(saveloc)
    plt.savefig(saveloc, dpi=300)
    plt.close()


def plot_flux_hist(ref_catalog_fluxes, decon_catalog_fluxes, saveloc):

    # get rid of zeros and negatives to allow for log scaling
    ref_catalog_fluxes = np.array(ref_catalog_fluxes)
    ref_catalog_fluxes = ref_catalog_fluxes[np.where(ref_catalog_fluxes>0)[0]]
    decon_catalog_fluxes = np.array(decon_catalog_fluxes)
    decon_catalog_fluxes = decon_catalog_fluxes[np.where(decon_catalog_fluxes>0)[0]]

    # get bin edges
    null, bin_edges = np.histogram(
        np.log(np.append(ref_catalog_fluxes, decon_catalog_fluxes)),
        bins='auto'
    )
    # histogram data
    ref_hist, bin_edges = np.histogram(
        np.log(ref_catalog_fluxes), bins=bin_edges
    )
    decon_hist, bin_edges = np.histogram(
        np.log(decon_catalog_fluxes), bins=bin_edges
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
    print 'Saving plot to {}'.format(saveloc)
    plt.savefig(saveloc, dpi=300)
    plt.close()


def plot_pos_offsets_scatter(ra_offsets, dec_offsets, match_radius, saveloc):

    plt.figure()
    plt.grid(True, zorder=0)
    plt.scatter(
        ra_offsets, dec_offsets, marker='o', s=1., color='black', zorder=10
    )
    plt.xlabel('RA offset (arcmin)')
    plt.ylabel('Dec offset (arcmin)')
    plt.title('Source Positional Offsets, Deconvolved - GLEAM')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([
        -match_radius*60., match_radius*60., -match_radius*60., match_radius*60.
    ])
    print 'Saving plot to {}'.format(saveloc)
    plt.savefig(saveloc, dpi=300)
    plt.close()


def plot_pos_offsets_vectors(
    matched_ref_catalog_ras, matched_ref_catalog_decs, ra_offsets, dec_offsets,
    obs_ra, obs_dec, match_radius, search_radius, saveloc
):

    arrow_scale_factor = .5
    plt.figure()
    for i in range(len(ra_offsets)):
        if matched_ref_catalog_ras[i]-obs_ra > search_radius:
            matched_ref_catalog_ras[i] = matched_ref_catalog_ras[i]-360.
        elif matched_ref_catalog_ras[i]-obs_ra < -search_radius:
            matched_ref_catalog_ras[i] = matched_ref_catalog_ras[i]+360.
    plt.quiver(matched_ref_catalog_ras, matched_ref_catalog_decs, ra_offsets,
               dec_offsets)
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.title('Source Positional Offsets')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis([obs_ra-search_radius-2., obs_ra+search_radius+2.,
              obs_dec-search_radius-2., obs_dec+search_radius+2.])
    print 'Saving plot to {}'.format(saveloc)
    plt.savefig(saveloc, dpi=300)
    plt.close()


if __name__=='__main__':
    match_catalogs_wrapper(
        '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Feb2019',
        '/Users/rubybyrne/diffuse_survey_qa/4pol_decon_Feb2019'
    )
    #write_params_to_csv(
    #    '/Users/rubybyrne/diffuse_survey_qa/test_params.csv',
    #    'testparam', [1131556544, 1130783824], [100, 1], overwrite=True)
