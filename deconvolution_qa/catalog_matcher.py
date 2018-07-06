#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np


def match_catalogs():

    fhd_run_path = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018'
    obsid = '1130773024'
    ref_catalog_path = '/Users/rubybyrne/FHD/catalog_data/GLEAM_plus_rlb2017.sav'
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

    print match_index_list
    match_index_list = np.asarray(match_index_list)
    unmatched_indices = np.where(match_index_list == -1.)[0]
    decon_catalog_match_ratio = 1.-float(len(unmatched_indices))/float(len(match_index_list))
    ref_catalog_match_ratio = float(
        len(ref_catalog_limited)-len(match_index_list)
        )/float(len(ref_catalog_limited))
    print len(match_index_list)
    print len(ref_catalog_limited)
    print decon_catalog_match_ratio
    print ref_catalog_match_ratio

    matched_decon_catalog_fluxes = [source['flux']['I'] for i, source in enumerate(decon_catalog_limited) if match_index_list[i] != -1.]
    matched_ref_catalog_fluxes = [ref_catalog_limited[int(ind)]['flux']['I'] for ind in match_index_list if int(ind) != -1]
    print len(matched_decon_catalog_fluxes)
    print len(matched_ref_catalog_fluxes)

    print 'DECON CATALOG FLUXES:    '
    print matched_decon_catalog_fluxes
    print 'REFERENCE CATALOG FLUXES:    '
    print matched_ref_catalog_fluxes

"""    plt.figure()
    plt.scatter(matched_ref_catalog_fluxes, matched_decon_catalog_fluxes, marker='o', s=1.)
    # plt.xlim(ra_plot_range[0],ra_plot_range[1])
    # plt.ylim(dec_plot_range[0],dec_plot_range[1])
    plt.xlabel('East/West Location (m)')
    plt.ylabel('North/South Location (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(
        '/Users/rubybyrne/array_simulation/'
        'antenna_hex_{}m.png'.format(int(antenna_spacing))
    )
    plt.close()"""



if __name__=='__main__':
    match_catalogs()
