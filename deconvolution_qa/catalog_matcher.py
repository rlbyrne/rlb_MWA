#!/usr/bin/python

import scipy.io
import matplotlib.pyplot as plt
import numpy as np


def match_catalogs():

    fhd_run_path = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018'
    obsid = '1130773024'
    ref_catalog_path = '/Users/rubybyrne/FHD/catalog_data/GLEAM_plus_rlb2017.sav'
    search_radius = 12.
    match_radius = .05  # a match is within 3 arcmin of the source

    obs_sav = scipy.io.readsav(
        '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
        )
    obs_ra = float(obs_sav['obs']['obsra'][0])
    obs_dec = float(obs_sav['obs']['obsdec'][0])

    decon_catalog = scipy.io.readsav(
        '{}/deconvolution/{}_fhd.sav'.format(fhd_run_path, obsid)
        )['source_array']
    decon_catalog_limited = []
    for source in decon_catalog:
        if min((float(source['ra'])-obs_ra)**2.,
               (float(source['ra'])-obs_ra+360.)**2.,
               (float(source['ra'])-obs_ra-360.)**2.
               ) + (float(source['dec'])-obs_dec)**2. < search_radius**2.:
            decon_catalog_limited.append(source)
    decon_catalog_limited.sort(key=lambda x: float(x['flux']['I']), reverse=True)

    decon_catalog_fluxes = [source['flux']['I'] for source in decon_catalog_limited]

    ref_catalog = scipy.io.readsav(ref_catalog_path)['catalog']
    ref_catalog_limited = []
    for source in ref_catalog:
        if min((float(source['ra'])-obs_ra)**2.,
               (float(source['ra'])-obs_ra+360.)**2.,
               (float(source['ra'])-obs_ra-360.)**2.
               ) + (float(source['dec'])-obs_dec)**2. < search_radius**2.:
            ref_catalog_limited.append(source)
    ref_catalog_limited.sort(key=lambda x: float(x['flux']['I']), reverse=True)
    ref_catalog_ras = [float(source['ra']) for source in ref_catalog_limited]
    ref_catalog_decs = [float(source['dec']) for source in ref_catalog_limited]

    ref_catalog_fluxes = [source['flux']['I'] for source in ref_catalog_limited]

    match_index_list = [-1.]*len(decon_catalog_limited)
    for i, decon_source in enumerate(decon_catalog_limited):
        ra = float(decon_source['ra'])
        dec = float(decon_source['dec'])
        for ref_index in range(len(ref_catalog_limited)):
            if (ref_catalog_ras[ref_index]-ra)**2. + (ref_catalog_decs[ref_index]-dec)**2. < match_radius**2.:
                match_index_list[i] = ref_index
                del ref_catalog_limited[ref_index]
                break

    match_index_list = np.array(match_index_list)
    unmatched_indices = np.where(match_index_list == -1.)[0]
    decon_catalog_match_ratio = 1.-float(len(unmatched_indices))/float(len(match_index_list))
    ref_catalog_match_ratio = float(
        len(ref_catalog_limited)-len(match_index_list)
        )/float(len(ref_catalog_limited))
    print len(match_index_list)
    print len(ref_catalog_limited)
    print decon_catalog_match_ratio
    print ref_catalog_match_ratio
    print match_index_list[0:200]


if __name__=='__main__':
    match_catalogs()
