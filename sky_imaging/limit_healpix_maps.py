#!/usr/bin/python

from astropy.io import fits
import scipy.io
import scipy
import numpy as np
import healpy as hp
import sys
import os
import healpix_utils
import plot_healpix_map

obs_list_1 = [
    '1131551744',
    '1130783824',
    '1131562544',
    '1131709912',
    '1130776864',
    '1131461496',
    '1130782264',
    #'1131454176', high power and systematics in Stokes V
    '1131715432',
    '1131733552',
    '1131542624',
    '1130773144',
    '1131461376',
    '1131557144',
    '1131454296',
    '1131731752',
    '1130778664',
    '1131470496',
    '1131559064',
    '1131717232',
    '1131463536',
    '1130773264',
    '1131463416',
    '1131717352',
    '1131713632',
    '1131478056',
    '1131468936',
    '1131468696',
    '1131535424',
    '1131463296',
    '1131465216',
    '1131710032',
    '1130776624',
    '1131456096',
    #'1131456216',
    '1131540824',
    '1131711952',
    '1131459576',
    '1131477936',
    '1131733672',
    '1131564464',
    '1130787784',
    #'1131475896',
    '1131461616',
    '1131558944',
    '1131470616',
    '1131549944',
    '1131553544',
    #'1131477816',
    '1131459696',
    '1130780464',
    '1131726352',
    #'1131715312',
    '1131470736',
    '1131548024',
    '1131710152',
    '1130785864',
    #'1131724672',
    '1131544424'
]
obs_list_2 = ['1131542504',
    #'1131717112',
    '1131733432',
    '1131735232',
    '1131553664',
    '1131724432',
    '1131542744',
    '1131455976',
    '1131719152',
    '1131454416',
    #'1131728032',
    '1130787544',
    '1130776744',
    #'1131726472',
    '1130780224',
    '1131551624',
    '1131722632',
    '1131547904',
    '1130776624',
    '1131562664',
    '1131550064',
    '1131537104',
    '1131555224',
    '1131467136',
    '1131539024',
    '1131555344',
    '1131546104',
    '1131548144',
    '1131472416',
    '1131558824',
    '1131544304',
    '1130789584',
    '1131476136',
    '1130789344',
    #'1131728272',
    '1131722872',
    '1130785744',
    '1131730072',
    '1131459816',
    '1131564584',
    '1131457776',
    '1131724552',
    '1130787664',
    '1130778424',
    '1131728152',
    '1131722752',
    '1131538904',
    '1131544544',
    '1130778544',
    '1131467016',
    '1131546344',
    '1130789464',
    '1131713512',
    '1131546224',
    '1131474336',
    '1130782144',
    '1131735472',
    '1130775064',
    '1130774824',
    '1131720832',
    '1130774944',
    '1131557264',
    '1130783944',
    #'1131713752',
    '1131472296',
    '1131465096',
    '1131457896',
    '1131555464',
    #'1131720712',
    #'1131711832',
    '1131562424',
    '1131551864',
    '1131540704',
    '1130780344',
    '1131731632',
    '1131468816',
    #'1131711712',
    '1131472536',
    #'1131729832',
    '1130773024',
    #'1131720952',
    #'1131718912',
    #'1131719032',
    '1131474096',
    '1131465336',
    '1131715552',
    '1131458016',
    '1131540944',
    '1131557024',
    '1131731872',
    '1131553424',
    '1131560864',
    '1130784064',
    '1131466896',
    '1130782024',
    '1131560624',
    '1131474216',
    '1131564344',
    '1131729952',
    '1131560744',
    '1130785624',
    '1131709432',
    '1131536624',
    '1131536384',
    '1131711112',
    '1131709192',
    '1131710992',
    #'1131709792', stripes in V
    '1131453456',
    '1131565304',
    '1131478776',
    '1131566504',
    '1131565184',
    '1131566624',
    '1131566744',
    '1131565064',
    '1131567944',
    '1131478656',
    '1131568544',
    #'1131740872', excess power
    '1131739432',
    '1130788504',
    '1130788264',
    '1131740752',
    #'1131735952', # maybe excess power
    #'1131739552', excess power
    '1131455736',
    '1131710392',
    '1131708952',
    '1131457176',
    '1131716512',
    #'1131713272', excess power
    '1131458976',
    '1131712192',
    '1131453936',
    '1131457536',
    '1131537704',
    '1131543584'
]

def find_obs_groups():

    field_center = hp.pixelfunc.ang2vec(0., -27., lonlat=True)
    field_radius = 13.5

    map = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_May2020/nsamples_map_more_obs.fits'
    )

    pix_use = []
    signal_use = []
    for ind, pix in enumerate(map.pix_arr):
        pix_vec = hp.pix2vec(map.nside, pix, nest=map.nest)
        dist_deg = hp.rotator.angdist(pix_vec, field_center)*180./np.pi
        if dist_deg <= field_radius:
            pix_use.append(pix)
            signal_use.append(map.signal_arr[ind])

    map.signal_arr = np.array(signal_use)
    map.pix_arr = np.array(pix_use)

    map.write_data_to_fits('/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/nsamples_eor0.fits')
    plot_healpix_map.plot_filled_pixels(
        map,
        '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/nsamples_eor0.png'
    )

    obs_use_list = []
    pix_use_list = []
    for obsid in obs_list_1+obs_list_2:
        if obsid in obs_list_1:
            path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020'
        elif obsid in obs_list_2:
            path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'

        obs_struct = scipy.io.readsav(
            '{}/metadata/{}_obs.sav'.format(path, obsid)
        )['obs']
        obs_vec = hp.pixelfunc.ang2vec(
            float(obs_struct['obsra']), float(obs_struct['obsdec']),
            lonlat=True
        )
        dist_deg = hp.rotator.angdist(obs_vec, field_center)*180./np.pi
        if dist_deg <= 30.:
            obs_reference_map = healpix_utils.load_map(
                '{}/output_data/{}_weighted_Residual_I_HEALPix.fits'.format(path, obsid)
            )
            if np.shape(np.intersect1d(
                obs_reference_map.pix_arr, map.pix_arr
            ))[0] > 0:
                obs_use_list.append(obsid)
                pix_use_list.append(obs_reference_map.pix_arr)

    least_sampled_pixels = map.pix_arr[
        np.where(map.signal_arr == np.min(map.signal_arr))[0]
    ]
    obs_use_least_sampled = []
    obs_to_match = [obsid for obsid in obs_use_list]
    for obsind, obsid in enumerate(obs_use_list):
        if np.shape(np.intersect1d(least_sampled_pixels, pix_use_list[obsind]))[0] > 0:
            obs_use_least_sampled.append(obsid)
            obs_to_match.remove(obsid)

    print np.shape(obs_use_least_sampled)

    obs_groups = []

    # Check for one-obs groups
    obs_to_remove = []
    for obsid in obs_use_least_sampled:
        if np.shape(np.intersect1d(
            pix_use_list[obs_use_list.index(obsid)], map.pix_arr
        ))[0] == np.shape(map.pix_arr)[0]:
            print 'Obsid {} is an obs group'.format(obsid)
            obs_groups.append([obsid])
            obs_to_remove.append(obsid)
    for obsid in obs_to_remove:
        obs_use_least_sampled.remove(obsid)

    # Check for two-obs groups
    if len(obs_use_least_sampled) > 0:
        obs_to_remove = []
        for obsid in obs_use_least_sampled:
            #see if two obs cover the whole field
            pix_list_1 = pix_use_list[obs_use_list.index(obsid)]
            no_match = True
            obsind = 0
            while no_match and obsind < len(obs_to_match):
                obsid_2 = obs_to_match[obsind]
                print obsid_2
                print obs_use_list.index(obsid_2)
                pix_list_2 = pix_use_list[obs_use_list.index(obsid_2)]
                if np.shape(np.intersect1d(
                    np.concatenate((pix_list_1, pix_list_2)), map.pix_arr
                ))[0] == np.shape(map.pix_arr)[0]:
                    print 'Obsids {} and {} are an obs group'.format(obsid, obsid_2)
                    obs_groups.append([obsid, obsid_2])
                    obs_to_match.remove(obsid_2)
                    obs_to_remove.append(obsid)
                    no_match = False
                obsind += 1
        for obsid in obs_to_remove:
            obs_use_least_sampled.remove(obsid)

    # Check for three-obs groups
    if len(obs_use_least_sampled) > 0:
        obs_to_remove = []
        for obsid in obs_use_least_sampled:
            #see if two obs cover the whole field
            pix_list_1 = pix_use_list[obs_use_list.index(obsid)]
            no_match = True
            obsind_2 = 0
            while no_match and obsind_2 < len(obs_to_match):
                obsid_2 = obs_to_match[obsind_2]
                pix_list_2 = pix_use_list[obs_use_list.index(obsid_2)]
                obsind_3 = 0
                while no_match and obsind_3 < len(obs_to_match):
                    if obsind_2 != obsind_3:
                        obsid_3 = obs_to_match[obsind_3]
                        pix_list_3 = pix_use_list[obs_use_list.index(obsid_3)]
                        if np.shape(np.intersect1d(
                            np.concatenate((pix_list_1, pix_list_2, pix_list_3)),
                            map.pix_arr
                        ))[0] == np.shape(map.pix_arr)[0]:
                            print 'Obsids {}, {}, and {} are an obs group'.format(
                                obsid, obsid_2, obsid_3
                            )
                            obs_groups.append([obsid, obsid_2, obsid_3])
                            obs_to_match.remove(obsid_2)
                            obs_to_match.remove(obsid_3)
                            obs_to_remove.append(obsid)
                            no_match = False
                    obsind_3 += 1
                obsind_2 += 1
        for obsid in obs_to_remove:
            obs_use_least_sampled.remove(obsid)

    print obs_groups


def make_combined_images():

    field_center = hp.pixelfunc.ang2vec(0., -27., lonlat=True)
    field_radius = 13.5

    map = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_May2020/nsamples_map_more_obs.fits'
    )

    pix_use = []
    signal_use = []
    for ind, pix in enumerate(map.pix_arr):
        pix_vec = hp.pix2vec(map.nside, pix, nest=map.nest)
        dist_deg = hp.rotator.angdist(pix_vec, field_center)*180./np.pi
        if dist_deg <= field_radius:
            pix_use.append(pix)
            signal_use.append(map.signal_arr[ind])

    map.signal_arr = np.array(signal_use)
    map.pix_arr = np.array(pix_use)

    obs_groups = [
        ['1131454296'],
        ['1131455736'],
        ['1131710392'],
        ['1131708952'],
        ['1131710032', '1131713632'],
        ['1131457176', '1131715432'],
        ['1131453936', '1131716512'],
        ['1131537104', '1131709912', '1131456096']
    ]

    if False:
        for pol in ['I', 'Q', 'U', 'V']:
            full_map = healpix_utils.load_map(
                '/Users/rubybyrne/diffuse_survey_plotting_May2020/Stokes{}_average_map_more_obs.fits'.format(pol)
            )
            full_map.signal_arr = np.array([
                full_map.signal_arr[
                    np.where(full_map.pix_arr == pix)[0]
                ][0] for pix in map.pix_arr
            ])
            full_map.pix_arr = map.pix_arr
            full_map.write_data_to_fits(
                '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/Stokes{}_eor0_full_map.fits'.format(pol)
            )
            if pol == 'I':
                colorbar_range = [-1e4, 1e4]
            else:
                colorbar_range = [-2e3, 2e3]
            plot_healpix_map.plot_filled_pixels(
                full_map,
                '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/Stokes{}_eor0_full_map.png'.format(pol),
                colorbar_range=colorbar_range
            )

    for group_ind, group in enumerate(obs_groups):

        list_1 = [obs for obs in group if obs in obs_list_1]
        list_2 = [obs for obs in group if obs in obs_list_2]
        obs_lists_together = [list_1, list_2]
        paths_list = [
            '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
            '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'
        ]
        if len(list_1) == 0:
            paths_list = paths_list[1]
            obs_lists_together = obs_lists_together[1]
        if len(list_2) == 0:
            paths_list = paths_list[0]
            obs_lists_together = obs_lists_together[0]
        combined_maps, weight_maps = healpix_utils.average_healpix_maps(
            paths_list,
            obs_lists = obs_lists_together,
            nside=128,
            cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
            weighting='weighted',
            apply_radial_weighting=True,
            apply_rm_correction=True,
            rm_file='/Users/rubybyrne/diffuse_survey_rm_empirical.csv'
        )
        pols = ['I', 'Q', 'U', 'V']
        for map_ind in range(len(combined_maps)):
            pol = pols[map_ind]
            combined_maps[map_ind].signal_arr = np.array([
                combined_maps[map_ind].signal_arr[
                    np.where(combined_maps[map_ind].pix_arr == pix)[0]
                ][0] for pix in map.pix_arr
            ])
            combined_maps[map_ind].pix_arr = map.pix_arr

            combined_maps[map_ind].write_data_to_fits(
                '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/Stokes{}_eor0_map{}_tapered_empirical_rm.fits'.format(pol, group_ind+1)
            )
            if pol == 'I':
                colorbar_range = [-1e4, 1e4]
            else:
                colorbar_range = [-5e3, 5e3]
            plot_healpix_map.plot_filled_pixels(
                combined_maps[map_ind],
                '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/Stokes{}_eor0_map{}_tapered_empirical_rm.png'.format(pol, group_ind+1),
                colorbar_range=colorbar_range
            )


def plot_individual_images():

    field_center = hp.pixelfunc.ang2vec(0., -27., lonlat=True)
    field_radius = 13.5

    map = healpix_utils.load_map(
        '/Users/rubybyrne/diffuse_survey_plotting_May2020/nsamples_map_more_obs.fits'
    )

    pix_use = []
    signal_use = []
    for ind, pix in enumerate(map.pix_arr):
        pix_vec = hp.pix2vec(map.nside, pix, nest=map.nest)
        dist_deg = hp.rotator.angdist(pix_vec, field_center)*180./np.pi
        if dist_deg <= field_radius:
            pix_use.append(pix)
            signal_use.append(map.signal_arr[ind])

    map.signal_arr = np.array(signal_use)
    map.pix_arr = np.array(pix_use)

    obs_groups = [
        ['1131454296'],
        ['1131455736'],
        ['1131710392'],
        ['1131708952'],
        ['1131710032', '1131713632'],
        ['1131457176', '1131715432'],
        ['1131453936', '1131716512'],
        ['1131537104', '1131709912', '1131456096']
    ]


    for group_ind, group in enumerate(obs_groups):

        if len(group) > 1:
            for obs_ind, obs in enumerate(group):
                pols = ['I', 'Q', 'U', 'V']
                if obs in obs_list_1:
                    path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
                else:
                    path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'

                combined_maps, weight_maps = healpix_utils.average_healpix_maps(
                    path,
                    obs_lists = [obs],
                    nside=128,
                    cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
                    weighting='weighted',
                    apply_radial_weighting=False,
                    apply_rm_correction=True
                )
                pols = ['I', 'Q', 'U', 'V']
                for map_ind in range(len(combined_maps)):
                    pol = pols[map_ind]
                    add_pixels = np.setdiff1d(
                        map.pix_arr, combined_maps[map_ind].pix_arr
                    )
                    combined_maps[map_ind].pix_arr = np.concatenate((
                        np.array(combined_maps[map_ind].pix_arr), add_pixels
                    ))
                    combined_maps[map_ind].signal_arr = np.concatenate((
                        np.array(combined_maps[map_ind].signal_arr),
                        np.zeros(np.shape(add_pixels)[0])
                    ))
                    combined_maps[map_ind].signal_arr = np.array([
                        combined_maps[map_ind].signal_arr[
                            np.where(combined_maps[map_ind].pix_arr == pix)[0]
                        ][0] for pix in map.pix_arr
                    ])
                    combined_maps[map_ind].pix_arr = map.pix_arr

                    if pol == 'I':
                        colorbar_range = [-1e4, 1e4]
                    else:
                        colorbar_range = [-5e3, 5e3]
                    plot_healpix_map.plot_filled_pixels(
                        combined_maps[map_ind],
                        '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots/Stokes{}_eor0_map{}_obs{}.png'.format(pol, group_ind+1, obs_ind+1),
                        colorbar_range=colorbar_range
                    )


if __name__ == '__main__':

    make_combined_images()
