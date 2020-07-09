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
import csv

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
c = 3.e8


def get_effective_rotation_angles(rms, start_freq_mhz, end_freq_mhz):

    wl_max = c/(start_freq_mhz*1.e6)
    wl_min = c/(end_freq_mhz*1.e6)
    fresS_min, fresC_min = scipy.special.fresnel(2*np.sqrt(rms/np.pi+0j)*wl_min)
    fresS_max, fresC_max = scipy.special.fresnel(2*np.sqrt(rms/np.pi+0j)*wl_max)
    cos_int = (
        np.cos(2.*rms*wl_min**2.)/wl_min
        - np.cos(2.*rms*wl_max**2.)/wl_max
        + 2*np.sqrt(np.pi*rms+0j)*(fresS_min-fresS_max)
    )
    sin_int = (
        np.sin(2.*rms*wl_min**2.)/wl_min
        - np.sin(2.*rms*wl_max**2.)/wl_max
        - 2*np.sqrt(np.pi*rms+0j)*(fresC_min-fresC_max)
    )
    rot_angles = np.arctan2(np.real(sin_int), np.real(cos_int))
    return rot_angles


def create_rm_lookup_table(start_freq_mhz, end_freq_mhz):

    min_rm = -4
    max_rm = 1
    stepsize = 1e-7
    rms = np.arange(min_rm, max_rm, stepsize)
    rot_angles = get_effective_rotation_angles(rms, start_freq_mhz, end_freq_mhz)
    return rot_angles, rms


def interpolate_rms(
    rms_lookup, rot_angles_lookup, orig_rm, rot_angle,
    start_freq_mhz, end_freq_mhz
):

    orig_rot_angle = get_effective_rotation_angles(
        orig_rm, start_freq_mhz, end_freq_mhz
    )
    rms_lookup_use = np.copy(rms_lookup)
    rot_angles_lookup_use = np.copy(rot_angles_lookup)
    rot_angle_use = rot_angle

    # Correct for the branch cut
    orig_ind = np.min(np.where(rms_lookup_use > orig_rm))
    low_angle_inds = np.where(
        rot_angles_lookup_use[orig_ind:] < orig_rot_angle
    )[0]
    if np.shape(low_angle_inds)[0] > 0:
        rot_angles_lookup_use[
            (np.min(low_angle_inds)+orig_ind):
        ] += 2*np.pi
    high_angle_inds = np.where(
        rot_angles_lookup_use[:orig_ind] > orig_rot_angle
    )[0]
    if np.shape(high_angle_inds)[0] > 0:
        rot_angles_lookup_use[
            :(np.max(high_angle_inds)+1)
        ] -= 2*np.pi

    # Cut down the search range
    min_use_index = np.max(np.concatenate((np.where(
        rot_angles_lookup_use[:orig_ind] < orig_rot_angle-np.pi
    )[0], [0])))
    max_use_index = np.min(np.concatenate((np.where(
        rot_angles_lookup_use[orig_ind:] > orig_rot_angle+np.pi
    )[0]+orig_ind, [np.shape(rot_angles_lookup_use)[0]-1])))
    rms_lookup_use = rms_lookup_use[min_use_index:max_use_index+1]
    rot_angles_lookup_use = rot_angles_lookup_use[min_use_index:max_use_index+1]

    # Modify rot_angle to be in the correct range
    if rot_angle_use > orig_rot_angle+np.pi:
        rot_angle_use -= 2*np.pi
    elif rot_angle_use < orig_rot_angle-np.pi:
        rot_angle_use += 2*np.pi

    interp_func = scipy.interpolate.interp1d(
        rot_angles_lookup_use,
        rms_lookup_use,
        kind='cubic', bounds_error=True
    )
    print 'RM: {}'.format(orig_rm)
    print 'rot angle: {}'.format(rot_angle_use)
    interp_rm_val = interp_func(rot_angle_use)
    return interp_rm_val


def main():

    n_iter = 5
    step_size = .7

    #obs_list_1 = [obs_list_1[0]]
    #obs_list_2 = []

    rm_file = '/Users/rubybyrne/diffuse_survey_rm_tot.csv'
    #rm_outfile = '/Users/rubybyrne/diffuse_survey_rm_empirical.csv'
    rm_outfile = '/Users/rubybyrne/diffuse_survey_rm_empirical_Jul2020.csv'
    start_freq_mhz = 167.
    end_freq_mhz = 198.

    # Get RMs
    rm_data = np.genfromtxt(
        rm_file, delimiter=',', dtype=None, names=True, encoding=None
    )
    rms_orig = np.array([
        rm_data['RM'][np.where(rm_data['ObsID'] == int(obsid))][0] for obsid in obs_list_1+obs_list_2
    ])

    # Create lookup table:
    rot_angles_lookup, rms_lookup = create_rm_lookup_table(
        start_freq_mhz, end_freq_mhz
    )

    rms_use = np.copy(rms_orig)
    for iter_ind in range(n_iter):
        rot_angle_deltas_list = np.zeros(len(obs_list_1)+len(obs_list_2))
        rot_angle_list = np.zeros(len(obs_list_1)+len(obs_list_2))

        # Create average maps
        if iter_ind == 0: # Use saved maps
            q_average_map_path = '/Users/rubybyrne/diffuse_survey_plotting_May2020/StokesQ_average_map_more_obs.fits'
            u_average_map_path = '/Users/rubybyrne/diffuse_survey_plotting_May2020/StokesU_average_map_more_obs.fits'
            q_average_map = healpix_utils.load_map(q_average_map_path)
            u_average_map = healpix_utils.load_map(u_average_map_path)
        else: # Recalculate average maps with new RM values
            combined_maps, weight_maps = healpix_utils.average_healpix_maps(
                ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020',
                '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'],
                obs_lists = [obs_list_1, obs_list_2],
                nside=128,
                cube_names=['Residual_I', 'Residual_Q', 'Residual_U', 'Residual_V'],
                weighting='weighted',
                apply_radial_weighting=True,
                apply_rm_correction=True,
                use_rms=rms_use
            )
            q_average_map = combined_maps[1]
            u_average_map = combined_maps[2]
        q_average_map.explicit_to_implicit_ordering()
        u_average_map.explicit_to_implicit_ordering()

        # Calculate the empirical rotation angles
        for obsind, obsid in enumerate(obs_list_1+obs_list_2):
            if obsid in obs_list_1:
                path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Feb2020'
            elif obsid in obs_list_2:
                path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_baseline_cut_optimal_weighting_Mar2020'

            q_map = healpix_utils.load_map(
                '{}/output_data/{}_weighted_Residual_Q_HEALPix.fits'.format(path, obsid),
                quiet=True
            )
            u_map = healpix_utils.load_map(
                '{}/output_data/{}_weighted_Residual_Q_HEALPix.fits'.format(path, obsid),
                quiet=True
            )

            # Apply RM correction to maps
            maps_rot = healpix_utils.rm_correction(
                obsid, [None, q_map, u_map, None], rm_file=None,
                start_freq_mhz=start_freq_mhz, end_freq_mhz=end_freq_mhz,
                use_single_freq_calc=False, use_rm=rms_use[obsind]
            )
            q_map_rot = maps_rot[1]
            u_map_rot = maps_rot[2]
            # Apply additional rotation
            q_map_rot_signal = (
                np.cos(2*rot_angle_list[obsind])*q_map_rot.signal_arr
                + np.sin(2*rot_angle_list[obsind])*u_map_rot.signal_arr
            )
            u_map_rot_signal = (
                -np.sin(2*rot_angle_list[obsind])*q_map_rot.signal_arr
                + np.cos(2*rot_angle_list[obsind])*u_map_rot.signal_arr
            )
            q_map_rot.signal_arr = q_map_rot_signal
            u_map_rot.signal_arr = u_map_rot_signal

            # Confirm that pixel ordering matches
            if np.sum(np.abs(q_map_rot.pix_arr-u_map_rot.pix_arr)) != 0:
                print 'ERROR: Different pixel ordering.'

            q_average_map_signal = np.array(
                [q_average_map.signal_arr[ind] for ind in q_map_rot.pix_arr]
            )
            u_average_map_signal = np.array(
                [u_average_map.signal_arr[ind] for ind in q_map_rot.pix_arr]
            )

            #Get radial weighting
            obs_struct = scipy.io.readsav(
                '{}/metadata/{}_obs.sav'.format(path, obsid)
            )['obs']
            obs_vec = hp.pixelfunc.ang2vec(
                float(obs_struct['obsra']), float(obs_struct['obsdec']),
                lonlat=True
            )
            rad_weights = np.ones(np.shape(q_map_rot.pix_arr)[0])
            for pixind, pix in enumerate(q_map_rot.pix_arr):
                pix_vec = hp.pix2vec(q_map_rot.nside, pix, nest=q_map_rot.nest)
                rad_weights[pixind] = healpix_utils.obs_radial_weighting_function(
                    hp.rotator.angdist(pix_vec, obs_vec)*180./np.pi
                )

            tangent_numerator = np.sum(rad_weights*(
                q_average_map_signal*u_map_rot.signal_arr - u_average_map_signal*q_map_rot.signal_arr
            ))
            tangent_denominator = np.sum(rad_weights*(
                q_average_map_signal*q_map_rot.signal_arr + u_average_map_signal*u_map_rot.signal_arr
            ))
            rot_angle_delta = np.arctan2(tangent_numerator, tangent_denominator)
            rot_angle_deltas_list[obsind] = rot_angle_delta

        # Ensure that the change in the rotation angles is mean-zero
        mean_angle = np.arctan2(
            np.sum(np.sin(rot_angle_deltas_list)),
            np.sum(np.cos(rot_angle_deltas_list))
        )
        rot_angle_deltas_list = rot_angle_deltas_list - mean_angle
        rot_angle_list += step_size*rot_angle_deltas_list

        eff_rot_angle_start = get_effective_rotation_angles(
            rms_use, start_freq_mhz, end_freq_mhz
        )
        eff_rot_angle = eff_rot_angle_start + rot_angle_list
        # Ensure that the rotation angles are within +/- pi
        eff_rot_angle = np.arctan2(np.sin(eff_rot_angle), np.cos(eff_rot_angle))
        print eff_rot_angle

        # Convert effective rotation angles to RMs
        for obsind in range(len(obs_list_1)+len(obs_list_2)):
            rms_use[obsind] = interpolate_rms(
                rms_lookup, rot_angles_lookup,
                rms_orig[obsind], eff_rot_angle[obsind],
                start_freq_mhz, end_freq_mhz
            )

    # Save RMs to a CSV file
    csv_outfile = open(rm_outfile, 'w')
    outfile_writer = csv.writer(csv_outfile)
    outfile_writer.writerow(['ObsID', 'RM'])
    for obsind, obsid in enumerate(obs_list_1+obs_list_2):
        outfile_writer.writerow([obsid, new_rm_vals[obsind]])
    csv_outfile.close()


if __name__=='__main__':
    main()
