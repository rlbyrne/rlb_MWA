#!/usr/bin/python

from astropy.io import fits
import scipy.io
import numpy as np
import healpy as hp
import sys
import math
import os


class HealpixMap:

    def __init__(self, signal_arr, pix_arr, nside, nest=False,
                 coords='equitorial', quiet=False):
        if len(pix_arr) == 0:
            if not quiet:
                print 'No pixel values supplied. Assuming implicit ordering.'
        else:
            if len(signal_arr) != len(pix_arr):
                print 'ERROR: Pixel index and data lengths do not match. Exiting.'
                sys.exit(1)
        if nest == 'nested' or nest == 'nest' or nest is True:
            nest = True
        elif nest == 'ring' or nest is False:
            nest = False
        else:
            print 'ERROR: Invalid nest parameter. Exiting.'
            sys.exit(1)
        self.signal_arr = np.array(signal_arr)
        self.pix_arr = np.array(pix_arr, dtype=int)
        self.nside = nside
        self.nest = nest
        self.coords = coords

    def get_ra_dec(self, ra_cut=270):
        if self.coords == 'galactic':
            theta_gal, phi_gal = hp.pixelfunc.pix2ang(
                self.nside, self.pix_arr, nest=self.nest
                )
            rot = hp.rotator.Rotator(coord=['G', 'C'])
            theta_eq, phi_eq = rot(theta_gal, phi_gal)
            ra_arr = phi_eq*180/math.pi
            dec_arr = 90.-theta_eq*180/math.pi
        elif self.coords == 'equitorial':
            ra_arr, dec_arr = hp.pixelfunc.pix2ang(
                self.nside, self.pix_arr, nest=self.nest, lonlat=True
                )
        else:
            print 'ERROR: Coordinates must be galactic or equitorial.'
            sys.exit(1)
        ra_arr[np.where(ra_arr > ra_cut)] = ra_arr[
            np.where(ra_arr > ra_cut)
        ] - 360.
        ra_arr[np.where(ra_arr < ra_cut-360.)] = ra_arr[
            np.where(ra_arr > ra_cut)
        ] + 360.
        self.ra_arr = ra_arr
        self.dec_arr = dec_arr

    def get_pixel_corners(self, ra_cut=270):
        corner_coords = hp.boundaries(
            self.nside, self.pix_arr, step=1, nest=self.nest
        )
        if self.coords == 'galactic':
            thetas_gal = []
            phis_gal = []
            for pixel in range(len(self.signal_arr)):
                thetas_gal_pix, phis_gal_pix = hp.pixelfunc.vec2ang(
                    np.transpose(corner_coords[pixel]), lonlat=False
                )
                thetas_gal.append(thetas_gal_pix)
                phis_gal.append(phis_gal_pix)
            thetas_gal = np.array(thetas_gal)
            phis_gal = np.array(phis_gal)
            rot = hp.rotator.Rotator(coord=['G', 'C'])
            thetas_eq, phis_eq = rot(thethas_gal, phis_gal)
            ras = phis_eq*180/math.pi
            decs = 90. - thetas_eq*180/math.pi
        elif self.coords == 'equitorial':
            ras = []
            decs = []
            for pixel in range(len(self.signal_arr)):
                ras_pix, decs_pix = hp.pixelfunc.vec2ang(
                    np.transpose(corner_coords[pixel]), lonlat=True
                )
                ras_pix[np.where(ras_pix > ra_cut)] = ras_pix[
                    np.where(ras_pix > ra_cut)
                ] - 360.
                ras_pix[np.where(ras_pix < ra_cut-360.)] = ras_pix[
                    np.where(ras_pix > ra_cut)
                ] + 360.
                ras_pix[np.where(ras_pix-np.amin(ras_pix) > 180.)] = ras_pix[
                    np.where(ras_pix-np.amin(ras_pix) > 180.)
                ] - 360.
                ras.append(ras_pix)
                decs.append(decs_pix)
            ras = np.array(ras)
            decs = np.array(decs)
        else:
            print 'ERROR: Coordinates must be galactic or equitorial.'
            sys.exit(1)
        self.pix_corner_ras_arr = ras
        self.pix_corner_decs_arr = decs

    def explicit_to_implicit_ordering(self):
        if len(self.pix_arr) > 0:
            signal_vals_implicit = np.full(12*self.nside**2, hp.pixelfunc.UNSEEN)
            signal_vals_implicit[self.pix_arr] = self.signal_arr
            self.signal_arr = signal_vals_implicit
            self.pix_arr = np.empty(0, dtype=int)
        elif len(self.signal_arr) != 12*self.nside**2:
            print 'ERROR: Wrong number of pixels for implicit ordering. Exiting.'
            sys.exit(1)

    def implicit_to_explicit_ordering(self):
        use_inds = np.where(self.signal_arr != hp.pixelfunc.UNSEEN)[0]
        if len(self.pix_arr) == 0:
            self.pix_arr = use_inds
        else:
            self.pix_arr = self.pix_arr[use_inds]
        self.signal_arr = self.signal_arr[use_inds]

    def reorder_nest_to_ring(self):
        if self.nest:
            self.explicit_to_implicit_ordering()
            self.signal_arr = hp.pixelfunc.reorder(filtered_map, n2r=True)
            self.nest = False
            self.implicit_to_explicit_ordering()

    def reorder_ring_to_nest(self):
        if not self.nest:
            self.explicit_to_implicit_ordering()
            self.signal_arr = hp.pixelfunc.reorder(self.signal_arr, r2n=True)
            self.nest = True
            self.implicit_to_explicit_ordering()

    def resample(self, nside, quiet=False):
        if not quiet:
            print 'Resampling map: nside {} to {}'.format(self.nside, nside)
        if self.nest:
            ordering = 'nested'
        else:
            ordering = 'ring'
        self.explicit_to_implicit_ordering()
        self.signal_arr = hp.pixelfunc.ud_grade(
            np.array(self.signal_arr), nside, pess=True, order_in=ordering
            )
        self.nside = nside
        self.implicit_to_explicit_ordering()

    def get_alm(self, lmax=None):
        self.explicit_to_implicit_ordering()
        if self.nest:
            signal_vals = hp.pixelfunc.reorder(self.signal_arr, n2r=True)
        else:
            signal_vals = self.signal_arr
        alm = hp.sphtfunc.map2alm(signal_vals, lmax=lmax)
        alm = np.copy(alm)
        lmax = hp.sphtfunc.Alm.getlmax(len(alm))
        l, m = hp.sphtfunc.Alm.getlm(lmax)
        lm = np.zeros([len(alm), 2])
        lm[:, 0] = l
        lm[:, 1] = m
        self.spherical_harmonics_amp = alm
        self.spherical_harmonics_lm = lm

    def filter_map(self, lmin=None, lmax=None, filter_width=0):
        # This function borrows from code by Miguel Morales
        print self.coords
        if lmax is not None:
            alm_limit = int(math.ceil(lmax + filter_width/2.))
        else:
            alm_limit = None
        self.get_alm(lmax=alm_limit)
        alm = self.spherical_harmonics_amp
        lm = self.spherical_harmonics_lm
        # Define tukey windowing function
        window_fn = np.ones(int(max(lm[:, 0])+1))
        if lmin is not None:
            for i in range(int(math.ceil(lmin-filter_width/2.))):
                window_fn[i] = 0
            for i in range(int(math.ceil(lmin-filter_width/2.)),
                           int(math.floor(lmin+filter_width/2.))+1
                           ):
                window_fn[i] = (
                    -math.cos((i-lmin+filter_width/2.)*math.pi/filter_width)+1
                    )/2.
        if lmax is not None:
            for i in range(int(math.ceil(lmax-filter_width/2.)),
                           int(math.floor(lmax+filter_width/2.))+1
                           ):
                window_fn[i] = (
                    math.cos((i-lmax+filter_width/2.)*math.pi/filter_width)+1
                    )/2.
            for i in range(int(math.floor(lmax+filter_width/2.))+1,
                           len(window_fn)
                           ):
                window_fn[i] = 0
        # Apply windowing function
        for i in range(len(window_fn)):
            alm[np.where(lm[:, 0] == i)] = alm[
                np.where(lm[:, 0] == i)
                ]*window_fn[i]
        filtered_map = HealpixMap(
            hp.sphtfunc.alm2map(alm, self.nside), np.array([]),
            nside=self.nside, nest=False, coords=self.coords
        )
        # If input was nested ordering, convert to nested
        if self.nest:
            filtered_map.reorder_ring_to_nest()
        else:
            filtered_map.implicit_to_explicit_ordering()
        filtered_map.implicit_to_explicit_ordering()
        self.signal_arr = filtered_map.signal_arr
        self.pix_arr = filtered_map.pix_arr

    def write_data_to_fits(self, save_filename):
        signal_column = fits.Column(
            name='SIGNAL',
            array=np.array(self.signal_arr),
            format='1E'
            )
        pixelnum_column = fits.Column(
            name='PIXEL',
            array=np.array(self.pix_arr),
            format='1J'
            )
        header = fits.Header()  # initialize header object
        header['nside'] = self.nside
        if self.nest:
            header['ordering'] = 'nested'
        else:
            header['ordering'] = 'ring'
        header['indxschm'] = 'explicit'
        hdu_0 = fits.PrimaryHDU()
        hdu_1 = fits.BinTableHDU.from_columns(
            [signal_column, pixelnum_column],
            header=header
            )
        hdu_list = fits.HDUList([hdu_0, hdu_1])
        print 'Saving data to {}'.format(save_filename)
        hdu_list.writeto(save_filename, overwrite=True)


def load_map(data_filename):
    # Load a HEALPix map formatted with FHD image conventions

    print 'Loading HEALPix map {}'.format(data_filename)
    contents = fits.open(data_filename)
    nside = int(contents[1].header['nside'])
    ordering = contents[1].header['ordering']
    data = contents[1].data
    contents.close()

    signal_vals = data.field('SIGNAL')
    pixel_vals = data.field('PIXEL')

    if ordering.lower() == 'ring':
        nest = False
    elif ordering.lower() == 'nested':
        nest = True
    else:
        print 'ERROR: Invalid ordering parameter.'
        print 'Ordering must be "ring" or "nested". Exiting.'
        sys.exit(1)

    healpix_map = HealpixMap(
        signal_vals, pixel_vals, nside, nest=nest, coords='equitorial'
        )
    return healpix_map


def load_global_map(data_filename):
    # Load a HEALPix map formatted with Haslam-style conventions

    contents = fits.open(data_filename)
    nside = int(contents[1].header['nside'])
    ordering = contents[1].header['ordering']
    data = contents[1].data
    contents.close()

    if ordering.lower() == 'ring':
        nest = False
    elif ordering.lower() == 'nested':
        nest = True
    else:
        print 'ERROR: Invalid ordering parameter.'
        print 'Ordering must be "ring" or "nested". Exiting.'
        sys.exit(1)

    signal_vals = data.field('TEMPERATURE')
    healpix_map = HealpixMap(
        signal_vals, [], nside, nest=nest, coords='galactic'
        )
    healpix_map.implicit_to_explicit_ordering()
    return healpix_map


def load_fhd_output_map(data_filename, cube='model', freq_index=0):
    # Load a HEALPix map formatted with FHD-for-eppsilon conventions

    if (
        cube != 'model' and cube != 'data'
        and cube != 'variance' and cube != 'weights'
    ):
        print 'ERROR: Invalid cube option.'
        print 'Cube must be "model", "dirty", "variance", or "weights". Exiting.'
    if cube == 'model':
        cube_name = 'model_cube'
    if cube == 'dirty':
        cube_name = 'dirty_cube'
    if cube == 'variance':
        cube_name = 'variance_cube'
    if cube == 'weights':
        cube_name = 'weights_cube'

    data = scipy.io.readsav(data_filename)
    cube_data = data[cube_name][freq_index]
    nside = int(data['nside'])
    hpx_inds = [int(val) for val in data['hpx_inds']]
    healpix_map = HealpixMap(
        cube_data, hpx_inds, nside, nest=False, coords='equitorial'
        )
    return healpix_map


def difference_healpix_maps(map1, map2):

    if map1.nside != map2.nside:
        print 'ERROR: Healpix map nsides do not match. Exiting.'
        sys.exit(1)
    if map1.nest != map2.nest:
        print 'WARNING: Healpix map orderings do not match. Reordering.'
        if map1.nest:
            map1.reorder_nest_to_ring()
        else:
            map2.reorder_nest_to_ring()
    if map1.pix_arr == map2.pix_arr:
        data_diff = list(np.array(map1.signal_arr) - np.array(map2.signal_arr))
        diff_map = HealpixMap(
            data_diff, map1.pix_arr, map1.nside, nest=map1.nest
            )
    else:
        pixelnum_use = list(set(map1.pix_arr).intersection(map2.pix_arr))
        data_diff = [
            map1.signal_arr[map1.pix_arr.index(pixelnum)]
            - map2.signal_arr[map2.pix_arr.index(pixelnum)]
            for pixelnum in pixelnum_use
            ]
        diff_map = HealpixMap(
            data_diff, pixelnum_use, map1.nside, nest=map1.nest
            )
    return diff_map


def multiply_healpix_maps(map1, map2):

    if map1.nside != map2.nside:
        print 'ERROR: Healpix map nsides do not match. Exiting.'
        sys.exit(1)
    if map1.nest != map2.nest:
        print 'WARNING: Healpix map orderings do not match. Reordering.'
        if map1.nest:
            map1.reorder_nest_to_ring()
        else:
            map2.reorder_nest_to_ring()
    if map1.pix_arr == map2.pix_arr:
        data_prod = list(np.array(map1.signal_arr) * np.array(map2.signal_arr))
        prod_map = HealpixMap(
            data_prod, map1.pix_arr, map1.nside, nest=map1.nest
            )
    else:
        pixelnum_use = list(set(map1.pix_arr).intersection(map2.pix_arr))
        data_prod = [
            map1.signal_arr[map1.pix_arr.index(pixelnum)]
            * map2.signal_arr[map2.pix_arr.index(pixelnum)]
            for pixelnum in pixelnum_use
            ]
        prod_map = HealpixMap(
            data_prod, pixelnum_use, map1.nside, nest=map1.nest
            )
    return prod_map


def divide_healpix_maps(map1, map2):

    if map1.nside != map2.nside:
        print 'ERROR: Healpix map nsides do not match. Exiting.'
        sys.exit(1)
    if map1.nest != map2.nest:
        print 'WARNING: Healpix map orderings do not match. Reordering.'
        if map1.nest:
            map1.reorder_nest_to_ring()
        else:
            map2.reorder_nest_to_ring()
    if map1.pix_arr == map2.pix_arr:
        data_div = list(np.array(map1.signal_arr) / np.array(map2.signal_arr))
        div_map = HealpixMap(
            data_div, map1.pix_arr, map1.nside, nest=map1.nest
            )
    else:
        pixelnum_use = list(set(map1.pix_arr).intersection(map2.pix_arr))
        data_div = [
            map1.signal_arr[map1.pix_arr.index(pixelnum)]
            / map2.signal_arr[map2.pix_arr.index(pixelnum)]
            for pixelnum in pixelnum_use
            ]
        div_map = HealpixMap(
            data_div, pixelnum_use, map1.nside, nest=map1.nest
            )
    return div_map


def average_healpix_maps_simple(maps_arr):
    # This function averages a set of preloaded maps.
    # The function returns the averaged map, the variance at each pixel, and
    # the number of samples that contribute to each pixel.

    nside = maps_arr[0].nside
    nest = maps_arr[0].nest
    for map in maps_arr:
        if map.nside != nside:
            print 'ERROR: Healpix map nsides do not match. Exiting.'
            sys.exit(1)
        if map.nest != nest:
            print 'ERROR: Healpix map orderings do not match. Exiting.'
            sys.exit(1)
    pixelnum_use = list(set(
        [pixelnum for map in maps_arr for pixelnum in map.pix_arr]
        ))
    pixel_vals = [[] for i in range(len(pixelnum_use))]
    for map in maps_arr:
        for i, pixelnum in enumerate(pixelnum_use):
            if pixelnum in map.pix_arr:
                pixel_vals[i].append(
                    map.signal_arr[map.pix_arr.index(pixelnum)]
                    )

    ave_map = HealpixMap(
        [np.average(signal_vals) for signal_vals in pixel_vals], pixelnum_use,
        nside, nest=nest
        )
    var_map = HealpixMap(
        [np.var(signal_vals) for signal_vals in pixel_vals], pixelnum_use,
        nside, nest=nest
        )
    nsamples_map = HealpixMap(
        [len(signal_vals) for signal_vals in pixel_vals], pixelnum_use,
        nside, nest=nest
        )

    return ave_map, var_map, nsamples_map


def average_healpix_maps(
    fhd_run_path, obs_list, obs_weights=None, nside=None,
    cube_names=['Residual_I'], apply_radial_weighting=False
):

    if obs_weights is None or len(obs_weights) == 0:
        print 'Observation weights not provided. Using equal weighting.'
        obs_weights = np.array([1.]*len(obs_list))
    else:
        if len(obs_weights) != len(obs_list):
            print 'ERROR: Invalid obs_weights. Exiting.'
            sys.exit(1)
        obs_weights = np.array(obs_weights)
    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
    exclude_obs_list = []
    for obs in obs_list:
        for cube in cube_names:
            if '{}_uniform_{}_HEALPix.fits'.format(obs, cube) not in data_files:
                print 'WARNING: File {}/output_data/{}_uniform_{}_HEALPix.fits not found. Excluding observation {} from the average.'.format(
                    fhd_run_path, obs, cube, obs
                )
                exclude_obs_list.append(obs)
                continue
    obs_weights = np.array([
        obs_weights[ind] for ind in range(len(obs_weights))
        if obs_list[ind] not in exclude_obs_list
    ])
    obs_list = [obs for obs in obs_list if obs not in exclude_obs_list]

    print 'Averaging {} observations'.format(len(obs_list))

    for obs_ind, obsid in enumerate(obs_list):
        maps = []
        for cube_ind, cube in enumerate(cube_names):
            map = load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
                fhd_run_path, obsid, cube
            ))
            if obs_ind == 0 and cube_ind == 0:  # Use first map for conventions
                if nside is None:
                    nside = map.nside
                    print 'nside is not specified. Using nside {}.'.format(nside)
                else:
                    if map.nside != nside:
                        map.resample(nside)
                nest = map.nest
                coords = map.coords
                # Initialize output arrays
                signal_array = np.full(
                    (len(cube_names), 12*nside**2), hp.pixelfunc.UNSEEN
                )
                weights_array = np.full(12*nside**2, hp.pixelfunc.UNSEEN)
            else:
                if map.nside != nside:
                    map.resample(nside)
                if map.nest != nest:
                    print 'Map nesting conventions do not match. Converting.'
                    if nest:
                        map.reorder_ring_to_nest()
                    else:
                        map.reorder_nest_to_ring()
                if map.coords != coords:
                    print 'ERROR: Map coordinates do not match. Exiting.'
                    sys.exit(1)
            maps.append(map)

        for pix in maps[0].pix_arr:  # Use the first cube for the pixel list
            if apply_radial_weighting:
                obs_struct = scipy.io.readsav(
                    '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
                )['obs']
                obs_vec = hp.pixelfunc.ang2vec(
                    float(obs_struct['obsra']), float(obs_struct['obsdec']),
                    lonlat=True
                )
                pix_vec = hp.pix2vec(nside, pix, nest=nest)
                rad_weight = obs_radial_weighting_function(
                    (pix_vec[0]-obs_vec[0])**2.
                    + (pix_vec[1]-obs_vec[1])**2.
                    + (pix_vec[2]-obs_vec[2])**2.
                )
                use_weight = rad_weight*obs_weights[obs_ind]
            else:
                use_weight = obs_weights[obs_ind]
            if signal_array[0, pix] == hp.pixelfunc.UNSEEN:
                for cube_ind in range(len(cube_names)):
                    signal_array[cube_ind, pix] = use_weight*maps[
                        cube_ind
                    ].signal_arr[
                        np.where(maps[cube_ind].pix_arr == pix)
                    ]
                    weights_array[pix] = use_weight
            else:
                for cube_ind in range(len(cube_names)):
                    signal_array[cube_ind, pix] += use_weight*maps[
                        cube_ind
                    ].signal_arr[
                        np.where(maps[cube_ind].pix_arr == pix)
                    ]
                    weights_array[pix] += use_weight

    weighted_ave_signal_array = signal_array/weights_array[None, :]
    weighted_ave_signal_array[
        :, np.where(weights_array == hp.pixelfunc.UNSEEN)
    ] = hp.pixelfunc.UNSEEN

    averaged_maps = []
    for cube_ind in range(len(cube_names)):
        map = HealpixMap(
            weighted_ave_signal_array[cube_ind, :], [], nside, nest=nest,
            coords=coords, quiet=True
        )
        map.implicit_to_explicit_ordering()
        averaged_maps.append(map)

    weights_map = HealpixMap(
        weights_array, [], nside, nest=nest, coords=coords, quiet=True
    )
    weights_map.implicit_to_explicit_ordering()

    return averaged_maps, weights_map


def calculate_variance_healpix_maps(
    fhd_run_path, obs_list, saved_averaged_maps=None, obs_weights=None,
    nside=None, cube_names=['Residual_I'], apply_radial_weighting=False
):

    if obs_weights is None or len(obs_weights) == 0:
        print 'Observation weights not provided. Using equal weighting.'
        obs_weights = np.array([1.]*len(obs_list))
    else:
        if len(obs_weights) != len(obs_list):
            print 'ERROR: Invalid obs_weights. Exiting.'
            sys.exit(1)
        obs_weights = np.array(obs_weights)
    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
    exclude_obs_list = []
    for obs in obs_list:
        for cube in cube_names:
            if '{}_uniform_{}_HEALPix.fits'.format(obs, cube) not in data_files:
                print 'WARNING: File {}/output_data/{}_uniform_{}_HEALPix.fits not found. Excluding observation {} from the average.'.format(
                    fhd_run_path, obs, cube, obs
                )
                exclude_obs_list.append(obs)
                continue
    obs_weights = np.array([
        obs_weights[ind] for ind in range(len(obs_weights))
        if obs_list[ind] not in exclude_obs_list
    ])
    obs_list = [obs for obs in obs_list if obs not in exclude_obs_list]

    if saved_averaged_maps is None:
        averaged_maps, null = average_healpix_maps(
            fhd_run_path, obs_list, obs_weights=obs_weights, nside=nside,
            cube_names=cube_names,
            apply_radial_weighting=apply_radial_weighting
        )
    else:
        if len(saved_averaged_maps) != len(cube_names):
            print 'WARNING: Invalid saved_averaged_maps paths provided. Recalculating data average.'
            averaged_maps, null = average_healpix_maps(
                fhd_run_path, obs_list, obs_weights=obs_weights, nside=nside,
                cube_names=cube_names,
                apply_radial_weighting=apply_radial_weighting
            )
        else:
            print 'Loading saved data averages.'
            averaged_maps = []
            for cube_ind in range(len(cube_names)):
                map = load_map(saved_averaged_maps[cube_ind])
                averaged_maps.append(map)
    for cube_ind in range(len(cube_names)):
        averaged_maps[cube_ind].explicit_to_implicit_ordering()

    print 'Calculating variance for {} observations'.format(len(obs_list))

    for obs_ind, obsid in enumerate(obs_list):
        maps = []
        for cube_ind, cube in enumerate(cube_names):
            map = load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
                fhd_run_path, obsid, cube
            ))
            if obs_ind == 0 and cube_ind == 0:  # Use first map for conventions
                if nside is None:
                    nside = map.nside
                    print 'nside is not specified. Using nside {}.'.format(nside)
                else:
                    if map.nside != nside:
                        map.resample(nside)
                nest = map.nest
                coords = map.coords
                # Initialize output arrays
                var_array = np.full(
                    (len(cube_names), 12*nside**2), hp.pixelfunc.UNSEEN
                )
                weights_array = np.full(12*nside**2, hp.pixelfunc.UNSEEN)
                nsamples_array = np.full(12*nside**2, hp.pixelfunc.UNSEEN)
            else:
                if map.nside != nside:
                    map.resample(nside)
                if map.nest != nest:
                    print 'Map nesting conventions do not match. Converting.'
                    if nest:
                        map.reorder_ring_to_nest()
                    else:
                        map.reorder_nest_to_ring()
                if map.coords != coords:
                    print 'ERROR: Map coordinates do not match. Exiting.'
                    sys.exit(1)
            maps.append(map)

        for pix in maps[0].pix_arr:  # Use the first cube for the pixel list
            if apply_radial_weighting:
                obs_struct = scipy.io.readsav(
                    '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
                )['obs']
                obs_vec = hp.pixelfunc.ang2vec(
                    float(obs_struct['obsra']), float(obs_struct['obsdec']),
                    lonlat=True
                )
                pix_vec = hp.pix2vec(nside, pix, nest=nest)
                rad_weight = obs_radial_weighting_function(
                    (pix_vec[0]-obs_vec[0])**2.
                    + (pix_vec[1]-obs_vec[1])**2.
                    + (pix_vec[2]-obs_vec[2])**2.
                )
                use_weight = rad_weight*obs_weights[obs_ind]
            else:
                use_weight = obs_weights[obs_ind]
            if nsamples_array[pix] == hp.pixelfunc.UNSEEN:
                for cube_ind in range(len(cube_names)):
                    var_array[cube_ind, pix] = use_weight*(
                        maps[cube_ind].signal_arr[
                            np.where(maps[cube_ind].pix_arr == pix)
                        ][0] - averaged_maps[cube_ind].signal_arr[pix]
                    )**2
                    weights_array[pix] = use_weight
                    nsamples_array[pix] = 1
            else:
                for cube_ind in range(len(cube_names)):
                    var_array[cube_ind, pix] += use_weight*(
                        maps[cube_ind].signal_arr[
                            np.where(maps[cube_ind].pix_arr == pix)
                        ] - averaged_maps[cube_ind].signal_arr[pix]
                    )**2
                    weights_array[pix] += use_weight
                    nsamples_array[pix] += 1

    weighted_ave_var_array = var_array/weights_array[None, :]

    weighted_ave_var_array[
        :, np.where(weights_array == hp.pixelfunc.UNSEEN)
    ] = hp.pixelfunc.UNSEEN
    weighted_ave_var_array[
        :, np.where(nsamples_array == 1)
    ] = 0

    variance_maps = []
    for cube_ind in range(len(cube_names)):
        map = HealpixMap(
            weighted_ave_var_array[cube_ind, :], [], nside, nest=nest,
            coords=coords, quiet=True
        )
        map.implicit_to_explicit_ordering()
        variance_maps.append(map)
        averaged_maps[cube_ind].implicit_to_explicit_ordering()

    weights_map = HealpixMap(
        weights_array, [], nside, nest=nest, coords=coords, quiet=True
    )
    weights_map.implicit_to_explicit_ordering()

    nsamples_map = HealpixMap(
        nsamples_array, [], nside, nest=nest, coords=coords, quiet=True
    )
    nsamples_map.implicit_to_explicit_ordering()

    return variance_maps, averaged_maps, weights_map, nsamples_map


def obs_radial_weighting_function(dist, max_dist=14., taper_width=4.):

    if dist < max_dist-taper_width:
        weight = 1.
    elif dist > max_dist:
        weight = 0.
    else:
        weight = 0.5*np.cos(np.pi*(dist+taper_width-max_dist)/taper_width)+0.5
    return weight


def combine_maps_nearest_data(
    fhd_run_path, obs_list_file=None, nside=None, cube_names=['Residual_I']
):
    # This function loads a set of maps from an FHD output directory (or a list
    # of observations) and returns a combined map.
    # Each pixel of the combined map is based on the closest observation to
    # that pixel.

    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    if obs_list_file is None:  # use all obs in the data directory
        data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
        data_files = [
            file for file in data_files
            if '_uniform_{}_HEALPix.fits'.format(cube_names[0]) in file
        ]
        obs_list = [file[0:10] for file in data_files]
    else:  # use the obs file list
        obs_list = open(obs_list_file, 'r').readlines()
        # strip newline characters
        obs_list = [obs.strip() for obs in obs_list]
        # remove duplicates
        obs_list = list(set(obs_list))

    print 'Combining {} observations'.format(len(obs_list))

    obs_centers = []
    for obsid in obs_list:
        obs_struct = scipy.io.readsav(
            '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
        )['obs']
        obs_vec = hp.pixelfunc.ang2vec(
            float(obs_struct['obsra']), float(obs_struct['obsdec']),
            lonlat=True
        )
        obs_centers.append(obs_vec)

    for obs_ind, obsid in enumerate(obs_list):
        maps = []
        for cube_ind, cube in enumerate(cube_names):
            map = load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
                fhd_run_path, obsid, cube
            ))
            if obs_ind == 0 and cube_ind == 0:  # Use first map for conventions
                if nside is None:
                    nside = map.nside
                    print 'nside is not specified. Using nside {}.'.format(nside)
                else:
                    if map.nside != nside:
                        map.resample(nside)
                nest = map.nest
                coords = map.coords
                # Initialize output arrays
                signal_array = np.full(
                    (len(cube_names), 12*nside**2), hp.pixelfunc.UNSEEN
                )
                pixels_used = np.array([], dtype=int)
            else:
                if map.nside != nside:
                    map.resample(nside)
                if map.nest != nest:
                    print 'Map nesting conventions do not match. Converting.'
                    if nest:
                        map.reorder_ring_to_nest()
                    else:
                        map.reorder_nest_to_ring()
                if map.coords != coords:
                    print 'ERROR: Map coordinates do not match. Exiting.'
                    sys.exit(1)
            map.filter_map(lmin=None, lmax=10, filter_width=2)
            maps.append(map)

        # Use the first cube for the pixel list (assume cube pixels match)
        overlapping_pix = np.intersect1d(
            maps[0].pix_arr, pixels_used, assume_unique=True
        )
        unique_pix = np.setdiff1d(
            maps[0].pix_arr, overlapping_pix, assume_unique=True
        )
        pixels_used = np.append(pixels_used, unique_pix)

        for pix in overlapping_pix:
            vec = hp.pix2vec(nside, pix, nest=nest)
            distances = [
                (vec[0]-obs_centers[obs_ind][0])**2.
                + (vec[1]-obs_centers[obs_ind][1])**2.
                + (vec[2]-obs_centers[obs_ind][2])**2.
                for obs_ind in range(obs_ind+1)
            ]
            if distances.index(min(distances)) == obs_ind:
                for cube_ind in range(len(cube_names)):
                    signal_array[cube_ind, pix] = maps[cube_ind].signal_arr[
                        np.where(maps[cube_ind].pix_arr == pix)
                    ]
        for pix in unique_pix:
            for cube_ind in range(len(cube_names)):
                signal_array[cube_ind, pix] = maps[cube_ind].signal_arr[
                    np.where(maps[cube_ind].pix_arr == pix)[0]
                ]

    combined_maps = []
    for cube_ind in range(len(cube_names)):
        map = HealpixMap(
            signal_array[cube_ind, :], [], nside, nest=nest, coords=coords,
            quiet=True
        )
        map.implicit_to_explicit_ordering()
        combined_maps.append(map)

    return combined_maps
