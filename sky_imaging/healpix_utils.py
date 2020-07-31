#!/usr/bin/python

from astropy.io import fits
import scipy.io
import scipy
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
        if (
            coords is not None and coords != 'equitorial'
            and coords != 'galactic' and coords != 'ecliptic'
        ):
            print 'ERROR: Invalid coords parameter.'
            print 'Valid options are equitorial, galactic, or ecliptic.'
            print 'Exiting.'
            sys.exit(1)
        self.signal_arr = np.array(signal_arr)
        self.pix_arr = np.array(pix_arr, dtype=int)
        self.nside = nside
        self.nest = nest
        self.coords = coords
        if self.coords == 'galactic':
            self.coords_healpy_conv = 'G'
        elif self.coords == 'equitorial':
            self.coords_healpy_conv = 'C'
        elif self.coords == 'ecliptic':
            self.coords_healpy_conv = 'E'
        else:
            self.coords_healpy_conv = None

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
            np.where(ra_arr < ra_cut-360.)
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
                where_above = np.where(ras_pix > ra_cut)
                if len(where_above[0]) > 0:
                    ras_pix[where_above] = ras_pix[where_above] - 360.
                where_below = np.where(ras_pix < ra_cut-360.)
                if len(where_below[0]) > 0:
                    ras_pix[where_below] = ras_pix[where_below] + 360.
                where_separated = np.where(ras_pix-np.amin(ras_pix) > 180.)
                if len(where_separated[0]) > 0:
                    ras_pix[where_separated] = ras_pix[where_separated] - 360.
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
        self.implicit_to_explicit_ordering()

    def filter_map(
        self, lmin=None, lmax=None, filter_width=0, return_full_map=False
    ):
        # This function borrows from code by Miguel Morales
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
            nside=self.nside, nest=False, coords=self.coords, quiet=True
        )
        # If input was nested ordering, convert to nested
        if self.nest:
            filtered_map.reorder_ring_to_nest()
        # Use only pixels from the initial map if return_full_map is unset
        if not return_full_map:
            filtered_map.explicit_to_implicit_ordering()
            for pix_ind, pixel in enumerate(self.pix_arr):
                self.signal_arr[pix_ind] = filtered_map.signal_arr[pixel]
        else:
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


def load_map(data_filename, quiet=False):
    # Load a HEALPix map formatted with FHD image conventions

    if not quiet:
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
    hpx_inds = data['hpx_inds']
    hpx_inds.astype(int)
    healpix_map = HealpixMap(
        cube_data, hpx_inds, nside, nest=False, coords='equitorial'
    )
    return healpix_map


def load_fhd_input_map(data_filename, cube_ind=0):
    # Load a HEALPix map formatted in a .sav file to be input to FHD as a model

    data = scipy.io.readsav(data_filename, python_dict=True)
    map_data = data['model_arr']
    if np.size(map_data[0]) > 1: # array contains multiple maps
        map_data = map_data[cube_ind]
    nside = int(data['nside'])
    hpx_inds = data['hpx_inds']
    hpx_inds.astype(int)
    healpix_map = HealpixMap(
        map_data, hpx_inds, nside, nest=False, coords='equitorial'
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
    if np.max(np.abs(map1.pix_arr-map2.pix_arr)) == 0:
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

    pixelnum_use = np.intersect1d(map1.pix_arr, map2.pix_arr)
    denominator = np.squeeze(np.array([
        map2.signal_arr[np.where(map2.pix_arr==pixelnum)[0]]
        for pixelnum in pixelnum_use
    ]))
    pixelnum_use = pixelnum_use[np.where(denominator != 0.)[0]]
    data_div = np.squeeze(np.array([
        map1.signal_arr[np.where(map1.pix_arr==pixelnum)[0]]
        / map2.signal_arr[np.where(map2.pix_arr==pixelnum)[0]]
        for pixelnum in pixelnum_use
    ]))
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
    fhd_run_paths, obs_lists=None, obs_weights_list=None, nside=None,
    cube_names=['Residual_I'], weighting='uniform', apply_radial_weighting=False,
    apply_rm_correction=False,
    use_rms=None, # List of RMs, only used if apply_rm_correction=True
    rm_file=None, # File for RM lookup, only used if apply_rm_correction=True and use_rms=None
    quiet=False
):

    if isinstance(fhd_run_paths, str): #check if string
        fhd_run_paths = [fhd_run_paths]
    if obs_lists is not None:
        if len(fhd_run_paths) == 1 and isinstance(obs_lists[0], str):
            obs_lists = [obs_lists]
        elif len(fhd_run_paths) != len(obs_lists):
            print 'ERROR: obs_list not provided for each path'
            sys.exit(1)
    if obs_weights_list is not None:
        if len(fhd_run_paths) == 1 and isinstance(obs_lists[0], str):
            obs_weights_list = [obs_weights_list]
        elif len(fhd_run_paths) != len(obs_weights_list):
            print 'ERROR: obs_weights_list not provided for each path'
            sys.exit(1)

    for path_ind, fhd_run_path in enumerate(fhd_run_paths):

        if obs_lists is not None:
            obs_list = obs_lists[path_ind]
        else:
            obs_list = None
        if obs_weights_list is not None:
            obs_weights = obs_weights_list[path_ind]
        else:
            obs_weights = None

        if fhd_run_path[-1] == '/':
            fhd_run_path = fhd_run_path[:-1]

        if obs_list is None:  # use all obs in the data directory
            if obs_weights is not None:
                print 'WARNING: No obs_list provided. Disregarding obs_weights and using equal weighting.'
                obs_weights = None
            data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
            for cube_ind, cube in enumerate(cube_names):
                data_files_cube = [
                    file for file in data_files
                    if '_{}_{}_HEALPix.fits'.format(weighting, cube) in file
                ]
                print data_files
                print '_{}_{}_HEALPix.fits'.format(weighting, cube)
                obs_list_cube = [file[0:10] for file in data_files_cube]
                if cube_ind == 0:
                    obs_list = obs_list_cube
                else:
                    obs_list = [
                        obs for obs in obs_list_cube if obs in obs_list
                    ]

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
                if '{}_{}_{}_HEALPix.fits'.format(obs, weighting, cube) not in data_files:
                    print 'WARNING: File {}/output_data/{}_{}_{}_HEALPix.fits not found. Excluding observation {} from the average.'.format(
                        fhd_run_path, obs, weighting, cube, obs
                    )
                    exclude_obs_list.append(obs)
                    continue
        obs_weights = np.array([
            obs_weights[ind] for ind in range(len(obs_weights))
            if obs_list[ind] not in exclude_obs_list
        ])
        obs_list = [obs for obs in obs_list if obs not in exclude_obs_list]

        if not quiet:
            print 'Averaging {} observations from {}'.format(len(obs_list), fhd_run_path)

        for obs_ind, obsid in enumerate(obs_list):
            if not quiet:
                print 'Loading observation {} of {}'.format(obs_ind+1, len(obs_list))
            maps = []
            for cube_ind, cube in enumerate(cube_names):
                map = load_map('{}/output_data/{}_{}_{}_HEALPix.fits'.format(
                    fhd_run_path, obsid, weighting, cube
                ), quiet=quiet)
                if obs_ind == 0 and cube_ind == 0 and path_ind == 0:  # Use first map for conventions
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

            if apply_rm_correction:
                if use_rms is None:
                    rm_val = None
                else:
                    rm_val = use_rms[obs_ind]
                maps = rm_correction(
                    obsid, maps, use_single_freq_calc=False, rm_file=rm_file,
                    use_rm=rm_val
                )

            if apply_radial_weighting:  # Restore obs structure if necessary
                obs_struct = scipy.io.readsav(
                    '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
                )['obs']
                obs_vec = hp.pixelfunc.ang2vec(
                    float(obs_struct['obsra']), float(obs_struct['obsdec']),
                    lonlat=True
                )

            for pix in maps[0].pix_arr:  # Use the first cube for the pixel list
                if apply_radial_weighting:
                    pix_vec = hp.pix2vec(nside, pix, nest=nest)
                    rad_weight = obs_radial_weighting_function(
                        hp.rotator.angdist(pix_vec, obs_vec)*180./np.pi
                    )
                    use_weight = rad_weight*obs_weights[obs_ind]
                else:
                    use_weight = obs_weights[obs_ind]
                if signal_array[0, pix] == hp.pixelfunc.UNSEEN:
                    for cube_ind in range(len(cube_names)):
                        signal_array[cube_ind, pix] = use_weight*maps[cube_ind].signal_arr[
                            np.where(maps[cube_ind].pix_arr == pix)
                        ]
                    weights_array[pix] = use_weight
                else:
                    for cube_ind in range(len(cube_names)):
                        signal_array[cube_ind, pix] += use_weight*maps[cube_ind].signal_arr[
                            np.where(maps[cube_ind].pix_arr == pix)
                        ]
                    weights_array[pix] += use_weight

    weighted_ave_signal_array = signal_array/weights_array[None, :]
    weighted_ave_signal_array[
        :, np.where(weights_array == hp.pixelfunc.UNSEEN)
    ] = hp.pixelfunc.UNSEEN
    weighted_ave_signal_array[
        :, np.where(weights_array == 0)
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


def rm_correction(
    obsid, maps, rm_file='/Users/rubybyrne/diffuse_survey_rm_tot.csv',
    start_freq_mhz=167., end_freq_mhz=198., use_single_freq_calc=False,
    use_rm=None
):

    if rm_file is None:
        rm_file='/Users/rubybyrne/diffuse_survey_rm_tot.csv'

    if len(maps) < 3:
        print 'ERROR: RM correction requires Stokes Q and U maps.'
        sys.exit(1)
    if use_rm is None: # Look up RM if none is provided explicitly
        rm_data = np.genfromtxt(
            rm_file, delimiter=',', dtype=None, names=True, encoding=None
        )
        if int(obsid) not in rm_data['ObsID']:
            print 'ERROR: Obsid {} not found in {}'.format(obsid, rm_file)
            sys.exit(1)
        rm = rm_data['RM'][np.where(rm_data['ObsID'] == int(obsid))][0]
    else:
        rm = use_rm

    c = 3.e8
    if start_freq_mhz == end_freq_mhz:
        use_single_freq_calc = True

    if use_single_freq_calc:
        reference_freq_mhz = np.mean([start_freq_mhz, end_freq_mhz])
        wavelength = c/(reference_freq_mhz*1.e6)
        rot_angle = rm * wavelength**2.
        new_q = np.cos(2*rot_angle)*maps[1].signal_arr + np.sin(2*rot_angle)*maps[2].signal_arr
        new_u = -np.sin(2*rot_angle)*maps[1].signal_arr + np.cos(2*rot_angle)*maps[2].signal_arr
    else:
        wl_max = c/(start_freq_mhz*1.e6)
        wl_min = c/(end_freq_mhz*1.e6)
        fresS_min, fresC_min = scipy.special.fresnel(2*np.sqrt(rm/np.pi+0j)*wl_min)
        fresS_max, fresC_max = scipy.special.fresnel(2*np.sqrt(rm/np.pi+0j)*wl_max)
        cos_int = (
            np.cos(2.*rm*wl_min**2.)/wl_min
            - np.cos(2.*rm*wl_max**2.)/wl_max
            + 2*np.sqrt(np.pi*rm+0j)*(fresS_min-fresS_max)
        )
        sin_int = (
            np.sin(2.*rm*wl_min**2.)/wl_min
            - np.sin(2.*rm*wl_max**2.)/wl_max
            - 2*np.sqrt(np.pi*rm+0j)*(fresC_min-fresC_max)
        )
        const = (1/wl_min - 1/wl_max)/(cos_int**2 + sin_int**2)
        new_q = np.real(const*(
            cos_int*maps[1].signal_arr + sin_int*maps[2].signal_arr
        ))
        new_u = np.real(const*(
            -sin_int*maps[1].signal_arr + cos_int*maps[2].signal_arr
        ))

    maps[1].signal_arr = new_q
    maps[2].signal_arr = new_u
    return maps


def calculate_variance_healpix_maps(
    fhd_run_paths, obs_lists=None, obs_weights_list=None,
    saved_averaged_maps=None, nside=None, cube_names=['Residual_I'],
    weighting='uniform', apply_radial_weighting=False, apply_rm_correction=False
):

    if isinstance(fhd_run_paths, str): #check if string
        fhd_run_paths = [fhd_run_paths]
    if obs_lists is not None:
        if len(fhd_run_paths) == 1 and isinstance(obs_lists[0], str):
            obs_lists = [obs_lists]
        elif len(fhd_run_paths) != len(obs_lists):
            print 'ERROR: obs_list not provided for each path'
            sys.exit(1)
    if obs_weights_list is not None:
        if len(fhd_run_paths) == 1 and isinstance(obs_lists[0], str):
            obs_weights_list = [obs_weights_list]
        elif len(fhd_run_paths) != len(obs_weights_list):
            print 'ERROR: obs_weights_list not provided for each path'
            sys.exit(1)

    if saved_averaged_maps is None:
        averaged_maps, null = average_healpix_maps(
            fhd_run_paths, obs_lists=obs_lists,
            obs_weights_list=obs_weights_list, nside=nside,
            cube_names=cube_names, weighting=weighting,
            apply_radial_weighting=apply_radial_weighting,
            apply_rm_correction=apply_rm_correction
        )
    else:
        if len(saved_averaged_maps) != len(cube_names):
            print 'WARNING: Invalid saved_averaged_maps paths provided. Recalculating data average.'
            averaged_maps, null = average_healpix_maps(
                fhd_run_paths, obs_lists=obs_lists,
                obs_weights_list=obs_weights_list, nside=nside,
                cube_names=cube_names, weighting=weighting,
                apply_radial_weighting=apply_radial_weighting,
                apply_rm_correction=apply_rm_correction
            )
        else:
            print 'Loading saved data averages.'
            averaged_maps = []
            for cube_ind in range(len(cube_names)):
                map = load_map(saved_averaged_maps[cube_ind])
                averaged_maps.append(map)
    for cube_ind in range(len(cube_names)):
        averaged_maps[cube_ind].explicit_to_implicit_ordering()

    for path_ind, fhd_run_path in enumerate(fhd_run_paths):

        if obs_lists is not None:
            obs_list = obs_lists[path_ind]
        else:
            obs_list = None
        if obs_weights_list is not None:
            obs_weights = obs_weights_list[path_ind]
        else:
            obs_weights = None

        if fhd_run_path[-1] == '/':
            fhd_run_path = fhd_run_path[:-1]

        if obs_list is None:  # use all obs in the data directory
            if obs_weights is not None:
                print 'WARNING: No obs_list provided. Disregarding obs_weights and using equal weighting.'
                obs_weights = None
            data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
            for cube_ind, cube in enumerate(cube_names):
                data_files_cube = [
                    file for file in data_files
                    if '_{}_{}_HEALPix.fits'.format(weighting, cube) in file
                ]
                obs_list_cube = [file[0:10] for file in data_files_cube]
                if cube_ind == 0:
                    obs_list = obs_list_cube
                else:
                    obs_list = [
                        obs for obs in obs_list_cube if obs in obs_list
                    ]

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
                if '{}_{}_{}_HEALPix.fits'.format(obs, weighting, cube) not in data_files:
                    print 'WARNING: File {}/output_data/{}_{}_{}_HEALPix.fits not found. Excluding observation {} from the average.'.format(
                        fhd_run_path, obs, weighting, cube, obs
                    )
                    exclude_obs_list.append(obs)
                    continue
        obs_weights = np.array([
            obs_weights[ind] for ind in range(len(obs_weights))
            if obs_list[ind] not in exclude_obs_list
        ])
        obs_list = [obs for obs in obs_list if obs not in exclude_obs_list]

        print 'Calculating variance for {} observations'.format(len(obs_list))

        for obs_ind, obsid in enumerate(obs_list):
            maps = []
            for cube_ind, cube in enumerate(cube_names):
                map = load_map('{}/output_data/{}_{}_{}_HEALPix.fits'.format(
                    fhd_run_path, obsid, weighting, cube
                ))
                if obs_ind == 0 and cube_ind == 0 and path_ind == 0:  # Use first map for conventions
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

            if apply_radial_weighting:  # Restore obs structure if necessary
                obs_struct = scipy.io.readsav(
                    '{}/metadata/{}_obs.sav'.format(fhd_run_path, obsid)
                )['obs']
                obs_vec = hp.pixelfunc.ang2vec(
                    float(obs_struct['obsra']), float(obs_struct['obsdec']),
                    lonlat=True
                )

            for pix in maps[0].pix_arr:  # Use the first cube for the pixel list
                if apply_radial_weighting:
                    pix_vec = hp.pix2vec(nside, pix, nest=nest)
                    rad_weight = obs_radial_weighting_function(
                        hp.rotator.angdist(pix_vec, obs_vec)*180./np.pi
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
    snr_maps = []
    for cube_ind in range(len(cube_names)):
        var_map = HealpixMap(
            np.squeeze(weighted_ave_var_array[cube_ind, :]), [], nside, nest=nest,
            coords=coords, quiet=True
        )
        var_map.implicit_to_explicit_ordering()
        variance_maps.append(var_map)
        averaged_maps[cube_ind].implicit_to_explicit_ordering()
        stddev_map = HealpixMap(
            np.squeeze(weighted_ave_var_array[cube_ind, :]), [], nside, nest=nest,
            coords=coords, quiet=True
        )
        stddev_map.implicit_to_explicit_ordering()
        stddev_map.signal_arr = np.sqrt(stddev_map.signal_arr)
        snr = divide_healpix_maps(averaged_maps[cube_ind], stddev_map)
        snr.signal_arr = np.abs(snr.signal_arr)
        snr_maps.append(snr)

    weights_map = HealpixMap(
        weights_array, [], nside, nest=nest, coords=coords, quiet=True
    )
    weights_map.implicit_to_explicit_ordering()

    nsamples_map = HealpixMap(
        nsamples_array, [], nside, nest=nest, coords=coords, quiet=True
    )
    nsamples_map.implicit_to_explicit_ordering()

    return averaged_maps, variance_maps, snr_maps, weights_map, nsamples_map


def obs_radial_weighting_function(
    dist, max_dist=10., taper_width=6.,
    mask_edges=False # if true, this will set some weights to zero
):

    if mask_edges:
        zero_point = 0.
    else:
        zero_point = 0.00001 # Don't zero things out completely
    if dist < max_dist-taper_width:
        weight = 1.
    elif dist > max_dist:
        weight = zero_point
    else:
        weight = (0.5-zero_point)*np.cos(
            np.pi*(dist+taper_width-max_dist)/taper_width
        )+0.5+zero_point
    return weight


def combine_maps_nearest_data(
    fhd_run_path, obs_list=None, obs_list_file=None, nside=None,
    cube_names=['Residual_I']
):
    # This function loads a set of maps from an FHD output directory (or a list
    # of observations) and returns a combined map.
    # Each pixel of the combined map is based on the closest observation to
    # that pixel.

    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    if obs_list is None:
        if obs_list_file is None:  # use all obs in the data directory
            data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
            for cube_ind, cube in enumerate(cube_names):
                data_files_cube = [
                    file for file in data_files
                    if '_uniform_{}_HEALPix.fits'.format(cube) in file
                ]
                obs_list_cube = [file[0:10] for file in data_files_cube]
                if cube_ind == 0:
                    obs_list = obs_list_cube
                else:
                    obs_list = [
                        obs for obs in obs_list_cube if obs in obs_list
                    ]
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


def write_data_to_standard_fits(maps, save_filename, history_str=''):

    if len(maps) != 4:
        print 'ERROR: Must provide Stokes I,Q,U,V maps.'
        sys.exit()

    data = np.zeros((len(maps[0].signal_arr), 1, 4))
    data[:, 0, 0] = maps[0].signal_arr

    for map_ind in range(1, 4):
        if maps[map_ind].nside != maps[0].nside:
            print 'ERROR: Map nsides do not match.'
            sys.exit()
        if maps[map_ind].nest != maps[0].nest:
            print 'ERROR: Map nest conventions do not match.'
            sys.exit()
        if maps[map_ind].coords != maps[0].coords:
            print 'ERROR: Map coordinate system conventions do not match.'
            sys.exit()
        if len(maps[map_ind].pix_arr) != len(maps[0].pix_arr):
            print 'ERROR: Maps have different numbers of pixels.'
            sys.exit()
        if not np.array_equal(maps[map_ind].pix_arr, maps[0].pix_arr):
            print 'WARNING: Pixel ordering does not match. Reordering.'
            new_signal_arr = np.zeros((len(maps[map_ind].pix_arr)))
            for pix_ind in range(len(maps[map_ind].pix_arr)):
                new_signal_arr[pix_ind] = maps[map_ind].signal_arr[
                    np.where(
                        maps[map_ind].pix_arr == maps[0].pix_arr[pix_ind]
                    )[0][0]
                ]
            maps[map_ind].signal_arr = new_signal_arr
            maps[map_ind].pix_arr = maps[0].pix_arr
        data[:, 0, map_ind] = maps[map_ind].signal_arr

    header = fits.Header()

    # Conforming to fits format
    header['SIMPLE'] = True
    header['BITPIX'] = 32
    header['NAXIS'] = 3
    header['NSIDE'] = maps[0].nside
    if maps[0].nest:
        ordering = 'nested'
    else:
        ordering = 'ring'
    header['ORDERING'] = ordering
    header['BUNIT'] = 'Jy/sr'
    header['COORDSYS'] = maps[0].coords_healpy_conv
    header['AUTHOR'] = 'Ruby Byrne'
    header['INSTRUME'] = 'MWA'
    header['HISTORY'] = history_str

    ax_nums = {'pixel': 1, 'freq': 2, 'pol': 3}

    # set up pixel axis
    header['CTYPE' + str(ax_nums['pixel'])] = \
        ('Pix_Ind', 'Index into pixel array in HPX_INDS extension.')
    header['CRVAL' + str(ax_nums['pixel'])] = 1
    header['CRPIX' + str(ax_nums['pixel'])] = 1
    header['CDELT' + str(ax_nums['pixel'])] = 1

    # set up frequency axis
    header['CTYPE' + str(ax_nums['freq'])] = 'FREQ'
    header['CUNIT' + str(ax_nums['freq'])] = 'MHz'
    header['CRVAL' + str(ax_nums['freq'])] = 182.
    header['CRPIX' + str(ax_nums['freq'])] = 1
    header['CDELT' + str(ax_nums['freq'])] = 30.

    # set up polarization axis
    header['CTYPE' + str(ax_nums['pol'])] = \
        ('STOKES', 'Polarization integers: I,Q,U,V=1,2,3,4')
    header['CRVAL' + str(ax_nums['pol'])] = 1
    header['CRPIX' + str(ax_nums['pol'])] = 1
    header['CDELT' + str(ax_nums['pol'])] = 1

    primary_hdu = fits.PrimaryHDU(data=data, header=header)
    hdulist = fits.HDUList([primary_hdu])

    # set up pixel array
    c1 = fits.Column(name='HPX_INDS', format='K', array=maps[0].pix_arr)
    coldefs = fits.ColDefs([c1])
    hpx_hdu = fits.BinTableHDU.from_columns(coldefs)
    hpx_hdu.header['EXTNAME'] = 'HPX_INDS'
    hdulist.append(hpx_hdu)

    hdulist.writeto(save_filename)
