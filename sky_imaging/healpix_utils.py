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
                 coords='equitorial'):
        if pix_arr == []:
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
        self.signal_arr = signal_arr
        self.pix_arr = pix_arr
        self.nside = nside
        self.nest = nest
        self.coords = coords

    def get_ra_dec(self, ra_cut=270):
        ra_arr = []
        dec_arr = []
        rot = hp.rotator.Rotator(coord=['G', 'C'])
        for pixel in self.pix_arr:
            if self.coords == 'galactic':
                theta_gal, phi_gal = hp.pixelfunc.pix2ang(
                    self.nside, pixel, nest=self.nest
                    )
                theta_eq, phi_eq = rot(theta_gal, phi_gal)
                ra = phi_eq*180/math.pi
                dec = 90.-theta_eq*180/math.pi
            elif self.coords == 'equitorial':
                ra, dec = hp.pixelfunc.pix2ang(
                    self.nside, pixel, nest=self.nest, lonlat=True
                    )
            else:
                print 'ERROR: Coordinates must be galactic or equitorial.'
                sys.exit(1)
            if ra > ra_cut:
                ra -= 360.
            elif ra < ra_cut-360.:
                ra += 360.
            ra_arr.append(ra)
            dec_arr.append(dec)
        self.ra_arr = ra_arr
        self.dec_arr = dec_arr

    def get_pixel_corners(self, ra_cut=270):
        pix_corner_ras_arr = []
        pix_corner_decs_arr = []
        rot = hp.rotator.Rotator(coord=['G', 'C'])
        for pixel in self.pix_arr:
            corner_coords = hp.boundaries(self.nside, pixel, step=1,
                                          nest=self.nest)
            if self.coords == 'galactic':
                thetas_gal, phis_gal = hp.pixelfunc.vec2ang(
                    np.transpose(corner_coords), lonlat=False)
                ras = [0]*len(thetas_gal)
                decs = [0]*len(thetas_gal)
                for corner in range(len(thetas_gal)):
                    theta_eq, phi_eq = rot(thetas_gal[corner],
                                           phis_gal[corner])
                    ras[corner] = phi_eq*180/math.pi
                    decs[corner] = 90. - theta_eq*180/math.pi
            elif self.coords == 'equitorial':
                ras, decs = hp.pixelfunc.vec2ang(np.transpose(corner_coords),
                                                 lonlat=True)
            else:
                print 'ERROR: Coordinates must be galactic or equitorial.'
                sys.exit(1)
            # Correct for the branch cut
            for ind, ra_val in enumerate(ras):
                if ra_val > ra_cut:
                    ras[ind] = ra_val-360.
                if ra_val < ra_cut-360.:
                    ras[ind] = ra_val+360.
                # Cluster pixel corners that are separated by the branch cut
                if ind > 0:
                    if abs(ras[ind]-ras[0]-360.) < abs(ras[ind]-ras[0]):
                        ras[ind] -= 360.
                    elif abs(ras[ind]-ras[0]+360.) < abs(ras[ind]-ras[0]):
                        ras[ind] += 360.
            pix_corner_ras_arr.append(ras)
            pix_corner_decs_arr.append(decs)
        self.pix_corner_ras_arr = pix_corner_ras_arr
        self.pix_corner_decs_arr = pix_corner_decs_arr

    def explicit_to_implicit_ordering(self):
        signal_vals_implicit = [hp.pixelfunc.UNSEEN]*12*self.nside**2
        for i in range(len(self.signal_arr)):
            signal_vals_implicit[self.pix_arr[i]] = self.signal_arr[i]
        self.signal_arr = signal_vals_implicit
        self.pix_arr = []

    def implicit_to_explicit_ordering(self):
        use_inds = [ind for ind in range(len(self.signal_arr))
                    if self.signal_arr[ind] != hp.pixelfunc.UNSEEN]
        if self.pix_arr == []:
            self.pix_arr = use_inds
        else:
            self.pix_arr = [self.pix_arr[ind] for ind in use_inds]
        self.signal_arr = [self.signal_arr[ind] for ind in use_inds]

    def reorder_nest_to_ring(self):
        if self.nest:
            self.explicit_to_implicit_ordering()
            self.signal_arr = hp.pixelfunc.reorder(filtered_map, n2r=True)
            self.nest = False
            self.implicit_to_explicit_ordering()

    def reorder_ring_to_nest(self):
        if not self.nest:
            self.explicit_to_implicit_ordering()
            self.signal_arr = hp.pixelfunc.reorder(filtered_map, r2n=True)
            self.nest = True
            self.implicit_to_explicit_ordering()

    def resample(self, nside):
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
        alm = hp.sphtfunc.map2alm(signal_vals_implicit, lmax=lmax)
        alm = np.copy(alm)
        lmax = hp.sphtfunc.Alm.getlmax(len(alm))
        l, m = hp.sphtfunc.Alm.getlm(lmax)
        lm = np.zeros([len(alm), 2])
        lm[:, 0] = l
        lm[:, 1] = m
        self.alm = spherical_harmonics_amp
        self.lm = spherical_harmonics_lm

    def filter_map(self, lmin=None, lmax=None, filter_width=0):
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
        filtered_map = hp.sphtfunc.alm2map(alm, self.nside)
        self.signal_arr = filtered_map
        self.pix_arr = []
        # If input was nested ordering, convert to nested
        if self.nest:
            self.nest = False
            self.reorder_ring_to_nest()
            filtered_map = hp.pixelfunc.reorder(filtered_map, r2n=True)
        else:
            self.implicit_to_explicit_ordering()

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


def average_healpix_maps(maps_arr):

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


def combine_maps_nearest_data(
    fhd_run_path, obs_list_file=None, nside=None, cube_name='Residual_I'
):

    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    if obs_list_file is None:  # use all obs in the data directory
        data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
        data_files = [
            file for file in data_files
            if '_uniform_{}_HEALPix.fits'.format(cube_name) in file
        ]
        obs_list = [file[0:10] for file in data_files]
    else:  # use the obs file list
        obs_list = open(obs_list_file, 'r').readlines()
        # strip newline characters
        obs_list = [obs.strip() for obs in obs_list]
        # remove duplicates
        obs_list = list(set(obs_list))

    print 'Combining {} observations'.format(len(obs_list))

    healpix_maps = []
    all_pixels = []
    for obsid in obs_list:
        map = load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
            fhd_run_path, obsid, cube_name
        ))
        if nside is not None:
            if map.nside != nside:
                map.resample(nside)
        healpix_maps.append(map)
        nside = map.nside  # use the nside of the first obs
        all_pixels.extend(map.pix_arr)

    all_pixels = list(set(all_pixels))  # remove duplicates
    if len(set([map.nest for map in healpix_maps])) > 1:
        print 'ERROR: Nest conventions do not match. Exiting.'
        sys.exit(1)
    nest = healpix_maps[0].nest
    coords = healpix_maps[0].coords

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

    signal_arr = []
    pix_arr = []
    for pix in all_pixels:
        vec = hp.pix2vec(nside, pix, nest=nest)
        distances = []
        for map_ind in range(len(healpix_maps)):
            dist = (
                (vec[0]-obs_centers[map_ind][0])**2.
                + (vec[1]-obs_centers[map_ind][1])**2.
                + (vec[2]-obs_centers[map_ind][2])**2.
            )
            distances.append(dist)
        map_indices = np.argsort(distances)
        use_map = 0
        while pix not in healpix_maps[map_indices[use_map]].pix_arr:
            use_map += 1
        signal_arr.append(healpix_maps[map_indices[use_map]].signal_arr[
            (healpix_maps[map_indices[use_map]].pix_arr).index(pix)
        ])
        pix_arr.append(pix)

    combined_map = HealpixMap(
        signal_arr, pix_arr, nside, nest=nest, coords=coords
    )

    return combined_map


def combine_maps_nearest_data_memory_efficient(
    fhd_run_path, obs_list_file=None, nside=None, cube_names=['Residual_I']
):

    if fhd_run_path[-1] == '/':
        fhd_run_path = fhd_run_path[:-1]

    if obs_list_file is None:  # use all obs in the data directory
        data_files = os.listdir('{}/output_data/'.format(fhd_run_path))
        data_files = [
            file for file in data_files
            if '_uniform_{}_HEALPix.fits'.format(cube_name) in file
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
        obs_maps = []
        for cube in cube_names:
            map = load_map('{}/output_data/{}_uniform_{}_HEALPix.fits'.format(
                fhd_run_path, obsid, cube
            ))
            if nside is not None:
                if map.nside != nside:
                    map.resample(nside)
            obs_maps.append(map)
        nside = obs_maps[0].nside  # use nside of the first map of first obs
        nest = obs_maps[0].nest
        coords = obs_maps[0].coords

        if obs_ind == 0:
            pix_array = obs_maps[0].pix_arr
            signal_array = [[]*len(cube_names)]
            signal_array[0] = obs_maps[0].signal_arr
            for cube_ind in range(1, len(cube_names)):
                signal_array[cube_ind].extend([
                    obs_maps[cube_ind].signal_arr[
                        (obs_maps[cube_ind].pix_arr).index(pix)
                    ] for pix in obs_maps[0].pix_arr
                ])
        else:
            add_pixels = [
                pix for pix in obs_maps[0].pix_arr if pix not in pix_array
            ]
            calculate_pixels = [
                pix for pix in obs_maps[0].pix_arr if pix in pix_array
            ]
            # Add pixels that aren't already represented in the map
            pix_arr.extend(add_pixels)
            for cube_ind in range(len(cube_names)):
                signal_array[cube_ind].extend([
                    obs_maps[cube_ind].signal_arr[
                        (obs_maps[cube_ind].pix_arr).index(pix)
                    ] for pix in add_pixels
                ])
            # Calculate obs center distances for pixels represented in the map
            for pix in calculate_pixels:
                vec = hp.pix2vec(nside, pix, nest=nest)
                distances = [
                    (vec[0]-obs_centers[obs_ind][0])**2.
                    + (vec[1]-obs_centers[obs_ind][1])**2.
                    + (vec[2]-obs_centers[obs_ind][2])**2.
                    for obs_ind in range(obs_ind+1)
                ]
                if distances.index(min(distances)) == obs_ind:
                    pixel_index = pix_array.index(pix)
                    for cube_ind in range(len(cube_names)):
                        signal_array[cube_ind][
                            pixel_index
                        ] = obs_maps[cube_ind].signal_arr[
                            (obs_maps[cube_ind].pix_arr).index(pix)
                        ]

    combined_maps = []
    for cube_ind in range(len(cube_names)):
        map = HealpixMap(
            signal_array[cube_ind], pix_array, nside, nest=nest, coords=coords
        )
        combined_maps.append(map)

    return combined_maps
