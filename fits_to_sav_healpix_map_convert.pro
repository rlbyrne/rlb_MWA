pro fits_to_sav_healpix_map_convert

    stokes_maps_paths = ['/Users/ruby/EoR/sky_maps/StokesI_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesQ_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesU_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesV_nearest_short_baselines.fits']
    outpath = '/Users/ruby/EoR/sky_maps/nearest_short_baselines_Aug2019.sav'

    model_arr = ptrarr(4)
    for ind=0,3 do begin
        read_fits_map, stokes_maps_paths[0], data_new, nside=nside_new
        if ind ne 0 then begin
            if ~array_equal(data_new[*,1], data[*,1]) or nside_new ne nside then begin
                print, 'ERROR: index mismatch'
                return
            endif
        endif
        data = data_new
        nside = nside_new
        model_arr[ind] = ptr_new(data[*,0])
    endfor
    hpx_inds = data[*,1]
    save, model_arr, hpx_inds, nside, filename=outpath

end