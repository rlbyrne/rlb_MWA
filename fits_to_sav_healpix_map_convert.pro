pro fits_to_sav_healpix_map_convert

    stokes_maps_paths = ['/Users/ruby/EoR/sky_maps/StokesI_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesQ_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesU_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesV_nearest_short_baselines.fits']
    obsfile = '/Users/ruby/Downloads/1130773144_obs.sav' ;Need an obs file for conversion from Jy/beam to Jy/pixel
    outpath = '/Users/ruby/EoR/sky_maps/nearest_short_baselines_Aug2019.sav'

    obs = getvar_savefile(obsfile, 'obs')
    beam_area = beam_width_calculate(obs, /area_return, /radians) ;Synthesized beam area in sr
    model_arr = ptrarr(4)
    for ind=0,3 do begin
        read_fits_map, stokes_maps_paths[ind], data_new, nside=nside_new
        if ind ne 0 then begin
            if ~array_equal(data_new[*,1], data[*,1]) or nside_new ne nside then begin
                print, 'ERROR: index mismatch'
                return
            endif
        endif
        data = data_new
        nside = nside_new
        pix_area = 4*!pi/(12*nside^2) ;Pixel area in sr
        model_arr[ind] = ptr_new(data[*,0]/beam_area*pix_area) ;Convert values to Jy/pixel
    endfor
    hpx_inds = long(data[*,1])
    save, model_arr, hpx_inds, nside, filename=outpath

end