pro fits_to_sav_healpix_map_convert, from_jy_per_beam=from_jy_per_beam, obsfile=obsfile

    stokes_maps_paths = ['/Users/ruby/EoR/sky_maps/StokesI_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesQ_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesU_nearest_short_baselines.fits',$
        '/Users/ruby/EoR/sky_maps/StokesV_nearest_short_baselines.fits']
    outpath = '/Users/ruby/EoR/sky_maps/nearest_short_baselines_altnorm_Aug2019.sav'
    
    if keyword_set(from_jy_per_beam) and ~keyword_set(obsfile) then begin
      print, 'ERROR: No obs file provided'
      return
    endif

    if keyword_set(from_jy_per_beam) then begin
      obs = getvar_savefile(obsfile, 'obs')
      beam_area = beam_width_calculate(obs, /area_return, /radians) ;Synthesized beam area in sr
    endif
    
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
        if keyword_set(from_jy_per_beam) then begin
            model_arr[ind] = ptr_new(data[*,0]/beam_area*pix_area) ;Convert from Jy/beam to Jy/pixel
        endif else begin
            model_arr[ind] = ptr_new(data[*,0]*pix_area) ;Convert from Jy/sr to Jy/pixel
        endelse
    endfor
    hpx_inds = long(data[*,1])
    save, model_arr, hpx_inds, nside, filename=outpath

end