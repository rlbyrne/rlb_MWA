pro fits_to_sav_healpix_map_convert, units_in=units_in, units_out=units_out, obsfile=obsfile

    from_jy_per_beam = 1
    if ~keyword_set(units_in) then units_in='Jy/sr'
    if ~keyword_set(units_out) then units_out='Jy/sr'
    
    obsfile = '/Users/ruby/Downloads/1130773144_obs.sav'
    ;stokes_maps_paths = ['/Users/ruby/EoR/sky_maps/StokesI_nearest_short_baselines.fits',$
    ;  '/Users/ruby/EoR/sky_maps/StokesQ_nearest_short_baselines.fits',$
    ;  '/Users/ruby/EoR/sky_maps/StokesU_nearest_short_baselines.fits',$
    ;  '/Users/ruby/EoR/sky_maps/StokesV_nearest_short_baselines.fits']
    ;outpath = '/Users/ruby/EoR/sky_maps/nearest_Jan2020norm_Aug2019.sav'
    stokes_maps_paths = ['/Users/ruby/EoR/diffuse_normalization/fhd_rlb_diffuse_model_limit_dimension_Jan2019/1061316296_weighted_Model_I_HEALPix.fits',$
      '/Users/ruby/EoR/diffuse_normalization/fhd_rlb_diffuse_model_limit_dimension_Jan2019/1061316296_weighted_Model_Q_HEALPix.fits',$
      '/Users/ruby/EoR/diffuse_normalization/fhd_rlb_diffuse_model_limit_dimension_Jan2019/1061316296_weighted_Model_U_HEALPix.fits',$
      '/Users/ruby/EoR/diffuse_normalization/fhd_rlb_diffuse_model_limit_dimension_Jan2019/1061316296_weighted_Model_V_HEALPix.fits']
    outpath = '/Users/ruby/EoR/diffuse_normalization/fhd_rlb_diffuse_model_limit_dimension_Jan2019/single_obs_stokes_maps.sav'

    case units_in of
      'Jy/cart_pixel':
      'Jy/hpx_pixel':
      'Jy/sr':
      'Jy/beam':
      else: begin
          print, 'ERROR: Invalid units_in parameter'
          return
        end
    endcase
    
    case units_out of
      'Jy/sr':
      'Jy/hpx_pixel':
      else: begin
          print, 'ERROR: Invalid units_out parameter'
          return
        end
    endcase
    
    if units_in eq 'Jy/beam' then begin
      if ~keyword_set(obsfile) then begin
        print, 'ERROR: No obs file provided'
        return
      endif
      obs = getvar_savefile(obsfile, 'obs')
      beam_area = beam_width_calculate(obs, /area_return, /radians) ;Synthesized beam area in sr
    endif
    
    if units_in eq 'Jy/cart_pixel' then begin
      if ~keyword_set(obsfile) then begin
        print, 'ERROR: No obs file provided'
        return
      endif
      obs = getvar_savefile(obsfile, 'obs')
      cart_pix_area = (obs.degpix*!DtoR)^2.
    endif
    
    model_arr = ptrarr(4)
    for ind=0,3 do begin
        read_fits_map, stokes_maps_paths[ind], data_new, nside=nside_new, ordering=ordering
        
        ;for FHD outputs, the first column is the pixels and the second column is the data
        ;third and fourth columns report N_obs and serror
        ;turn this off if the output is not directly out of FHD
        data_reordered = make_array(n_elements(data_new[*,0]), 2, /float)
        data_reordered[*,0] = data_new[*,1]
        data_reordered[*,1] = data_new[*,0]
        data_new = data_reordered
        
        if ordering eq 'NESTED' then begin
            print, 'Reordering Healpix map: nested to ring'
            ;Need to convert to implicit ordering
            data_implicit = make_array(12*nside_new^2, /float, value=-1.6375e+30)
            for pix=0,n_elements(data_new[*,0])-1 do begin
                data_implicit[long(data_new[pix,1])] = data_new[pix,0]
            endfor
            data_implicit = reorder(data_implicit, /n2r)
            keep_pixels = where(data_implicit ne -1.6375e+30)
            data_new[*,0] = data_implicit[keep_pixels]
            data_new[*,1] = keep_pixels
        endif
        if ind ne 0 then begin
          if ~array_equal(data_new[*,1], data[*,1]) or nside_new ne nside then begin
            print, 'ERROR: index mismatch'
            return
          endif
        endif
        data = data_new
        nside = nside_new
        pix_area = 4*!pi/(12*nside^2.) ;Pixel area in sr
        
        case units_in of
          'Jy/cart_pixel': model_arr[ind] = ptr_new(data[*,0]/cart_pix_area) ;Convert from Jy/cart_pixel to Jy/sr
          'Jy/hpx_pixel': model_arr[ind] = ptr_new(data[*,0]/pix_area) ;Convert from Jy/hpx_pixel to Jy/sr
          'Jy/sr': model_arr[ind] = ptr_new(data[*,0]) ;Keep in Jy/sr
          'Jy/beam': model_arr[ind] = ptr_new(data[*,0]*beam_area) ;Convert from Jy/beam to Jy/sr
        endcase
        
        if units_out eq 'Jy/hpx_pixel' then *model_arr[ind]*=pix_area

    endfor
    hpx_inds = long(data[*,1])
    print, 'Saving map to '+outpath
    save, model_arr, hpx_inds, nside, filename=outpath

end