pro fits_to_sav_healpix_map_convert, units_in=units_in, units_out=units_out, obsfile=obsfile

    from_jy_per_beam = 1
    if ~keyword_set(units_in) then units_in='Jy/sr'
    if ~keyword_set(units_out) then units_out='Jy/sr'
    
    ;path = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    ;stokes_maps_paths = [path+'/StokesI_average_map_no_rm_correction_more_obs.fits', path+'/StokesQ_average_map_no_rm_correction_more_obs.fits', $
    ;  path+'/StokesU_average_map_no_rm_correction_more_obs.fits', path+'/StokesV_average_map_no_rm_correction_more_obs.fits']
    ;outpath = path+'/averaged_diffuse_no_rm_correction_May2020.sav'
    
    ;stokes_maps_paths = ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/explore_diffuse_uv_planes/fhd_rlb_model_diffuse_May2020/output_data/1061316296_weighted_Model_I_HEALPix.fits',$
    ;    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/explore_diffuse_uv_planes/fhd_rlb_model_diffuse_May2020/output_data/1061316296_weighted_Model_Q_HEALPix.fits']
    ;outpath = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/explore_diffuse_uv_planes/fhd_rlb_model_diffuse_May2020/diffuse_roundtrip_test_May2020.sav'
      
    ;path = '/Users/rubybyrne/diffuse_survey_plotting_May2020'
    ;stokes_maps_paths = [path+'/StokesI_average_map_more_obs.fits', path+'/StokesQ_average_map_more_obs.fits', $
    ;  path+'/StokesU_average_map_more_obs.fits', path+'/StokesV_average_map_more_obs.fits']
    ;outpath = path+'/averaged_diffuse_May2020.sav'

    path = '/Users/rubybyrne/diffuse_survey_plotting_May2020/eor0_plots'
    for map_ind =1,8 do begin
      map_ind_str = strtrim(string(map_ind),2)
      stokes_maps_paths = [path+'/StokesI_eor0_map'+map_ind_str+'_tapered.fits', path+'/StokesQ_eor0_map'+map_ind_str+'_tapered.fits', $
        path+'/StokesU_eor0_map'+map_ind_str+'_tapered.fits', path+'/StokesV_eor0_map'+map_ind_str+'_tapered.fits']
      outpath = path+'/eor0_map'+map_ind_str+'_May2020.sav'
  
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
      
      model_arr = ptrarr(n_elements(stokes_maps_paths))
      for ind=0,n_elements(stokes_maps_paths)-1 do begin
          read_fits_map, stokes_maps_paths[ind], data_new, nside=nside_new, ordering=ordering
          
          ;for FHD outputs, the first column is the pixels and the second column is the data
          ;third and fourth columns report N_obs and serror
          ;turn this off if the output is not directly out of FHD
          if 0 then begin
            data_reordered = make_array(n_elements(data_new[*,0]), 2, /float)
            data_reordered[*,0] = data_new[*,1]
            data_reordered[*,1] = data_new[*,0]
            data_new = data_reordered
          endif
          
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
      ;print, 'Saving map to '+outpath
      ;save, model_arr, hpx_inds, nside, filename=outpath
      
      ;save polarizations individually
      model_arr_filled = model_arr  
      pol_names = ['I','Q','U','V']
      for pol_ind=0,3 do begin
        model_arr = ptrarr(4)
        for ind=0,3 do model_arr[ind] = ptr_new(fltarr(n_elements(hpx_inds)))
        *model_arr[pol_ind] = *model_arr_filled[pol_ind]
        save, model_arr, hpx_inds, nside, filename=path+'/eor0_map'+map_ind_str+'_Stokes'+pol_names[pol_ind]+'_May2020.sav'
      endfor
    
    endfor

end