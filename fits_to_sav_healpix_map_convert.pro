pro fits_to_sav_healpix_map_convert, units_in=units_in, units_out=units_out, obsfile=obsfile, $
  use_fhd_output_convention=use_fhd_output_convention, use_stokes_v=use_stokes_v

    if ~keyword_set(units_in) then units_in='Jy/sr'
    if ~keyword_set(units_out) then units_out='Jy/sr'
    if n_elements(use_fhd_output_convention) eq 0 then use_fhd_output_convention=0
    if n_elements(use_stokes_v) eq 0 then use_stokes_v=1

    ;path = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020'
    ;stokes_maps_paths = [path+'/StokesI_average_map_empirical_rm_in_eor0.fits', path+'/StokesQ_average_map_1131454296_rm_undone.fits', $
    ;  path+'/StokesU_average_map_1131454296_rm_undone.fits', path+'/StokesV_average_map_empirical_rm_in_eor0.fits']
    ;outpath = path+'/average_map_I.sav'
    
    csv_data = read_csv('/Users/rubybyrne/mwa_2013_rm.csv', header=csv_header)
    obsids = csv_data.field1
    for obs_ind = 0, n_elements(obsids)-1 do begin
    
      path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/fhd_rlb_model_diffuse_Aug2020_2/output_data'
      stokes_maps_paths = ['/Users/rubybyrne/diffuse_survey_plotting_Aug2020/StokesI_average_map_empirical_rm_in_eor0.fits', $
        '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/rm_corrected/fits_files/'+strtrim(obsids[obs_ind],2)+'_RMcorrected_StokesQ.fits', $
        '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/rm_corrected/fits_files/'+strtrim(obsids[obs_ind],2)+'_RMcorrected_StokesU.fits']
      outpath = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/rm_corrected/sav_files/diffuse_rm_corrected_'+strtrim(obsids[obs_ind],2)+'.sav'
  
  
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
      if use_stokes_v then read_n_maps=4 else read_n_maps=3
      for ind=0,read_n_maps-1 do begin
          read_fits_map, stokes_maps_paths[ind], data_new, nside=nside_new, ordering=ordering
          
          ;for FHD outputs, the first column is the pixels and the second column is the data
          ;third and fourth columns report N_obs and serror
          if use_fhd_output_convention then begin
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
      if ~use_stokes_v then model_arr[3] = ptr_new(make_array(n_elements(data[*,1]), /float, value=0.))
      hpx_inds = long(data[*,1])
      print, 'Saving map to '+outpath
      save, model_arr, hpx_inds, nside, filename=outpath
  
    endfor

end