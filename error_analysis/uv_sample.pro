Pro uv_sample, uv_img_clip=uv_img_clip, start_sample_index_x=start_sample_index_x, $
    start_sample_index_y=start_sample_index_y, diff_cubes=diff_cubes
    
  if uv_img_clip mod 1 ne 0 then begin
    uv_img_clip = round(uv_img_clip)
    print, 'ERROR: uv_img_clip must be an integer. Proceeding with uv_img_clip = ' + uv_img_clip
  endif
  uv_img_clip = number_formatter(uv_img_clip)
  ;sample_factors = [.0002,.0005,.001,.005,.01,.05,.1,.5,1,5]
  sample_factors = [.0002,.0005,.001,.005,.01,.05,.5,1,5]
  sample_factors = number_formatter(sample_factors)
  choose_terms = [1]
  polarizations = ['xx']
  if keyword_set(diff_cubes) then obsids = [''] else obsids = ['1061316176','1061316296']
  
  for obs = 0, n_elements(obsids)-1 do begin
    for sample = 0, n_elements(sample_factors)-1 do begin
      for term = 0, n_elements(choose_terms)-1 do begin
        for pol = 0, n_elements(polarizations)-1 do begin
        
          ;format sample factor for filename by replacing periods with "p"s:
          pos = strpos(sample_factors[sample],'.')
          if pos ne -1 then filename_sample = strmid(sample_factors[sample],0,pos) + 'p' + $
            strmid(sample_factors[sample], pos+1, strlen(sample_factors[sample])-pos) else $
            filename_sample = sample_factors[sample]
            
          ;check that the source file exists:
          if keyword_set(diff_cubes) then begin
            file_loc = '/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/'
            filename = 'cross_diff_noise_sim_flatUV'+filename_sample $
              +'_'+polarizations[pol]+'_term'+number_formatter(choose_terms[term])+'.sav'
          endif else begin
            file_loc = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsids[obs] + '_' $
              + sample_factors[sample] + '/ps/'
            filename = obsids[obs] + '_gridded_uvf__even_odd_joint_model_' + polarizations[pol] + $
              '_bh_kcube.idlsave'
          endelse
          test = file_test(file_loc+filename)
          if test ne 1 then begin
            print, 'ERROR: file '+file_loc + filename+' not found'
            continue
          endif
          
          ;restore cubes:
          print, 'RESTORING: ' + file_loc+filename
          restore, file_loc+filename
          
          ;sample and save cubes:
          sample_array_x = indgen(n_elements(kx_mpc))
          sample_array_y = indgen(n_elements(ky_mpc))
          if keyword_set(sample_start_index_x) then sample_array_x = sample_array_x + sample_start_index_x
          if keyword_set(sample_start_index_y) then sample_array_y = sample_array_y + sample_start_index_y
          keep_indices_x = where((sample_array_x mod uv_img_clip) eq 0)
          keep_indices_y = where((sample_array_y mod uv_img_clip) eq 0)
          kx_mpc = kx_mpc[keep_indices_x]
          ky_mpc = ky_mpc[keep_indices_y]
          
          save_loc = file_loc
          if keyword_set(diff_cubes) then begin
            save_filename = 'cross_diff_noise_sim_flatUV'+filename_sample $
              +'_uvsample'+uv_img_clip+'_'+polarizations[pol]+'_term'+number_formatter(choose_terms[term])+'.sav'
            ACube = (ACube[keep_indices_x,*,*])[*,keep_indices_y,*]
            ACube_sigma2 = (ACube_sigma2[keep_indices_x,*,*])[*,keep_indices_y,*]
            ACube_sim = (ACube_sim[keep_indices_x,*,*])[*,keep_indices_y,*]
            BCube = (BCube[keep_indices_x,*,*])[*,keep_indices_y,*]
            BCube_sigma2 = (BCube_sigma2[keep_indices_x,*,*])[*,keep_indices_y,*]
            BCube_sim = (BCube_sim[keep_indices_x,*,*])[*,keep_indices_y,*]
            cube_cross = (cube_cross[keep_indices_x,*,*])[*,keep_indices_y,*]
            sigma2 = (sigma2[keep_indices_x,*,*])[*,keep_indices_y,*]
            sim_cube_cross = (sim_cube_cross[keep_indices_x,*,*])[*,keep_indices_y,*]
            print, 'SAVING SAMPLED CUBES HERE: ' + save_loc + save_filename
            save, filename=save_loc+save_filename, kperp_lambda_conv, delay_params, kx_mpc, $
              ky_mpc, kz_mpc, ACube_sigma2, BCube_sigma2, $
              ACube, BCube, cube_cross, ACube_sim, BCube_sim, sim_cube_cross, sigma2
          endif else begin
            save_filename = obsids[obs] + '_gridded_uvf__even_odd_joint_uvsample'+uv_img_clip+'_model_' + $
              polarizations[pol] + '_bh_kcube.idlsave'
            data_diff_1 = (data_diff_1[keep_indices_x,*,*])[*,keep_indices_y,*]
            data_diff_2 = (data_diff_2[keep_indices_x,*,*])[*,keep_indices_y,*]
            data_sum_1 = (data_sum_1[keep_indices_x,*,*])[*,keep_indices_y,*]
            data_sum_2 = (data_sum_2[keep_indices_x,*,*])[*,keep_indices_y,*]
            n_freq_contrib = (n_freq_contrib[keep_indices_x,*])[*,keep_indices_y]
            sigma2_1 = (sigma2_1[keep_indices_x,*,*])[*,keep_indices_y,*]
            sigma2_2 = (sigma2_2[keep_indices_x,*,*])[*,keep_indices_y,*]
            sim_noise_diff_1 = (sim_noise_diff_1[keep_indices_x,*,*])[*,keep_indices_y,*]
            sim_noise_diff_2 = (sim_noise_diff_2[keep_indices_x,*,*])[*,keep_indices_y,*]
            wt_meas_ave = (wt_meas_ave[keep_indices_x,*])[*,keep_indices_y]
            wt_meas_min = (wt_meas_min[keep_indices_x,*])[*,keep_indices_y]
            save, filename=saveloc+save_filename, ave_power_freq, ave_power_uvf, ave_weights, $
              data_diff_1, data_diff_2, data_sum_1, data_sum_2, delay_params, git_hashes, hubble_param, $
              kperp_lambda_conv, kx_mpc, ky_mpc, kz_mpc, n_freq_contrib, sigma2_1, sigma2_2, $
              sim_noise_diff_1, sim_noise_diff_2, sim_noise_sum_1, sim_noise_sum_2, t_sys_meas, $
              vs_mean, vs_name, window_int, wt_ave_power_freq, wt_meas_ave, wt_meas_min
          endelse
          
        endfor
      endfor
    endfor
  endfor
  
end