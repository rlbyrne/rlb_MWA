PRO cross_diff_var_plot, savefile = savefile, $
    kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, $
    where_0 = where_0, recalc = recalc, uvf_input = uvf_input, uv_img_clip = uv_img_clip, manual_uv_img_clip = manual_uv_img_clip
    
  if n_elements(kperp_wavelength_max) lt 1 then kperp_wavelength_max = [0]
  if n_elements(kperp_wavelength_min) lt 1 then kperp_wavelength_min = [0]
  plot_max = 20.
  
  ;; for simulated noise cubes with uniform UV coverage
  obsid1 = '1061316176'
  obsid2 = '1061316296'
  ;sample_factors = [.0002,.0005,.001,.005,.01,.05,.1,.5,1,5]
  sample_factors = [5]
  choose_terms = [1]
  polarizations = ['xx']
  data_range = [0,1e3]
  sample_factors_str = num_formatter_filename(sample_factors)
  cube_names = 'UVsim' + sample_factors_str
  note_part = 'crossed sim noise with ' + number_formatter(sample_factors) + ' UV coverage, single obs'
  
  if keyword_set(uvf_input) then begin
    cube_names = cube_names + '_UVFinput'
    note_part = note_part + ', UVF input'
    if keyword_set(uv_img_clip) then begin
      cube_names = cube_names + '_UVimgclip' + number_formatter(uv_img_clip)
      note_part = note_part + ', UV img clip ' + number_formatter(uv_img_clip)
    endif
    if keyword_set(manual_uv_img_clip) then begin
      cube_names = cube_names + '_manualUVimgclip' + number_formatter(manual_uv_img_clip)
      note_part = note_part + ', Manual UV img clip ' + number_formatter(manual_uv_img_clip)
    endif
  endif else cube_names = cube_names + '_Heal'
  
  for sample = 0, n_elements(sample_factors) -1 do begin
    for term = 0, n_elements(choose_terms) - 1 do begin
      for pol = 0, n_elements(polarizations) - 1 do begin
      
        filename1 = '/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_' + number_formatter(sample_factors[sample]) + '/ps/'+ obsid1
        filename2 = '/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_' + number_formatter(sample_factors[sample]) + '/ps/'+ obsid2
        if keyword_set(uvf_input) then begin
          if ~keyword_set(uv_img_clip) then begin
            filename1 = filename1 + '_gridded_uvf__even_odd_joint_model_' + polarizations[pol] + '_bh_kcube.idlsave'
            filename2 = filename2 + '_gridded_uvf__even_odd_joint_model_' + polarizations[pol] + '_bh_kcube.idlsave'
          endif else begin
            filename1 = filename1 + '_gridded_uvf__even_odd_joint_uvimgclip'+ num_formatter_filename(uv_img_clip) + '_model_' + polarizations[pol] + '_bh_kcube.idlsave'
            filename2 = filename2 + '_gridded_uvf__even_odd_joint_uvimgclip' + num_formatter_filename(uv_img_clip) + '_model_' + polarizations[pol] + '_bh_kcube.idlsave'
          endelse
        endif else begin
          filename1 = filename1 + '_cubeXX__even_odd_joint_model_' + polarizations[pol] + '_bh_kcube.idlsave'
          filename2 = filename2 + '_cubeXX__even_odd_joint_model_' + polarizations[pol] + '_bh_kcube.idlsave'
          undefine, uv_img_clip
        endelse
        
        check_file1 = file_test('/data3' + filename1)
        check_file2 = file_test('/data3' + filename2)
        if check_file1 eq 1 and check_file2 eq 1 then begin
          filename1 = '/data3' + filename1
          filename2 = '/data3' + filename2
        endif else begin
          print, '***ERROR*** : files not found'
          continue
        endelse
        
        if keyword_set(savefile) then output_loc = '/home/rlbyrne/error_analysis_plots/varratio_plot_'+cube_names[sample]+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
        note = note_part[sample] + ', ' + polarizations[pol] +', term ' + STRTRIM(STRING(choose_terms[term]), 2)
        
        if keyword_set(recalc) then recalc_use = 1 else recalc_use = 0
        CROSS_DIFFERENCE, choose_term = choose_terms[term], filename1 = filename1, filename2 = filename2, $
          kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
          cube_cross = diff_cross, sigma2 = sigma2, polarization = polarizations[pol], $
          obsid1 = obsid1, obsid2 = obsid2, /calculate_var, sample_factor = sample_factors_str[sample], recalc = recalc_use, $
          uvf_input = uvf_input, uv_img_clip
          
        var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
        sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
        
        ratio_var_diff_cross = var_diff_cross / sigma2_kperp
        wh_sig0 = WHERE(sigma2_kperp EQ 0, count_sig0)
        IF count_sig0 GT 0 THEN ratio_var_diff_cross[wh_sig0] = 0
        
        stop
        
        QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, DATA_RANGE = [0,plot_max], $
          TITLE = 'Measured / Expected Variance', XTITLE = 'kx', YTITLE = 'ky', $
          NOTE = note, SAVEFILE = output_loc, PNG=savefile, /log, START_MULTI_PARAMS = {nrow:1, ncol:3}, $
          multi_pos = multi_pos, no_ps_close=savefile, noerase=savefile
        QUICK_IMAGE, var_diff_cross, kx_mpc, ky_mpc, DATA_RANGE = [0,plot_max], $
          TITLE = 'Measured Variance', XTITLE = 'kx', YTITLE = 'ky', $
          NOTE = note, SAVEFILE = output_loc, PNG=savefile, /log, multi_pos = multi_pos[*,1], $
          noerase=savefile, no_ps_close=savefile
        QUICK_IMAGE, sigma2_kperp, kx_mpc, ky_mpc, DATA_RANGE = [0,plot_max], $
          TITLE = 'Expected Variance', XTITLE = 'kx', YTITLE = 'ky', $
          NOTE = note, SAVEFILE = output_loc, PNG=savefile, /log, multi_pos = multi_pos[*,2]
          
      endfor
    endfor
  endfor
  
END