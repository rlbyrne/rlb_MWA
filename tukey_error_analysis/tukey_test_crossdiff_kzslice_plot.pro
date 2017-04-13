pro tukey_test_crossdiff_kzslice_plot, tukey=tukey

  if n_elements(kperp_wavelength_max) lt 1 then kperp_wavelength_max = [0]
  if n_elements(kperp_wavelength_min) lt 1 then kperp_wavelength_min = [0]
  
  if n_elements(xplotrange) lt 1 then xplotrange = [1e5, 1e30]
  if n_elements(yplotrange) lt 1 then yplotrange = [1e5, 1e30]
  
  obsid1 = '1061316176'
  obsid2 = '1061316296'
  use_cube = 'res'
  choose_terms = [1]
  polarizations = ['xx']
  data_range = [1e-6,1e-2]
  
  if keyword_set(tukey) then begin
    note_part = use_cube+' cubes, obsids '+obsid1+' and '+obsid2+', tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endif else begin
    note_part = use_cube+' cubes, obsids '+obsid1+' and '+obsid2+', no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endelse
  
  for term = 0, n_elements(choose_terms) - 1 do begin
    for pol = 0, n_elements(polarizations) - 1 do begin
      print, 'calculating cross difference'
      
      filepath = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/ps/'
      if keyword_set(tukey) then begin
        filename1 = filepath + '1061316176_cubeXX__even_odd_joint_tk_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_tk_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endif else begin
        filename1 = filepath + 'Combined_obs_1061316176_cubeXX__even_odd_joint_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endelse
      
      check_file1 = file_test(filename1)
      check_file2 = file_test(filename2)
      if check_file1 eq 0 or check_file2 eq 0 then begin
        print, '***ERROR*** : files not found'
      endif
      
      if keyword_set(tukey) then begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/crossdiff_ratio_kzslice_tukey_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endif else begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/crossdiff_ratio_kzslice_tophat_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endelse
      
      CROSS_DIFFERENCE, choose_term = choose_terms[term], filename1 = filename1, filename2 = filename2,$
        kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
        cube_cross = diff_cross, sigma2 = var_expected, polarization = polarizations[pol], $
        obsid1 = obsid1, obsid2 = obsid2, /calculate_var, sample_factor = 0, recalc = 1, $
        uvf_input = uvf_input, uv_img_clip, save_cubes = 0
        
      var_measured = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
      var_expected_mean = MEAN(var_expected, DIMENSION = 3, /NAN)
      
      var_ratio = var_measured / var_expected_mean
      where_0 = WHERE(var_expected_mean EQ 0, count_0)
      IF count_0 GT 0 THEN var_ratio[where_0] = 0
      
      QUICK_IMAGE, var_ratio, kx_mpc, ky_mpc, $
        DATA_ASPECT = .5,$
        TITLE = 'Crossed Even/Odd Diff Cubes, Measured/Expected Variance', $
        XTITLE = 'kx', YTITLE = 'ky',$
        note = note_part, /log, $
        data_range = [1e-8,1e-1], savefile = output_loc
        
    ;if keyword_set(tukey) then begin
    ;  output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/crossdiff_varmeasured_tukey_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
    ;endif else begin
    ;  output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/crossdiff_varmeasured_tophat_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
    ;endelse
    ;
    ; QUICK_IMAGE, var_measured, kx_mpc, ky_mpc, $
    ;   DATA_ASPECT = .5,$
    ;   TITLE = 'Crossed Even/Odd Diff Cubes, Measured Variance', $
    ;   XTITLE = 'kx', YTITLE = 'ky',$
    ;   note = note_part, /log, $
    ;   data_range = [1e19,1e28], $
    ;   savefile = output_loc
        
    ;var_ratio = (var_measured - var_expected_mean) / (var_measured + var_expected_mean)^2
    ;where_0 = WHERE((var_measured + var_expected_mean) EQ 0, count_0)
    ;IF count_0 GT 0 THEN var_ratio[where_0] = 0
        
    ;QUICK_IMAGE, var_ratio, kx_mpc, ky_mpc, $
    ;  DATA_ASPECT = .5,$
    ;  TITLE = 'Crossed Even/Odd Diff Cubes, (M-E)/(M+E)^2 Variance', $
    ;  XTITLE = 'kx', YTITLE = 'ky',$
    ;  note = note_part,$
    ;  ;data_range = [-1e-20,1e-20], $
    ;  savefile = output_loc, /log
        
    endfor
  endfor
  
end