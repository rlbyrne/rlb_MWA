PRO cross_diff_weights_plot, savefile = savefile, uv_img_clip = uv_img_clip, uvf_input = uvf_input

  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF N_ELEMENTS(kperp_wavelength_min) LT 1 THEN kperp_wavelength_min = 0
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  
  cube_name = 'UVsim0p0005'
  note_part = 'crossed simulated noise cubes with uniform 0.0005 UV coverage, single obs'
  obsid1 = '1061316176'
  obsid2 = '1061316296'
  density = 1
  
  filename1 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_' + number_formatter(density) + '/ps/'+ obsid1
  filename2 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_' + number_formatter(density) + '/ps/'+ obsid2
  if keyword_set(uvf_input) then begin
    if ~keyword_set(uv_img_clip) then begin
      filename1 = filename1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = filename2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    endif else begin
      filename1 = filename1 + '_gridded_uvf__even_odd_joint_uvimgclip'+ num_formatter_filename(uv_img_clip) + '_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = filename2 + '_gridded_uvf__even_odd_joint_uvimgclip' + num_formatter_filename(uv_img_clip) + '_model_' + polarization + '_bh_kcube.idlsave'
    endelse
  endif else begin
    filename1 = filename1 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    filename2 = filename2 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
  endelse
  
  data_range = [0, 4e-1]
  ACube_weights = getvar_savefile(filename1, 'WT_MEAS_AVE')
  BCube_weights = getvar_savefile(filename2, 'WT_MEAS_AVE')
  ;ACube_sigma2 = getvar_savefile(filename1, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  ;ACube_sigma2_kperp = MEAN(ACube_sigma2, DIMENSION = 3, /NAN)
  ;BCube_sigma2 = getvar_savefile(filename2, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  ;BCube_sigma2_kperp = MEAN(BCube_sigma2, DIMENSION = 3, /NAN)
  
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  
  ;uvf_weights = getvar_savefile('/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/1061316176_even_cubeXX_weights_uvf.idlsave', 'WEIGHTS_CUBE')
  ;uvf_weights = mean(uvf_weights, dimension = 3, /nan)
  
  ;weights_ptr = getvar_savefile('/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/1061316176_even_gridded_uvf.sav', 'weights_uv_arr')
  ;weights = make_array(1200,1200,192,/complex)
  ;pol = 0
  ;for freq = 0,191 do begin
  ;  weights[*, *, freq] = *weights_ptr[pol,freq]
  ;endfor
  
  
  
  
  ;  IF KEYWORD_SET(savefile) THEN BEGIN
  ;    filename = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis_plots/varratio_' + cube_name + '_' + polarization + '_term' + STRTRIM(STRING(choose_term),2)
  ;    note = note_part + ', ' + polarization +', term ' + STRTRIM(STRING(choose_term), 2) ;+ ': Mean = ' + STRMID(STRTRIM(STRING(MEAN(ratio_var_diff_cross, /NAN)), 2), 0, 4)
  ;    IF kperp_wavelength_max NE 0 THEN filename = filename + '_kperpmax' + STRTRIM(STRING(kperp_wavelength_max),2)
  ;    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,plot_max], $
  ;      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky', $
  ;      NOTE = note, SAVEFILE = filename, /PNG
  ;  ENDIF ELSE BEGIN
  ;    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,plot_max], $
  ;      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky', $
  ;      NOTE = note
  ;  ENDELSE
  
  output_path = '/home/rlbyrne/error_analysis_plots/'
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, output_path + 'varcompare_scatter_'+cube_name+'_'+polarization+'_term' + STRTRIM(STRING(choose_term), 2), $
    /FONT, XSIZE = 8*(n_elements(kperp_wavelength_max)+1), YSIZE = 7
    
  QUICK_IMAGE, ACube_weights, kx_mpc, ky_mpc,  $
    DATA_ASPECT = 0.5, TITLE = 'Weights, Obs ' + obsid1, $
    XTITLE = 'kx', YTITLE = 'ky', data_range = data_range, $
    START_MULTI_PARAMS = {nrow:3, ncol:1}, MULTI_POS = multi_pos, $
    note = 'Flat UV coverage ' + number_formatter(density)
    
  quick_image, BCube_weights, kx_mpc, ky_mpc,  $
    DATA_ASPECT = 0.5, TITLE = 'Weights, Obs ' + obsid2, $
    XTITLE = 'kx', YTITLE = 'ky', data_range = data_range, $
    multi_pos = multi_pos[*,1], /noerase
    
  quick_image, ACube_weights - BCube_weights, kx_mpc, ky_mpc, $
    DATA_ASPECT = 0.5, TITLE = 'Difference', $
    XTITLE = 'kx', YTITLE = 'ky', data_range = [-0.01,0.01], $
    multi_pos = multi_pos[*,2], /noerase
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
END