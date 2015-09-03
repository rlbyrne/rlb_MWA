PRO cross_diff_weights_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, $
    polarization = polarization, use_cubes = use_cubes
    
  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF N_ELEMENTS(kperp_wavelength_min) LT 1 THEN kperp_wavelength_min = 0
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  plot_max = 0.0005
  
  if n_elements(use_cubes) eq 0 then use_cubes = 5
  
  CASE use_cubes OF
    4: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0002
      cube_name = 'UVsim0p0002'
      note_part = 'crossed simulated noise cubes with uniform 0.0002 UV coverage, single obs'
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid1 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid2 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    5: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0005
      cube_name = 'UVsim0p0005'
      note_part = 'crossed simulated noise cubes with uniform 0.0005 UV coverage, single obs'
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid1 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid2 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
    END
  ENDCASE
  
  ACube_weights = getvar_savefile(filename1, 'WT_MEAS_AVE')
  BCube_weights = getvar_savefile(filename2, 'WT_MEAS_AVE')
  ACube_sigma2 = getvar_savefile(filename1, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  ACube_sigma2_kperp = MEAN(ACube_sigma2, DIMENSION = 3, /NAN)
  BCube_sigma2 = getvar_savefile(filename2, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  BCube_sigma2_kperp = MEAN(BCube_sigma2, DIMENSION = 3, /NAN)
  
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  
  uvf_weights = getvar_savefile('/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/1061316176_even_cubeXX_weights_uvf.idlsave', 'WEIGHTS_CUBE')
  uvf_weights = mean(uvf_weights, dimension = 3, /nan)
  
  weights_ptr = getvar_savefile('/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/1061316176_even_gridded_uvf.sav', 'weights_uv_arr')
  weights = make_array(1200,1200,192,/complex)
  pol = 0
  for freq = 0,191 do begin
    weights[*, *, freq] = *weights_ptr[pol,freq]
  endfor
  
  stop
  
  
  
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
  
  quick_image, acube_weights, kx_mpc, ky_mpc, window = 0
  quick_image, bcube_weights, kx_mpc, ky_mpc, window = 1
  quick_image, abs(uvf_weights), kx_mpc, ky_mpc, window = 2
  stop
  
END