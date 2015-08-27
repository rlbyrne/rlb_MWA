PRO kz_var_2dplot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, day = day

  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1 ;default term 1
  IF N_ELEMENTS(day) EQ 0 THEN day = 1 ;default Aug23
  IF N_ELEMENTS(kperp_wavelength_max) EQ 0 THEN kperp_wavelength_max = 0 ;default all kperp
  IF day EQ 1 THEN day_name = 'Aug23' ELSE BEGIN
    IF day EQ 2 THEN day_name = 'Aug27' ELSE BEGIN
      PRINT, 'invalid day input'
      RETURN
    ENDELSE
  ENDELSE
  
  IF KEYWORD_SET(savefile) THEN savefile_name = '/nfs/eor-00/h1/rbyrne/error_analysis/kz_var_2dplot_log_' + day_name + '_term' + STRTRIM(STRING(choose_term), 2)
  
  CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, $
    Aug23_sigma2 = Aug23_sigma2, Aug27_sigma2 = Aug27_sigma2, $
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc
    
    
  IF day EQ 2 THEN sigma2 = Aug27_sigma2 ELSE sigma2 = Aug23_sigma2
  sigma2_ratio = sigma2
  FOR j = 0, N_ELEMENTS(kz_mpc) - 1 DO BEGIN
    sigma2_ratio_slice = sigma2[*,*,j] / MEAN(sigma2, DIMENSION = 3)
    wh_sig0 = WHERE(MEAN(sigma2, DIMENSION = 3) EQ 0, count_sig0)
    IF count_sig0 GT 0 THEN sigma2_ratio_slice[wh_sig0] = 0
    sigma2_ratio[*,*,j] = sigma2_ratio_slice
  ENDFOR
  
  sigma2_ratio_2d = kspace_rebinning_2d(sigma2_ratio, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, kperp_bin = kperp_bin, kpar_bin = kpar_bin)
  
  delay_params = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
    'DELAY_PARAMS')
    
  ;QUICK_IMAGE, sigma2_ratio_2d, kperp_edges_mpc, kpar_edges_mpc, /xlog, /ylog, $
  ;  DATA_ASPECT = 1, DATA_RANGE = [0,2], TITLE = 'Difference Cube, Variance / Mean Variance in kz', $
  ;  XTITLE = 'kperp', YTITLE = 'kpar', PNG = savefile, SAVEFILE = savefile_name, $
  ;  NOTE = 'Long Run, Zenith Pointing ' + day_name + ', Term' + STRTRIM(STRING(choose_term), 2)
    
  kpower_2d_plots, power = sigma2_ratio_2d, delay_params = delay_params, $
    kperp_edges = kperp_edges_mpc, kpar_edges = kpar_edges_mpc, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
    plotfile = savefile_name, PNG = savefile, FULL_TITLE = 'Difference Cube, Variance / Mean Variance in kz', $
    NOTE = 'Long Run, Zenith Pointing ' + day_name + ', Term' + STRTRIM(STRING(choose_term), 2), DATA_RANGE = [0,2]
  
  
END