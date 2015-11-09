PRO cross_diff_var_measured_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, use_cubes = use_cubes, $
    kperp_wavelength_min = kperp_wavelength_min, polarization = polarization
    
  ;CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, $
  ;  use_cubes = use_cubes, kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
  ;  diff_cross = diff_cross, sigma2_correct = sigma2, polarization = polarization, calculate_var = 0, $
  ;  ACube_sigma2 = ACube_sigma2, BCube_sigma2 = BCube_sigma2
    
  ACube_sigma2 = getvar_savefile('/data4/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_1061316176_0.0002/ps/1061316176_gridded_uvf__even_odd_joint_model_xx_bh_kcube.idlsave', $
    'sigma2_1')
  BCube_sigma2 = getvar_savefile('/data4/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_1061316296_0.0002/ps/1061316296_gridded_uvf__even_odd_joint_model_xx_bh_kcube.idlsave', $
    'sigma2_1')
    
  ;var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  ;var_diff_cross = MEAN(sigma2, DIMENSION = 3)
  ; var_diff_cross_0 = var_diff_cross
  ; var_diff_cross_0[*,*] = 0
  ; var_diff_cross_0[where(diff_cross ne 0)] = 1
    
  ACube_sigma2_kperp = MEAN(ACube_sigma2, DIMENSION = 3, /NAN)
  ACube_sigma2_kperp_0 = ACube_sigma2_kperp
  ACube_sigma2_kperp_0[*,*] = 0
  ACube_sigma2_kperp_0[where(ACube_sigma2_kperp ne 0)] = 1
  BCube_sigma2_kperp = MEAN(BCube_sigma2, DIMENSION = 3, /NAN)
  BCube_sigma2_kperp_0 = BCube_sigma2_kperp
  BCube_sigma2_kperp_0[*,*] = 0
  BCube_sigma2_kperp_0[where(BCube_sigma2_kperp ne 0)] = 1
  sigma2_kperp_0 = ACube_sigma2_kperp_0 + BCube_sigma2_kperp_0
  
  IF KEYWORD_SET(savefile) THEN BEGIN
    filename = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis_plots/varmeasured_crossdiff_UVsim0p0002_UVFinput
    IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN filename = filename + '_kperpmax' + STRTRIM(STRING(kperp_wavelength_max),2)
    QUICK_IMAGE, var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, $
      TITLE = 'Measured Variance > 0, Crossed Difference Cubes', XTITLE = 'kx', YTITLE = 'ky', DATA_RANGE = [0, 1], $
      NOTE = 'crossed simulated noise cubes with uniform 0.0002 UV coverage, single obs, UVF input', $
      SAVEFILE = filename, /PNG
  ENDIF ELSE BEGIN
    QUICK_IMAGE, BCube_sigma2_kperp_0, kx_mpc, ky_mpc, WINDOW = 1, $
      TITLE = 'Measured Variance, Crossed Difference Cubes', XTITLE = 'kx', YTITLE = 'ky';, DATA_RANGE = [0, 1]
  ENDELSE
  stop
  
END