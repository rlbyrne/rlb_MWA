PRO cross_diff_var_ratio_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, polarization = polarization, noise_sim = noise_sim

  cross_difference, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, $ ;;choose_term defaults to 1
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    diff_cross = diff_cross, sigma2_correct = sigma2, polarization = polarization, noise_sim = noise_sim
    
  var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
  
  ratio_var_diff_cross = var_diff_cross / sigma2_kperp
  wh_sig0 = WHERE(sigma2_kperp EQ 0, count_sig0)
  IF count_sig0 GT 0 THEN ratio_var_diff_cross[wh_sig0] = 0
  
  IF KEYWORD_SET(savefile) THEN BEGIN
    filename = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis/noise_sim_plots/varratio_crosssim_' + polarization + '_term' + STRTRIM(STRING(choose_term),2)
    IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN filename = filename + '_kperpmax' + STRTRIM(STRING(kperp_wavelength_max),2)
    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,2], $
      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky', $
      ;NOTE = 'Long Run Aug23 and Aug27, Term ' + STRTRIM(STRING(choose_term),2) + ': Mean = ' + STRMID(STRTRIM(STRING(MEAN(ratio_var_diff_cross, /NAN)), 2), 0, 4), $
      SAVEFILE = filename, /PNG
  ENDIF ELSE BEGIN
    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,2], $
      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky';, $
      ;NOTE = 'Long Run Aug23 and Aug27, Term ' + STRTRIM(STRING(choose_term),2) + ': Mean = ' + STRMID(STRTRIM(STRING(MEAN(ratio_var_diff_cross, /NAN)), 2), 0, 4)
  ENDELSE
  
  
END