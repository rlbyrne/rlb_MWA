PRO cross_diff_var_measured_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile

  cross_difference, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, $ ;;choose_term defaults to 1
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    diff_cross = diff_cross, sigma2_correct = sigma2
    
  var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  
  IF KEYWORD_SET(savefile) THEN BEGIN
    filename = '/nfs/eor-00/h1/rbyrne/error_analysis/varmeasured_crossdiff_term' + STRTRIM(STRING(choose_term),2)
    IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN filename = filename + '_kperpmax' + STRTRIM(STRING(kperp_wavelength_max),2)
    QUICK_IMAGE, var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, $
      TITLE = 'Measured Variance, Crossed Difference Cubes', XTITLE = 'kx', YTITLE = 'ky', DATA_RANGE = [0, 2*STDDEV(var_diff_cross)], $
      NOTE = 'Long Run Aug23 and Aug27, Term ' + STRTRIM(STRING(choose_term),2), $
      SAVEFILE = filename, /PNG
  ENDIF ELSE BEGIN
    QUICK_IMAGE, var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, $
      TITLE = 'Measured Variance, Crossed Difference Cubes', XTITLE = 'kx', YTITLE = 'ky', DATA_RANGE = [0, 2*STDDEV(var_diff_cross)],$
      NOTE = 'Long Run Aug23 and Aug27, Term ' + STRTRIM(STRING(choose_term),2)
  ENDELSE
  
  
END