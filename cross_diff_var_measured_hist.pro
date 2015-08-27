PRO cross_diff_var_measured_hist, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, all_kperp = all_kperp

  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1 ;; choose_term defaults to 1
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/varmeasured_crossdiff_hist_term' + STRTRIM(STRING(choose_term), 2), /FONT, XSIZE = 10, YSIZE = 7
  colors = ['black', 'red', 'blue', 'dark green', 'purple', 'maroon']
  legend_colors = []
  legend_titles = []
  
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF KEYWORD_SET(all_kperp) AND N_ELEMENTS(WHERE(kperp_wavelength_max EQ 0, /NULL)) LT 1 THEN kperp_wavelength_max = [0,kperp_wavelength_max]
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
    IF kperp_wavelength_max[i] EQ 0 THEN title = 'all kperp' ELSE title = 'kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
      diff_cross = diff_cross, sigma2_correct = sigma2
      
    var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
    
    IF i EQ 0 THEN do_overplot = 0 ELSE do_overplot = 1
    RLB_HISTPLOT, var_diff_cross, XSTYLE = 1, YSTYLE = 1,  OVERPLOT = do_overplot, $
      COLOR = colors[i], XTITLE = 'Measured Variance, Crossed Difference Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,10], RANGE = [1E15, 1E26], BINSIZE = binsize
    legend_colors = [legend_colors, colors[i]]
    legend_titles = [legend_titles, title]
  ENDFOR
  
  CGLEGEND, TITLE = legend_titles, COLOR = legend_colors, CHARSIZE = 1, LOCATION = [.8, .85], $
    ALIGNMENT = 4
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
  
END