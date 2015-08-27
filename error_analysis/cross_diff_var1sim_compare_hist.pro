PRO cross_diff_var1sim_compare_hist, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, all_kperp = all_kperp

  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1 ;; choose_term defaults to 1
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/varcompare_var1sim_crossdiff_hist_term' + STRTRIM(STRING(choose_term), 2), /FONT, XSIZE = 10, YSIZE = 7
  colors_measured = ['black', 'blue', 'dark green']
  colors_calc = ['gray', 'cornflower blue', 'green']
  legend_colors = []
  legend_titles = []
  
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF KEYWORD_SET(all_kperp) AND N_ELEMENTS(WHERE(kperp_wavelength_max EQ 0, /NULL)) LT 1 THEN kperp_wavelength_max = [0,kperp_wavelength_max]
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
    IF kperp_wavelength_max[i] EQ 0 THEN title_measured = 'normalized variance, all kperp' ELSE title_measured = 'normalized variance, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    IF kperp_wavelength_max[i] EQ 0 THEN title_calc = 'simulated normalized variance, all kperp' ELSE title_calc = 'simulated normalized variance, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], /NORMALIZE, $
      kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
      diff_cross = diff_cross, sigma2_correct = sigma2 ;; cross_difference script is altered so that the difference cubes are divided by their variance cubes
      
    CROSS_DIFF_VAR1SIM, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
      diff_cross = diff_cross_sim
      
    var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
    print, mean(var_diff_cross)
    title_measured = title_measured + ', mean = ' + STRTRIM(CGNUMBER_FORMATTER(MEAN(var_diff_cross, /NAN), DECIMAL = 2),2)
    var_diff_cross_sim = VARIANCE(REAL_PART(diff_cross_sim), DIMENSION = 3)
    print, mean(var_diff_cross_sim)
    title_calc = title_calc + ', mean = ' + STRTRIM(CGNUMBER_FORMATTER(MEAN(var_diff_cross_sim, /NAN), DECIMAL = 2),2)
    
    IF i EQ 0 THEN do_overplot = 0 ELSE do_overplot = 1
    RLB_HISTPLOT, var_diff_cross, XSTYLE = 1, YSTYLE = 1,  OVERPLOT = do_overplot, $
      COLOR = colors_measured[i], XTITLE = 'Variance, Crossed Difference Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,20], RANGE = [.01,5], BINSIZE = .1
    RLB_HISTPLOT, var_diff_cross_sim, XSTYLE = 1, YSTYLE = 1,  /OVERPLOT, $
      COLOR = colors_calc[i], XTITLE = 'Variance, Crossed Difference Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,20], RANGE = [.01,5], BINSIZE = .1
    legend_colors = [legend_colors, colors_measured[i], colors_calc[i]]
    legend_titles = [legend_titles, title_measured, title_calc]
  ENDFOR
  
  CGLEGEND, TITLE = legend_titles, COLOR = legend_colors, CHARSIZE = 1, LOCATION = [.6, .85], $
    ALIGNMENT = 4
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
  
END