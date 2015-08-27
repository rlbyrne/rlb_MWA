;; produces a histogram of the measured, calculated, and simulated variances

PRO cross_diff_simvar_compare_hist, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, all_kperp = all_kperp, $
  polarization = polarization, noise_sim = noise_sim

  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1 ;; choose_term defaults to 1
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/MWA/error_analysis/noise_sim_plots/varcompare_crosssim_hist_' + polarization + '_term' + STRTRIM(STRING(choose_term), 2), /FONT, XSIZE = 10, YSIZE = 7
  colors_measured = ['black','navy', 'red']
  colors_calc = ['slate gray','dodger blue','maroon']
  colors_sim = ['medium gray','cyan', 'deep pink']
  legend_colors = []
  legend_titles = []
  
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF KEYWORD_SET(all_kperp) AND N_ELEMENTS(WHERE(kperp_wavelength_max EQ 0, /NULL)) LT 1 THEN kperp_wavelength_max = [0,kperp_wavelength_max]
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
    IF kperp_wavelength_max[i] EQ 0 THEN title_measured = 'measured variance, all kperp' ELSE title_measured = 'measured variance, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    IF kperp_wavelength_max[i] EQ 0 THEN title_sim = 'simulated variance, all kperp' ELSE title_sim = 'simulated variance, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    IF kperp_wavelength_max[i] EQ 0 THEN title_calc = 'calculated variance, all kperp' ELSE title_calc = 'calculated variance, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      diff_cross = diff_cross, sigma2_correct = sigma2, polarization = polarization, noise_sim = noise_sim
      
    CROSS_DIFF_SIM, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      diff_cross = diff_cross_sim, polarization = polarization, noise_sim = noise_sim
      
    var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
    var_diff_cross_sim = VARIANCE(REAL_PART(diff_cross_sim), DIMENSION = 3)
    sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
    
    binsize = 0.08
    IF i EQ 0 THEN do_overplot = 0 ELSE do_overplot = 1
    RLB_HISTPLOT, var_diff_cross, XSTYLE = 1, YSTYLE = 1,  OVERPLOT = do_overplot, $
      COLOR = colors_measured[i], XTITLE = 'Measured Variance, Crossed Simulated Noise Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,12], RANGE = [1E15, 1E26], BINSIZE = binsize
    RLB_HISTPLOT, var_diff_cross_sim, XSTYLE = 1, YSTYLE = 1,  /OVERPLOT, $
      COLOR = colors_sim[i], XTITLE = 'Measured Variance, Crossed Simulated Noise Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,10], RANGE = [1E15, 1E26], BINSIZE = binsize
    RLB_HISTPLOT, sigma2_kperp, XSTYLE = 1, YSTYLE = 1,  /OVERPLOT, $
      COLOR = colors_calc[i], XTITLE = 'Measured Variance, Crossed Simulated Noise Cubes', /LOGDATA, $
      YTITLE = 'Histogram Count (%)', YRANGE = [0,10], RANGE = [1E15, 1E26], BINSIZE = binsize
    legend_colors = [legend_colors, colors_measured[i], colors_sim[i], colors_calc[i]]
    legend_titles = [legend_titles, title_measured, title_sim, title_calc]
  ENDFOR
  
  CGLEGEND, TITLE = legend_titles, COLOR = legend_colors, CHARSIZE = 1, LOCATION = [.7, .85], $
    ALIGNMENT = 4
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
  
END