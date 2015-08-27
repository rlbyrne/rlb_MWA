;; creates a 2D histogram of the measured and calculated variances

PRO var_compare_scatter_plot, choose_term = choose_term, $
    savefile = savefile, polarization = polarization, use_cubes = use_cubes
    
  kperp_wavelength_max = [0,50,0]
  kperp_wavelength_min = [0,0,50]
  xplotrange = [1e17, 1e28]
  yplotrange = [1e17, 1e28]
  
  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF N_ELEMENTS(kperp_wavelength_min) LT 1 THEN kperp_wavelength_min = 0
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  
  CASE use_cubes OF
    0: BEGIN
      cube_name = 'zenith'
      note_part = 'crossed difference cubes, Aug23 and Aug27 zenith pointings'
    END
    1: BEGIN
      cube_name = 'obssim'
      note_part = 'crossed simulated noise cubes, single obs'
    END
    2: BEGIN
      cube_name = '3hr'
      note_part = 'crossed difference cubes, Aug23 and Aug27 3hr'
    END
    3: BEGIN
      cube_name = 'obs'
      note_part = 'crossed difference cubes, single obs'
    END
  ENDCASE
  
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/MWA/error_analysis_plots/varcompare_scatter_'+cube_name+'_'+polarization+'_term' + STRTRIM(STRING(choose_term), 2), $
    /FONT, XSIZE = 30, YSIZE = 7
    
  ;;SIMULATION:
    
  CROSS_DIFFERENCE, choose_term = choose_term, use_cubes = use_cubes, /diff_sim,$
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    cube_cross = diff_cross_sim, sigma2_correct = sigma2_sim, polarization = polarization, $
    obsid1 = obsid1, obsid2 = obsid2
    
  note = note_part + ', ' + polarization +', term ' + STRTRIM(STRING(choose_term), 2)
  
  sigma2_kperp_sim = MEAN(sigma2_sim, DIMENSION = 3, /NAN)
  var_diff_cross_sim = VARIANCE(REAL_PART(diff_cross_sim), DIMENSION = 3)
  keep_indices = WHERE(sigma2_kperp_sim GT 0 AND var_diff_cross_sim GT 0)
  sigma2_kperp_sim = sigma2_kperp_sim[keep_indices]
  var_diff_cross_sim = var_diff_cross_sim[keep_indices]
  sigma2_kperp_sim_log = ALOG10(sigma2_kperp_sim)
  var_diff_cross_sim_log = ALOG10(var_diff_cross_sim)
  
  IF N_ELEMENTS(xplotrange) GT 0 THEN BEGIN
    min_sigma2 = ALOG10(xplotrange[0])
    max_sigma2 = ALOG10(xplotrange[1])
  ENDIF ELSE BEGIN
    min_sigma2 = MIN(sigma2_kperp_sim_log)
    max_sigma2 = MAX(sigma2_kperp_sim_log)
  ENDELSE
  IF N_ELEMENTS(yplotrange) GT 0 THEN BEGIN
    min_var_diff = ALOG10(yplotrange[0])
    max_var_diff = ALOG10(yplotrange[1])
  ENDIF ELSE BEGIN
    min_var_diff = MIN(var_diff_cross_sim_log)
    max_var_diff = MAX(var_diff_cross_sim_log)
  ENDELSE
  bin_sigma2 = (max_sigma2 - min_sigma2)/100.
  bin_var_diff = (max_var_diff - min_var_diff)/100.
  
  hist_vals_sim = HIST_2D(sigma2_kperp_sim_log, var_diff_cross_sim_log, MIN1 = min_sigma2, MAX1 = max_sigma2, MIN2 = min_var_diff, MAX2 = max_var_diff, $
    BIN1 = bin_sigma2, BIN2 = bin_var_diff)
  ;hist_vals_sim = FLOAT(hist_vals_sim) / TOTAL(hist_vals_sim) * 100
  sigma2_locs = 10.^((FINDGEN((SIZE(hist_vals_sim))[1]) * bin_sigma2) + min_sigma2)
  var_diff_locs = 10.^((FINDGEN((SIZE(hist_vals_sim))[2]) * bin_var_diff) + min_var_diff)
  histvals_max = MAX(hist_vals_sim)
  
  QUICK_IMAGE, hist_vals_sim, sigma2_locs, var_diff_locs, /XLOG, /YLOG, XRANGE = xrange, YRANGE = yrange,  $
    DATA_ASPECT = 1, DATA_RANGE = [0,histvals_max], TITLE = 'Simulation, all kperp', $
    XTITLE = 'calculated variance', YTITLE = 'measured variance', WINDOW = 1,$
    START_MULTI_PARAMS = {nrow:1, ncol:N_ELEMENTS(kperp_wavelength_max) + 1}, MULTI_POS = multi_pos, /LOG, $
    note = note;, $
  ;SAVEFILE = savefile_name, PNG = savefile
  CGPLOT, xrange, xrange, LINESTYLE = 2, /OVERPLOT, /XLOG, /YLOG, XSTYLE = 1, YSTYLE = 1, XRANGE = xrange, YRANGE = yrange
  
  
  ;;DATA:
  
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
  
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], kperp_wavelength_min = kperp_wavelength_min[i], $
      use_cubes = use_cubes, kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
      diff_cross = diff_cross, sigma2_correct = sigma2, polarization = polarization
      
    sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
    var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
    keep_indices = WHERE(sigma2_kperp GT 0 AND var_diff_cross GT 0)
    sigma2_kperp = sigma2_kperp[keep_indices]
    var_diff_cross = var_diff_cross[keep_indices]
    sigma2_kperp_log = ALOG10(sigma2_kperp)
    var_diff_cross_log = ALOG10(var_diff_cross)
    
    hist_vals = HIST_2D(sigma2_kperp_log, var_diff_cross_log, MIN1 = min_sigma2, MAX1 = max_sigma2, MIN2 = min_var_diff, MAX2 = max_var_diff,$
      BIN1 = bin_sigma2, BIN2 = bin_var_diff)
    ;hist_vals = FLOAT(hist_vals) / TOTAL(hist_vals) * 100
      
    IF kperp_wavelength_max[i] EQ 0 THEN BEGIN
      IF kperp_wavelength_min[i] EQ 0 THEN title = 'Data, all kperp' ELSE title = 'Data, kperp > ' + STRTRIM(STRING(kperp_wavelength_min[i]), 2)
    ENDIF ELSE BEGIN
      IF kperp_wavelength_min[i] EQ 0 THEN title = 'Data, kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2) ELSE $
        title = 'Data, ' + STRTRIM(STRING(kperp_wavelength_min[i]), 2) + ' < kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    ENDELSE
    
    QUICK_IMAGE, hist_vals, sigma2_locs, var_diff_locs, /XLOG, /YLOG, XRANGE = xrange, YRANGE = yrange, $
      DATA_ASPECT = 1, DATA_RANGE = [0,histvals_max], TITLE = title, $
      XTITLE = 'calculated variance', YTITLE = 'measured variance', $
      MULTI_POS = multi_pos[*,i+1], /NOERASE, /LOG;, SAVEFILE = savefile_name, PNG = savefile
    CGPLOT, xrange, xrange, LINESTYLE = 2, /OVERPLOT, /XLOG, /YLOG, XSTYLE = 1, YSTYLE = 1, XRANGE = xrange, YRANGE = yrange
    
  ENDFOR
  
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
END