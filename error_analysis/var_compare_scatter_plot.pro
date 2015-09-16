;; creates a 2D histogram of the measured and calculated variances

PRO var_compare_scatter_plot, choose_term = choose_term, $
    savefile = savefile, polarization = polarization, use_cubes = use_cubes, $
    kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, xplotrange = xplotrange, yplotrange = yplotrange, $
    diff_cross = diff_cross, sigma2_all = sigma2_all
    
  if n_elements(kperp_wavelength_max) lt 1 then kperp_wavelength_max = [0]
  if n_elements(kperp_wavelength_min) lt 1 then kperp_wavelength_min = [0]
  if n_elements(xplotrange) lt 1 then xplotrange = [1e20, 1e25]
  if n_elements(yplotrange) lt 1 then yplotrange = [1e20, 1e25]
  
  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
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
    4: BEGIN
      cube_name = 'UVsim0p0002'
      note_part = 'crossed simulated noise cubes with uniform 0.0002 UV coverage, single obs'
    END
    5: BEGIN
      cube_name = 'UVsim0p0005'
      note_part = 'crossed simulated noise cubes with uniform 0.0005 UV coverage, single obs'
    END
    6: BEGIN
      cube_name = 'UVsim0p0002_UVFinput'
      note_part = 'crossed simulated noise cubes with uniform 0.0002 UV coverage, single obs, UVF input'
    END
  ENDCASE
  
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/MWA/error_analysis_plots/varcompare_scatter_'+cube_name+'_'+polarization+'_term' + STRTRIM(STRING(choose_term), 2), $
    /FONT, XSIZE = 8*(n_elements(kperp_wavelength_max)+1), YSIZE = 7
    
  ;;SIMULATION:
  if n_elements(sigma2_all) gt 1 then calculate_var = 0 else calculate_var = 1
  CROSS_DIFFERENCE, choose_term = choose_term, use_cubes = use_cubes, /diff_sim,$
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    cube_cross = diff_cross_sim, sigma2_correct = sigma2_output, polarization = polarization, $
    obsid1 = obsid1, obsid2 = obsid2, calculate_var = calculate_var
  if calculate_var eq 1 then sigma2_all = sigma2_output
  
  note = note_part + ', ' + polarization +', term ' + STRTRIM(STRING(choose_term), 2)
  
  sigma2_kperp_all = MEAN(sigma2_all, DIMENSION = 3, /NAN)
  var_diff_cross_sim = VARIANCE(REAL_PART(diff_cross_sim), DIMENSION = 3)
  keep_indices = WHERE(sigma2_kperp_all GT 0 AND var_diff_cross_sim GT 0)
  sigma2_kperp_sim = sigma2_kperp_all[keep_indices]
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
  stop
  
  ;;DATA:
  
  if n_elements(diff_cross) eq 0 then begin
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = 0, kperp_wavelength_min = 0, $
      use_cubes = use_cubes, diff_cross = diff_cross, polarization = polarization, calculate_var = 0
  endif
  var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
  
    IF kperp_wavelength_max NE 0 THEN BEGIN
      kx_indices = WHERE(ABS(kx_mpc) LE kperp_wavelength_max[i]/kperp_lambda_conv AND ABS(kx_mpc) GT kperp_wavelength_min[i]/kperp_lambda_conv, /NULL)
      ky_indices = WHERE(ABS(ky_mpc) LE kperp_wavelength_max[i]/kperp_lambda_conv AND ABS(kx_mpc) GT kperp_wavelength_min[i]/kperp_lambda_conv, /NULL)
    ENDIF ELSE BEGIN
      kx_indices = WHERE(ABS(kx_mpc) GT kperp_wavelength_min[i]/kperp_lambda_conv, /NULL)
      ky_indices = WHERE(ABS(ky_mpc) GT kperp_wavelength_min[i]/kperp_lambda_conv, /NULL)
    ENDELSE
    
    sigma2_kperp_use = sigma2_kperp_all[kx_indices, ky_indices, *]
    var_diff_cross_use = var_diff_cross[kx_indices, ky_indices, *]
    keep_indices = WHERE(sigma2_kperp_use GT 0 AND var_diff_cross_use GT 0)
    sigma2_kperp_use = sigma2_kperp_use[keep_indices]
    var_diff_cross_use = var_diff_cross_use[keep_indices]
    sigma2_kperp_log = ALOG10(sigma2_kperp_use)
    var_diff_cross_log = ALOG10(var_diff_cross_use)
    
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
  
  stop
  
END