;; creates 2D density plots of where in kperp the ratio measured / calculated variance falls within a given value range

PRO kperp_var_ratio, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, $
    xplotrange = xplotrange, yplotrange = yplotrange, savefile = savefile
    
  thresholds = [0,0.2,0.4,1,1000]
  
  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF N_ELEMENTS(kperp_wavelength_min) LT 1 THEN kperp_wavelength_min = 0
  
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/kperp_var_ratio_term' + STRTRIM(STRING(choose_term), 2), $
    /FONT, XSIZE = 50, YSIZE = 20
    
  CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, $
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    diff_cross = diff_cross, sigma2_correct = sigma2
    
  sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
  var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  sigma2_ratio = var_diff_cross / sigma2_kperp
  wh_sig0 = WHERE(sigma2_kperp EQ 0, count_sig0)
  IF count_sig0 GT 0 THEN sigma2_ratio[wh_sig0] = 0
  
  binsize = MAX(kx_mpc)/50.
  
  FOR interval = 0, N_ELEMENTS(thresholds) - 2 DO BEGIN
  
    position_matrix = MAKE_ARRAY(N_ELEMENTS(kx_mpc), N_ELEMENTS(ky_mpc), /FLOAT, VALUE = 0)
    position_matrix[WHERE(sigma2_ratio GT thresholds[interval] AND sigma2_ratio LT thresholds[interval+1], /NULL)] = 1.
    
    kx_hist = HISTOGRAM(kx_mpc, BINSIZE = binsize, REVERSE_INDICES = ri_kx, LOCATIONS = kx_locs)
    ky_hist = HISTOGRAM(ky_mpc, BINSIZE = binsize, REVERSE_INDICES = ri_ky, LOCATIONS = ky_locs)
    plotvals = MAKE_ARRAY(N_ELEMENTS(kx_locs), N_ELEMENTS(ky_locs), VALUE = 0)
    FOR bin_x = 0, N_ELEMENTS(kx_locs) - 1 DO BEGIN
      FOR bin_y = 0, N_ELEMENTS(ky_locs) - 1 DO BEGIN
        plotvals[bin_x, bin_y] = TOTAL(position_matrix[ri_kx[ri_kx[bin_x]:ri_kx[bin_x+1]-1],ri_ky[ri_ky[bin_y]:ri_ky[bin_y+1]-1],0])
      ENDFOR
    ENDFOR
    
    IF interval EQ 0 THEN start_multi_params = {ncol:2, nrow:2} ELSE multi_pos = multi_pos_init[*,interval]
    IF interval EQ 0 THEN title_exten = '< ' + STRTRIM(STRING(thresholds[interval + 1], FORMAT = '(F5.2)'), 2) ELSE BEGIN
      IF interval EQ N_ELEMENTS(thresholds) - 2 THEN title_exten = '> ' + STRTRIM(STRING(thresholds[interval], FORMAT = '(F5.2)'), 2) ELSE $
        title_exten = STRTRIM(STRING(thresholds[interval], FORMAT = '(F5.2)'), 2) + ' to ' + STRTRIM(STRING(thresholds[interval+1], FORMAT = '(F5.2)'), 2)
    ENDELSE
        
    QUICK_IMAGE, plotvals, kx_locs, ky_locs, $
      TITLE = 'Measured / Calculated Variance ' + title_exten, XTITLE = 'kx', YTITLE = 'ky', $
      START_MULTI_PARAMS = start_multi_params, MULTI_POS = multi_pos, NOERASE = no_erase, DATA_RANGE = data_range
      
    
    IF interval EQ 0 THEN BEGIN
      multi_pos_init = multi_pos
      UNDEFINE, start_multi_params
      no_erase = 1
    ENDIF
    
  ENDFOR  
  
  
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
END