;; creates a 2D histogram of the measured and calculated variances

PRO tukey_test_uvf_crossdiff_scatter_plot, $
    savefile = savefile, $
    kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, xplotrange = xplotrange, yplotrange = yplotrange, $
    where_0 = where_0, recalc = recalc, tukey = tukey
    
  if n_elements(kperp_wavelength_max) lt 1 then kperp_wavelength_max = [0]
  if n_elements(kperp_wavelength_min) lt 1 then kperp_wavelength_min = [0]
  
  if n_elements(xplotrange) lt 1 then xplotrange = [1e7, 1e30]
  if n_elements(yplotrange) lt 1 then yplotrange = [1e7, 1e30]
  
  obsid1 = '1061316176'
  obsid2 = '1061316296'
  choose_terms = [1]
  polarizations = ['yy']
  data_range = [0,1e2]
  if keyword_set(tukey) then begin
    note_part = 'tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endif else begin
    note_part = 'no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endelse
  
  enterprise = 0
  
  for term = 0, n_elements(choose_terms) - 1 do begin
    for pol = 0, n_elements(polarizations) - 1 do begin
      print, 'calculating cross difference'
      
      filepath = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/ps/'
      if keyword_set(tukey) then begin
        filename1 = filepath + '1061316176_cubeXX__even_odd_joint_tk_model_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_tk_model_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endif else begin
        filename1 = filepath + 'Combined_obs_1061316176_cubeXX__even_odd_joint_model_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_model_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endelse
      
      check_file1 = file_test(filename1)
      check_file2 = file_test(filename2)
      if check_file1 eq 0 or check_file2 eq 0 then begin
        print, '***ERROR*** : files not found'
      endif
      
      if keyword_set(tukey) then begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/varcompare_tukey_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endif else begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/varcompare_tophat_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endelse
      
      IF KEYWORD_SET(savefile) THEN CGPS_OPEN, output_loc, $
        /FONT, XSIZE = 8*(n_elements(kperp_wavelength_max)+1), YSIZE = 7
        
      if keyword_set(recalc) then recalc_use = 1 else recalc_use = 0
      CROSS_DIFFERENCE, choose_term = choose_terms[term], filename1 = filename1, filename2 = filename2, /diff_sim,$
        kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
        sim_cube_cross = diff_cross_sim, cube_cross = diff_cross, sigma2 = sigma2_all, polarization = polarizations[pol], $
        obsid1 = obsid1, obsid2 = obsid2, /calculate_var, sample_factor = 0, recalc = recalc_use, $
        uvf_input = uvf_input, uv_img_clip, save_cubes = 0
        
      note = note_part
      
      ;; PLOT SIMULATION:
      
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
      sigma2_locs = 10.^((FINDGEN((SIZE(hist_vals_sim))[1]) * bin_sigma2) + min_sigma2)
      var_diff_locs = 10.^((FINDGEN((SIZE(hist_vals_sim))[2]) * bin_var_diff) + min_var_diff)
      histvals_max = MAX(hist_vals_sim)
      
      if n_elements(data_range) eq 0 then begin
        if keyword_set(where_0) then begin
          hist_vals_sim[where(hist_vals_sim ne 0, /null)] = 1
          data_range = [0,1]
        endif else begin
          data_range = [0, histvals_max]
        endelse
      endif
      
      QUICK_IMAGE, hist_vals_sim, sigma2_locs, var_diff_locs, /XLOG, /YLOG, XRANGE = xrange, YRANGE = yrange,  $
        DATA_ASPECT = 1, DATA_RANGE = data_range, TITLE = 'Simulation, all kperp', $
        XTITLE = 'expected variance', YTITLE = 'measured variance', WINDOW = 1,$
        START_MULTI_PARAMS = {nrow:1, ncol:N_ELEMENTS(kperp_wavelength_max) + 1}, MULTI_POS = multi_pos, /LOG, $
        note = note
      CGPLOT, xrange, xrange, LINESTYLE = 2, /OVERPLOT, /XLOG, /YLOG, XSTYLE = 1, YSTYLE = 1, XRANGE = xrange, YRANGE = yrange
      
      ;;PLOT DATA:
      
      var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
      
      FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
      
        IF kperp_wavelength_max[i] NE 0 THEN BEGIN
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
        
        if keyword_set(where_0) then hist_vals[where(hist_vals ne 0, /null)] = 1
        
        QUICK_IMAGE, hist_vals, sigma2_locs, var_diff_locs, /XLOG, /YLOG, XRANGE = xrange, YRANGE = yrange, $
          DATA_ASPECT = 1, DATA_RANGE = data_range, TITLE = title, $
          XTITLE = 'expected variance', YTITLE = 'measured variance', $
          MULTI_POS = multi_pos[*,i+1], /NOERASE, /LOG;, SAVEFILE = savefile_name, PNG = savefile
        CGPLOT, xrange, xrange, LINESTYLE = 2, /OVERPLOT, /XLOG, /YLOG, XSTYLE = 1, YSTYLE = 1, XRANGE = xrange, YRANGE = yrange
        
      ENDFOR
      
      IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
      undefine, multi_pos
      
    endfor
  endfor
  
END