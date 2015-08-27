PRO var_prop_error, sigma2_1, sigma2_2, var_errors, kperp_wavelength_max = kperp_wavelength_max

  ;; How bad is your variance propagation estimate anyway?
  ;; If you're multiplying A*CONJ(B) and propagating the variance like (VarA + VarB)^2 / 2
  ;; but VarA NE VarB, then your variance propagation is wrong! This procedure graphically shows
  ;; by how much.



  RESTORE, '~/error_analysis/var_prop_corrections.sav'
  xvals = var_ratio
  var_errors = correction_factor
  
  
  IF N_ELEMENTS(var_errors) EQ 0 THEN BEGIN
    xvals = (FINDGEN(100)+.5)/100
    var_ref = 1
    var_errors = MAKE_ARRAY(N_ELEMENTS(xvals), VALUE = 0.)
    FOR i = 0, N_ELEMENTS(xvals) - 1 DO BEGIN
      cross_difference_debug, var_ref, var_ref/(1/xvals[i] - 1), error, kx_mpc, ky_mpc, kz_mpc
      var_errors[i] = error
    ENDFOR
  ENDIF
  
  ymax = 3
  
  CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/var_prop_error', /FONT, XSIZE = 10, YSIZE = 7
  
  CGPLOT, xvals, var_errors, XTITLE = 'Variance 2 / (Variance 1 + Variance 2)', $
    YTITLE = 'Estimated Variance / Actual Variance', XRANGE = [0,1], YRANGE = [1,ymax], $
    TITLE = 'Error in Propagated Variances'
    
  legend_title = ['Error']
  legend_color = ['black']
  
  
  IF N_ELEMENTS(sigma2_1) EQ 0 AND N_ELEMENTS(sigma2_2) EQ 0 THEN BEGIN
    PRINT, 'Using variance cubes from Aug23 and Aug27, long run'
    sigma2_1 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'SIGMA2_1')
    sigma2_2 = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'SIGMA2_1')
  ENDIF
  
  sigma2_diff = sigma2_2 / (sigma2_1 + sigma2_2)
  
  ;PRINT, 'MEAN:' + STRING(MEAN(sigma2_diff, /NAN))
  ;PRINT, 'CORRECTION FACTOR (MEAN):' + STRING(var_prop_correct(1 - MEAN(sigma2_diff, /nan), MEAN(sigma2_diff, /nan)))
  ;PRINT, 'STANDARD DEV:' + STRING(STDDEV(sigma2_diff, /NAN))
  
  binsize = 1/200.
  histvals = HISTOGRAM(sigma2_diff, BINSIZE = binsize, MIN = 0, MAX = 1, /NAN)
  histvals = histvals / TOTAL(histvals) * 100
  hist_x = FINDGEN(N_ELEMENTS(histvals)) * binsize
  histvals_scaled = histvals / MAX(histvals) * .8 * (ymax - 1) + 1
  
  CGOPLOT, hist_x, histvals_scaled, PSYM = 10, COLOR = 'red'
  
  new_legend_title = ['Variance Distribution (Difference Cubes, Aug23 and Aug27 Long Run)']
  legend_title = [legend_title, new_legend_title]
  new_legend_color = ['red']
  legend_color = [legend_color, new_legend_color]
  
  
  IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN BEGIN
  
    kperp_lambda_conv = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KPERP_LAMBDA_CONV')
    kx_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KX_MPC')
    ky_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KY_MPC')
      
    kx_indices = WHERE(ABS(kx_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
    ky_indices = WHERE(ABS(ky_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
    ky_indices = ky_indices[1:*] ;remove first bin with zero values
    
    sigma2_1 = sigma2_1[kx_indices, ky_indices, *]
    sigma2_2 = sigma2_2[kx_indices, ky_indices, *]
    sigma2_diff = sigma2_2 / (sigma2_1 + sigma2_2)
    
    ;PRINT, 'MEAN:' + STRING(MEAN(sigma2_diff, /NAN))
    ;PRINT, 'CORRECTION FACTOR (MEAN):' + STRING(var_prop_correct(1 - MEAN(sigma2_diff, /nan), MEAN(sigma2_diff, /nan)))
    ;PRINT, 'STANDARD DEV:' + STRING(STDDEV(sigma2_diff, /NAN))
    
    binsize = 1/200.
    histvals = HISTOGRAM(sigma2_diff, BINSIZE = binsize, MIN = 0, MAX = 1, /NAN)
    histvals = histvals / TOTAL(histvals) * 100
    hist_x = FINDGEN(N_ELEMENTS(histvals)) * binsize
    histvals_scaled = histvals / MAX(histvals) * .8 * (ymax - 1) + 1
    
    CGOPLOT, hist_x, histvals_scaled, PSYM = 10, COLOR = 'blue'
    
    new_legend_title = 'Variance Distribution, kperp < ' + STRTRIM(STRING(kperp_wavelength_max), 2) + ' wavelengths'
    legend_title = [legend_title, new_legend_title]
    new_legend_color = ['blue']
    legend_color = [legend_color, new_legend_color]
    
  ENDIF
  
  CGLEGEND, TITLE = legend_title, COLOR = legend_color, CHARSIZE = 1, LOCATION = [.57, .85], $
    ALIGNMENT = 4
    
  CGPS_CLOSE, /PNG, /DELETE_PS
  
END