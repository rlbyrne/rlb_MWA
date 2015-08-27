PRO kz_var_hist, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, Aug23 = Aug23, Aug27 = Aug27

  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/kz_var_hist_term' + STRTRIM(STRING(choose_term), 2), /FONT, XSIZE = 10, YSIZE = 7
  colors = ['black', 'blue', 'dark green']
  legend_colors = []
  legend_titles = []
  
  IF N_ELEMENTS(kperp_wavelength_max) EQ 0 THEN kperp_wavelength_max = 0
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
  
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      Aug23_sigma2 = Aug23_sigma2, Aug27_sigma2 = Aug27_sigma2, $
      kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc
      
    ratio_Aug23_sigma2 = Aug23_sigma2
    ratio_Aug27_sigma2 = Aug27_sigma2
    
    FOR j = 0, (SIZE(Aug23_sigma2))[3] - 1 DO BEGIN
    
      Aug23_ratio_slice = Aug23_sigma2[*,*,j] / MEAN(Aug23_sigma2, DIMENSION = 3)
      wh_sig0 = WHERE(MEAN(Aug23_sigma2, DIMENSION = 3) EQ 0, count_sig0)
      IF count_sig0 GT 0 THEN Aug23_ratio_slice[wh_sig0] = 0
      ratio_Aug23_sigma2[*,*,j] = Aug23_ratio_slice
      
      Aug27_ratio_slice = Aug27_sigma2[*,*,j] / MEAN(Aug27_sigma2, DIMENSION = 3)
      wh_sig0 = WHERE(MEAN(Aug27_sigma2, DIMENSION = 3) EQ 0, count_sig0)
      IF count_sig0 GT 0 THEN Aug27_ratio_slice[wh_sig0] = 0
      ratio_Aug27_sigma2[*,*,j] = Aug27_ratio_slice
      
    ENDFOR
    
    print, MEAN(ratio_Aug23_sigma2, /NAN)
    print, MEAN(ratio_Aug27_sigma2, /NAN)
    
    IF KEYWORD_SET(Aug23) THEN hist_input = ratio_Aug23_sigma2 ELSE BEGIN
      IF KEYWORD_SET(Aug27) THEN hist_input = ratio_Aug27_sigma2 ELSE hist_input = [ratio_Aug23_sigma2, ratio_Aug27_sigma2]
    ENDELSE
    
    IF i EQ 0 THEN do_overplot = 0 ELSE do_overplot = 1
    RLB_HISTPLOT, hist_input, XSTYLE = 1, YSTYLE = 1,  OVERPLOT = do_overplot, RANGE = [.8,1.2], BINSIZE = binsize, $
      YRANGE = [0, 45], COLOR = colors[i], XTITLE = 'Variance / Mean Variance in kz', $
      YTITLE = 'Histogram Count (%)'
      
    legend_colors = [legend_colors, colors[i]]
    IF kperp_wavelength_max[i] EQ 0 THEN new_legend_title = 'all kperp' ELSE new_legend_title = 'kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    legend_titles = [legend_titles, new_legend_title]
    
  ENDFOR
  
  CGLEGEND, TITLE = legend_titles, COLOR = legend_colors, CHARSIZE = 1, LOCATION = [.8, .85], $
    ALIGNMENT = 4
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
  
END