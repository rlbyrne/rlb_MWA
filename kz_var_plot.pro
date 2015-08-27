PRO kz_var_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, Aug23 = Aug23, Aug27 = Aug27

  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  
  IF KEYWORD_SET(savefile) THEN CGPS_OPEN, '/nfs/eor-00/h1/rbyrne/error_analysis/kz_var_plot_term' + STRTRIM(STRING(choose_term), 2), /FONT, XSIZE = 10, YSIZE = 7
  colors = ['black', 'blue', 'dark green']
  legend_colors = []
  legend_titles = []
  
  IF N_ELEMENTS(kperp_wavelength_max) EQ 0 THEN kperp_wavelength_max = 0
  FOR i = 0, N_ELEMENTS(kperp_wavelength_max) - 1 DO BEGIN
  
    CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max[i], $
      Aug23_sigma2 = Aug23_sigma2, Aug27_sigma2 = Aug27_sigma2, $
      kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc
      
    kz_Aug23_sigma2 = kz_mpc
    kz_Aug27_sigma2 = kz_mpc
    
    FOR j = 0, N_ELEMENTS(kz_mpc) - 1 DO BEGIN
    
      Aug23_ratio_slice = Aug23_sigma2[*,*,j] / MEAN(Aug23_sigma2, DIMENSION = 3)
      wh_sig0 = WHERE(MEAN(Aug23_sigma2, DIMENSION = 3) EQ 0, count_sig0)
      IF count_sig0 GT 0 THEN Aug23_ratio_slice[wh_sig0] = 0
      kz_Aug23_sigma2[j] = MEAN(Aug23_ratio_slice)
      
      Aug27_ratio_slice = Aug27_sigma2[*,*,j] / MEAN(Aug27_sigma2, DIMENSION = 3)
      wh_sig0 = WHERE(MEAN(Aug27_sigma2, DIMENSION = 3) EQ 0, count_sig0)
      IF count_sig0 GT 0 THEN Aug27_ratio_slice[wh_sig0] = 0
      kz_Aug27_sigma2[j] = MEAN(Aug27_ratio_slice)
      
    ENDFOR
    
    IF KEYWORD_SET(Aug23) THEN y_vals = kz_Aug23_sigma2 ELSE BEGIN
      IF KEYWORD_SET(Aug27) THEN y_vals = kz_Aug27_sigma2 ELSE y_vals = (kz_Aug23_sigma2 + kz_Aug27_sigma2) / 2
    ENDELSE
    
    IF i EQ 0 THEN do_overplot = 0 ELSE do_overplot = 1
    CGPLOT, kz_mpc, y_vals, XSTYLE = 1, YSTYLE = 1,  OVERPLOT = do_overplot, COLOR = colors[i], XTITLE = 'kz', $
      YTITLE = 'Mean, Variance / Expected Variance', YRANGE = [0,2]
      
    legend_colors = [legend_colors, colors[i]]
    IF kperp_wavelength_max[i] EQ 0 THEN new_legend_title = 'all kperp' ELSE new_legend_title = 'kperp < ' + STRTRIM(STRING(kperp_wavelength_max[i]), 2)
    legend_titles = [legend_titles, new_legend_title]
    
  ENDFOR
  
  CGLEGEND, TITLE = legend_titles, COLOR = legend_colors, CHARSIZE = 1, LOCATION = [.8, .85], $
    ALIGNMENT = 4
    
  IF KEYWORD_SET(savefile) THEN CGPS_CLOSE, /PNG, /DELETE_PS
  
  
END