PRO cross_diff_var_ratio_plot, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, savefile = savefile, $
    polarization = polarization, use_cubes = use_cubes, use_var_correct = use_var_correct
    
  IF N_ELEMENTS(choose_term) EQ 0 THEN choose_term = 1
  IF N_ELEMENTS(kperp_wavelength_max) LT 1 THEN kperp_wavelength_max = 0
  IF N_ELEMENTS(kperp_wavelength_min) LT 1 THEN kperp_wavelength_min = 0
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  plot_max = 0.0005
  
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
    
  CROSS_DIFFERENCE, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, $
    use_cubes = use_cubes, kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    diff_cross = diff_cross, sigma2_correct = sigma2, polarization = polarization
    
  var_diff_cross = VARIANCE(REAL_PART(diff_cross), DIMENSION = 3)
  sigma2_kperp = MEAN(sigma2, DIMENSION = 3, /NAN)
  
  ratio_var_diff_cross = var_diff_cross / sigma2_kperp
  wh_sig0 = WHERE(sigma2_kperp EQ 0, count_sig0)
  IF count_sig0 GT 0 THEN ratio_var_diff_cross[wh_sig0] = 0
  
  stop
  
  IF KEYWORD_SET(savefile) THEN BEGIN
    filename = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis_plots/varratio_' + cube_name + '_' + polarization + '_term' + STRTRIM(STRING(choose_term),2)
    note = note_part + ', ' + polarization +', term ' + STRTRIM(STRING(choose_term), 2) ;+ ': Mean = ' + STRMID(STRTRIM(STRING(MEAN(ratio_var_diff_cross, /NAN)), 2), 0, 4)
    IF kperp_wavelength_max NE 0 THEN filename = filename + '_kperpmax' + STRTRIM(STRING(kperp_wavelength_max),2)
    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,plot_max], $
      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky', $
      NOTE = note, SAVEFILE = filename, /PNG
  ENDIF ELSE BEGIN
    QUICK_IMAGE, ratio_var_diff_cross, kx_mpc, ky_mpc, WINDOW = 1, DATA_RANGE = [0,plot_max], $
      TITLE = 'Measured / Calculated Variance, Crossed Simulated Noise Cubes', XTITLE = 'kx', YTITLE = 'ky', $
      NOTE = note
  ENDELSE
  
  
  stop
  
END