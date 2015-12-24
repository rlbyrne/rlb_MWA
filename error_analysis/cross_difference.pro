;; restores cubes and constants; returns the product of the difference cubes; saves difference cube products

PRO cross_difference, choose_term = choose_term, polarization = polarization, kperp_wavelength_max = kperp_wavelength_max, $
    kperp_wavelength_min = kperp_wavelength_min, uv_img_clip,  $
    obsid1 = obsid1, obsid2 = obsid2, filename1 = filename1, filename2 = filename2, sample_factor = sample_factor, $
    ;;PREFERENCES:
    diff_cross = diff_cross, $           ;calculate cross of difference cubes, defaults 1
    diff_sim = diff_sim, $               ;calculate cross of simulated difference cubes, defaults 0
    calculate_var = calculate_var, $     ;calculate variance of crossed difference cubes, defaults 1
    use_var_correct = use_var_correct, $ ;corrects the calculated variance for nonequal input variances, defaults 1
    normalize = normalize, $             ;normalized variances to 1, defaults 0
    recalc = recalc, $                   ;does not use saved values, defaults 0
    uvf_input = uvf_input, $             ;uses gridded UVF cubes, defaults 0 (Healpix outputs)
    ;;OUTPUTS:
    kperp_lambda_conv=kperp_lambda_conv, delay_params=delay_params, kx_mpc=kx_mpc, ky_mpc=ky_mpc, kz_mpc=kz_mpc, $
    ACube_sigma2=ACube_sigma2, BCube_sigma2=BCube_sigma2, ACube=ACube, BCube=BCube, cube_cross=cube_cross, $
    ACube_sim=ACube_sim, BCube_sim=BCube_sim, sim_cube_cross=sim_cube_cross, approx_sigma2=approx_sigma2, sigma2=sigma2
    
  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  IF N_ELEMENTS(diff_cross) LT 1 then diff_cross = 1
  IF N_ELEMENTS(calculate_var) LT 1 THEN calculate_var = 1
  IF N_ELEMENTS(use_var_correct) LT 1 THEN use_var_correct = 1
  
  IF N_ELEMENTS(obsid1) LT 1 THEN obsid1 = ''
  IF N_ELEMENTS(obsid2) LT 1 THEN obsid2 = ''
  
  save_loc = '/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/'
  if keyword_set(uvf_input) then begin
    if keyword_set(uvf_img_clip) then save_filename = save_loc + 'cross_diff_noise_sim_flatUV' + sample_factor + '_UVFclip' + num_formatter_filename(uvf_img_clip) + '_' + polarization + '_term' + number_formatter(choose_term) + '.sav' $
    else save_filename = save_loc + 'cross_diff_noise_sim_flatUV' + sample_factor + '_' + polarization + '_term' + number_formatter(choose_term) + '.sav'
  endif else begin
    if keyword_set(uvf_img_clip) then save_filename = save_loc + 'cross_diff_noise_sim_Heal_flatUV' + sample_factor + '_UVFclip' + num_formatter_filename(uvf_img_clip) + '_' + polarization + '_term' + number_formatter(choose_term) + '.sav' $
    else save_filename = save_loc + 'cross_diff_noise_sim_Heal_flatUV' + sample_factor + '_' + polarization + '_term' + number_formatter(choose_term) + '.sav'
  endelse
  if ~file_test(save_filename) then recalc = 1
  
  ;; get data and variances:
  if ~keyword_set(recalc) then begin
    kperp_lambda_conv = getvar_savefile(save_filename, 'KPERP_LAMBDA_CONV')
    delay_params = getvar_savefile(save_filename, 'DELAY_PARAMS')
    kx_mpc = getvar_savefile(save_filename, 'KX_MPC')
    ky_mpc = getvar_savefile(save_filename, 'KY_MPC')
    kz_mpc = getvar_savefile(save_filename, 'KZ_MPC')
    ACube_sigma2 = getvar_savefile(save_filename, 'ACube_sigma2')
    BCube_sigma2 = getvar_savefile(save_filename, 'BCube_sigma2')
  endif else begin
    kperp_lambda_conv = getvar_savefile(filename1, 'KPERP_LAMBDA_CONV')
    delay_params = getvar_savefile(filename1, 'DELAY_PARAMS')
    kx_mpc = getvar_savefile(filename1, 'KX_MPC')
    ky_mpc = getvar_savefile(filename1, 'KY_MPC')
    kz_mpc = getvar_savefile(filename1, 'KZ_MPC')
    ACube_sigma2 = getvar_savefile(filename1, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
    BCube_sigma2 = getvar_savefile(filename2, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  endelse
  
  IF n_elements(kperp_wavelength_max) GT 0 THEN BEGIN
    IF kperp_wavelength_max NE 0 THEN BEGIN
      kx_indices_max = WHERE(ABS(kx_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices_max = WHERE(ABS(ky_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices_max = ky_indices_max[1:*] ;remove first bin with zero values
      kx_mpc = kx_mpc[kx_indices_max]
      ky_mpc = ky_mpc[ky_indices_max]
      ACube_sigma2 = ACube_sigma2[kx_indices_max, ky_indices_max, *]
      BCube_sigma2 = BCube_sigma2[kx_indices_max, ky_indices_max, *]
    ENDIF
  ENDIF
  IF n_elements(kperp_wavelength_min) GT 0 THEN BEGIN
    kx_indices_min = WHERE(ABS(kx_mpc) GT kperp_wavelength_min/kperp_lambda_conv, /NULL)
    ky_indices_min = WHERE(ABS(ky_mpc) GT kperp_wavelength_min/kperp_lambda_conv, /NULL)
    kx_mpc = kx_mpc[kx_indices_min]
    ky_mpc = ky_mpc[ky_indices_min]
    ACube_sigma2 = ACube_sigma2[kx_indices_min, ky_indices_min, *]
    BCube_sigma2 = BCube_sigma2[kx_indices_min, ky_indices_min, *]
  ENDIF
  
  if keyword_set(diff_cross) then begin
    if ~keyword_set(recalc) then begin
      ACube = getvar_savefile(save_filename, 'ACube')
      BCube = getvar_savefile(save_filename, 'BCube')
      cube_cross = getvar_savefile(save_filename, 'cube_cross')
    endif else begin
      ACube = getvar_savefile(filename1, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
      BCube = getvar_savefile(filename2, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
      if n_elements(kx_indices_min) gt 0 and n_elements(ky_indices_min) gt 0 then begin
        ACube = ACube[kx_indices_min, ky_indices_min, *]
        BCube = ACube[kx_indices_min, ky_indices_min, *]
      endif
      if n_elements(kx_indices_max) gt 0 and n_elements(ky_indices_max) gt 0 then begin
        ACube = ACube[kx_indices_max, ky_indices_max, *]
        BCube = ACube[kx_indices_max, ky_indices_max, *]
      endif
      IF KEYWORD_SET(normalize) THEN BEGIN
        ACube = ACube / SQRT(ACube_sigma2)
        wh_sig0 = WHERE(ACube_sigma2 EQ 0, count_sig0)
        IF count_sig0 GT 0 THEN ACube[wh_sig0] = 0
        BCube = BCube / SQRT(BCube_sigma2)
        wh_sig0 = WHERE(BCube_sigma2 EQ 0, count_sig0)
        IF count_sig0 GT 0 THEN BCube[wh_sig0] = 0
      ENDIF
      cube_cross = ACube * CONJ(BCube)
    endelse
  endif else begin
    ACube = 0
    BCube = 0
    cube_cross = 0
  endelse
  
  IF KEYWORD_SET(diff_sim) THEN BEGIN
    if ~keyword_set(recalc) then begin
      ACube_sim = getvar_savefile(save_filename, 'ACube_sim')
      BCube_sim = getvar_savefile(save_filename, 'BCube_sim')
      sim_cube_cross = getvar_savefile(save_filename, 'sim_cube_cross')
    endif else begin
      cube_size = SIZE(ACube_sigma2)
      IF KEYWORD_SET(normalize) then begin
        ACube_sigma2_use = make_array(cube_size[1], cube_size[2], cube_size[3], value = 1.)
        BCube_sigma2_use = ACube_sigma2_use
      endif else begin
        ACube_sigma2_use = ACube_sigma2
        BCube_sigma2_use = BCube_sigma2
      endelse
      ACube_sim_real = RANDOMN(seed, cube_size[1], cube_size[2], cube_size[3]) * SQRT(ACube_sigma2_use)
      ACube_sim_imag = RANDOMN(seed, cube_size[1], cube_size[2], cube_size[3]) * SQRT(ACube_sigma2_use)
      ACube_sim = COMPLEX(ACube_sim_real, ACube_sim_imag)
      BCube_sim_real = RANDOMN(seed, cube_size[1], cube_size[2], cube_size[3]) * SQRT(BCube_sigma2_use)
      BCube_sim_imag = RANDOMN(seed, cube_size[1], cube_size[2], cube_size[3]) * SQRT(BCube_sigma2_use)
      BCube_sim = COMPLEX(BCube_sim_real, BCube_sim_imag)
      sim_cube_cross = ACube_sim * CONJ(BCube_sim)
    endelse
  ENDIF else begin
    ACube_sim = 0
    BCube_sim = 0
    sim_cube_cross = 0
  endelse
  
  IF calculate_var EQ 1 THEN BEGIN
    if ~keyword_set(recalc) then sigma2 = getvar_savefile(save_filename, 'sigma2') else begin
      approx_sigma2 = (ACube_sigma2 + BCube_sigma2)^2/2
      IF use_var_correct EQ 1 THEN BEGIN
        correction_factor = VAR_PROP_CORRECT(ACube_sigma2, BCube_sigma2)
        sigma2 = approx_sigma2 / correction_factor
        wh_sig0 = WHERE(correction_factor EQ 0 OR ~FINITE(correction_factor), count_sig0)
        IF count_sig0 GT 0 THEN sigma2[wh_sig0] = 0
      ENDIF ELSE sigma2 = approx_sigma2
    endelse
  ENDIF
  
  if keyword_set(recalc) then begin
    print, 'Saving crossed difference cubes here: ' + save_filename
    file_mkdir, save_loc
    save, filename = save_filename, kperp_lambda_conv, delay_params, kx_mpc, ky_mpc, kz_mpc, ACube_sigma2, BCube_sigma2, $
      ACube, BCube, cube_cross, ACube_sim, BCube_sim, sim_cube_cross, sigma2
  endif
  
END