;; produces simulated noise cubes from the given variance cubes, assuming Gaussian distributed values; returns the product of the simulated noise cubes

PRO cross_diff_sim, choose_term = choose_term, kperp_wavelength_max = kperp_wavelength_max, kperp_wavelength_min = kperp_wavelength_min, $
    normalize = normalize, noise_sim = noise_sim, polarization = polarization, $
    ;;OUTPUTS:
    ACube = ACube, ACube_sigma2 = ACube_sigma2, BCube = BCube, BCube_sigma2 = BCube_sigma2, $
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params,$
    cube_cross = cube_cross, sigma2_approx = sigma2_approx, sigma2_correct = sigma2_correct, $
    obsid1 = obsid1, obsid2 = obsid2, $
    ;;FOR HISTORICAL REASONS:
    Aug23_diff = Aug23_diff, Aug23_sigma2 = Aug23_sigma2, Aug27_diff = Aug27_diff, Aug27_sigma2 = Aug27_sigma2, diff_cross_sim = diff_cross_sim
    
  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1
  IF KEYWORD_SET(noise_sim) THEN BEGIN ;using simulated noise cubes
    IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
    IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
    obsid1 = '1061316176'
    obsid2 = '1061316296'
    filename1 = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis/'+ obsid1 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    filename2 = '/nfs/eor-00/h1/rbyrne/MWA/error_analysis/'+ obsid2 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
  ENDIF ELSE BEGIN ;using Long Run Aug23 and Aug27 zenith pointing cubes
    obsid1 = ''
    obsid2 = ''
    filename1 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave'
    filename2 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave'
  ENDELSE
  
  ;; get data and variances:
  ACube = getvar_savefile(filename1, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  ACube_sigma2 = getvar_savefile(filename1, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  BCube = getvar_savefile(filename2, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  BCube_sigma2 = getvar_savefile(filename2, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  
  ;; get constants:
  kperp_lambda_conv = getvar_savefile(filename1, 'KPERP_LAMBDA_CONV')
  delay_params = getvar_savefile(filename1, 'DELAY_PARAMS')
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  kz_mpc = getvar_savefile(filename1, 'KZ_MPC')
  
  
  IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN BEGIN
    IF kperp_wavelength_max NE 0 THEN BEGIN
      kx_indices = WHERE(ABS(kx_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = WHERE(ABS(ky_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = ky_indices[1:*] ;remove first bin with zero values
      kx_mpc = kx_mpc[kx_indices]
      ky_mpc = ky_mpc[ky_indices]
      ACube_sigma2 = ACube_sigma2[kx_indices, ky_indices, *]
      BCube_sigma2 = BCube_sigma2[kx_indices, ky_indices, *]
    ENDIF
  ENDIF
  
  IF SIZE(kperp_wavelength_min, /N_ELEMENTS) GT 0 THEN BEGIN
    kx_indices = WHERE(ABS(kx_mpc) GT kperp_wavelength_max/kperp_lambda_conv, /NULL)
    ky_indices = WHERE(ABS(ky_mpc) GT kperp_wavelength_max/kperp_lambda_conv, /NULL)
    kx_mpc = kx_mpc[kx_indices]
    ky_mpc = ky_mpc[ky_indices]
    ACube = ACube[kx_indices, ky_indices, *]
    ACube_sigma2 = ACube_sigma2[kx_indices, ky_indices, *]
    BCube = BCube[kx_indices, ky_indices, *]
    BCube_sigma2 = BCube_sigma2[kx_indices, ky_indices, *]
  ENDIF
  
  ;;simulate difference cubes
  ACube_size = SIZE(ACube_sigma2)
  ACube_sim_real = RANDOMN(seed, ACube_size[1], ACube_size[2], ACube_size[3]) * SQRT(ACube_sigma2)
  ACube_sim_imag = RANDOMN(seed, ACube_size[1], ACube_size[2], ACube_size[3]) * SQRT(ACube_sigma2)
  ACube_sim = COMPLEX(ACube_sim_real, ACube_sim_imag)
  BCube_size = SIZE(BCube_sigma2)
  BCube_sim_real = RANDOMN(seed, BCube_size[1], BCube_size[2], BCube_size[3]) * SQRT(BCube_sigma2)
  BCube_sim_imag = RANDOMN(seed, BCube_size[1], BCube_size[2], BCube_size[3]) * SQRT(BCube_sigma2)
  BCube_sim = COMPLEX(BCube_sim_real, BCube_sim_imag)
  
  ;; calculate cross difference
  cube_cross = ACube_sim * CONJ(BCube_sim)
  sigma2_approx = (ACube_sigma2 + BCube_sigma2)^2/2
  sigma2_correct = sigma2_approx / VAR_PROP_CORRECT(ACube_sigma2, BCube_sigma2)
  
  ;; for historical reasons (variable name changes):
  Aug23_diff = ACube
  Aug23_sigma2 = ACube_sigma2
  Aug27_diff = BCube
  Aug27_sigma2 = BCube_sigma2
  diff_cross_sim = cube_cross
  
END