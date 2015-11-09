;; restores cubes and constants; returns the product of the difference cubes

; use_cubes values:
; 0: Long Run, Aug23 and Aug27 zenith pointings (integrated 30 min) -- default behavior
; 1: Simulated noise cubes based on Long Run obsids 1061316176 and 1061316296 (Aug23 near zenith) (single obsid) -- also called with keyword noise_sim
; 2: Long Run, Aug23 and Aug27 full 3 hour integration -- also called with keyword three_hr
; 3: Long Run, obsids 1061316176 and 1061316296 (Aug23 near zenith) (single obsid)
; 4: Simulated noise cubes with uniform 0.0002 UV coverage based on Long Run obsids 1061316176 and 1061316296 (Aug23 near zenith) (single obsid)
; 5: Simulated noise cubes with uniform 0.0005 UV coverage based on Long Run obsids 1061316176 and 1061316296 (Aug23 near zenith) (single obsid)
; diff_sim: produces simulated difference cubes based on the variances determined by use_cubes

PRO cross_difference, use_cubes = use_cubes, choose_term = choose_term, polarization = polarization, kperp_wavelength_max = kperp_wavelength_max, $
    kperp_wavelength_min = kperp_wavelength_min, normalize = normalize, diff_sim = diff_sim, use_var_correct = use_var_correct, calculate_var = calculate_var, $
    ;;OUTPUTS:
    ACube = ACube, ACube_sigma2 = ACube_sigma2, BCube = BCube, BCube_sigma2 = BCube_sigma2, $
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params,$
    cube_cross = cube_cross, sigma2_approx = sigma2_approx, sigma2_correct = sigma2_correct, $
    obsid1 = obsid1, obsid2 = obsid2, filename1 = filename1, filename2 = filename2, $
    ;;FOR HISTORICAL REASONS:
    Aug23_diff = Aug23_diff, Aug23_sigma2 = Aug23_sigma2, Aug27_diff = Aug27_diff, Aug27_sigma2 = Aug27_sigma2, diff_cross = diff_cross, $
    noise_sim = noise_sim, three_hr = three_hr
    
  IF N_ELEMENTS(choose_term) LT 1 THEN choose_term = 1
  IF N_ELEMENTS(use_cubes) LT 1 THEN use_cubes = 0
  IF KEYWORD_SET(noise_sim) THEN use_cubes = 1
  IF KEYWORD_SET(three_hr) THEN use_cubes = 2
  IF N_ELEMENTS(polarization) LT 1 THEN polarization = 'xx'
  IF N_ELEMENTS(use_var_correct) LT 1 THEN use_var_correct = 1
  IF N_ELEMENTS(calculate_var) LT 1 THEN calculate_var = 1
  IF polarization NE 'xx' AND polarization NE 'yy' THEN polarization = 'xx'
  
  CASE use_cubes OF
    0: BEGIN ;Long Run Aug23 and Aug27 zenith pointings
      polarization = 'xx'
      filename1 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave'
      filename2 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456532_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave'
    END
    1: BEGIN ;Simulated noise cubes based on Long Run obsids 1061316176 and 1061316296
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/'+ obsid1 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/'+ obsid2 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
    END
    2: BEGIN ;Long Run, Aug23 and Aug27 full 3 hour integration
      filename1 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_rlb_aug23_3hr_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_rlb_aug27_3hr_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
    END
    3: BEGIN ;Long Run obsids 1061316176 and 1061316296
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/' + obsid1 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/' + obsid2 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
    END
    4: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0002
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid1 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid2 + '_cubeXX__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    5: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0005
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid1 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid2 + '_cubeXX__even_odd_joint_res_' + polarization + '_bh_kcube.idlsave'
    END
    6: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0002, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0002/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    7: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.0005, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.0005/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    8: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.001, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.001/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.001/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    9: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.005, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.005/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.005/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    10: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.01, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.01/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/nfs/eor-00/h1/rbyrne/MWA/sim_cubes/flatUV_0.01/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    11: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.05, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_0.05/ps/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_0.05/ps/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END 
    12: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.1, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_0.1/ps/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_0.1/ps/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    13: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 0.5, uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_0.5/ps/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_0.5/ps/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
    14: BEGIN ;Simulated flat UV coverage noise cubes based on Long Run obsids 1061316176 and 1061316296, UV coverage 1., uvf_input keyword set
      obsid1 = '1061316176'
      obsid2 = '1061316296'
      filename1 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid1 + '_1/ps/'+ obsid1 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
      filename2 = '/data3/MWA/FHD_Aug23/fhd_rlb_noise_sim_flatUV_'+ obsid2 + '_1/ps/'+ obsid2 + '_gridded_uvf__even_odd_joint_model_' + polarization + '_bh_kcube.idlsave'
    END
  ENDCASE
  
  IF N_ELEMENTS(obsid1) LT 1 THEN obsid1 = ''
  IF N_ELEMENTS(obsid2) LT 1 THEN obsid2 = ''
  
  ;; get data and variances:
  
  ACube_sigma2 = getvar_savefile(filename1, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  BCube_sigma2 = getvar_savefile(filename2, 'SIGMA2_' + STRTRIM(STRING(choose_term),2))
  IF KEYWORD_SET(diff_sim) THEN BEGIN ;diff_sim keyword is an alternative to calling cross_diff_sim
    ACube_size = SIZE(ACube_sigma2)
    ACube_real = RANDOMN(seed, ACube_size[1], ACube_size[2], ACube_size[3]) * SQRT(ACube_sigma2)
    ACube_imag = RANDOMN(seed, ACube_size[1], ACube_size[2], ACube_size[3]) * SQRT(ACube_sigma2)
    ACube = COMPLEX(ACube_real, ACube_imag)
    BCube_size = SIZE(BCube_sigma2)
    BCube_real = RANDOMN(seed, BCube_size[1], BCube_size[2], BCube_size[3]) * SQRT(BCube_sigma2)
    BCube_imag = RANDOMN(seed, BCube_size[1], BCube_size[2], BCube_size[3]) * SQRT(BCube_sigma2)
    BCube = COMPLEX(BCube_real, BCube_imag)
  ENDIF ELSE BEGIN
    ACube = getvar_savefile(filename1, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
    BCube = getvar_savefile(filename2, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  ENDELSE
  
  ;; get constants:
  kperp_lambda_conv = getvar_savefile(filename1, 'KPERP_LAMBDA_CONV')
  delay_params = getvar_savefile(filename1, 'DELAY_PARAMS')
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  kz_mpc = getvar_savefile(filename1, 'KZ_MPC')
  
  IF n_elements(kperp_wavelength_max) GT 0 THEN BEGIN
    IF kperp_wavelength_max NE 0 THEN BEGIN
      kx_indices = WHERE(ABS(kx_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = WHERE(ABS(ky_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = ky_indices[1:*] ;remove first bin with zero values
      kx_mpc = kx_mpc[kx_indices]
      ky_mpc = ky_mpc[ky_indices]
      ACube = ACube[kx_indices, ky_indices, *]
      ACube_sigma2 = ACube_sigma2[kx_indices, ky_indices, *]
      BCube = BCube[kx_indices, ky_indices, *]
      BCube_sigma2 = BCube_sigma2[kx_indices, ky_indices, *]
    ENDIF
  ENDIF
  
  IF n_elements(kperp_wavelength_min) GT 0 THEN BEGIN
    kx_indices = WHERE(ABS(kx_mpc) GT kperp_wavelength_min/kperp_lambda_conv, /NULL)
    ky_indices = WHERE(ABS(ky_mpc) GT kperp_wavelength_min/kperp_lambda_conv, /NULL)
    kx_mpc = kx_mpc[kx_indices]
    ky_mpc = ky_mpc[ky_indices]
    ACube = ACube[kx_indices, ky_indices, *]
    ACube_sigma2 = ACube_sigma2[kx_indices, ky_indices, *]
    BCube = BCube[kx_indices, ky_indices, *]
    BCube_sigma2 = BCube_sigma2[kx_indices, ky_indices, *]
  ENDIF
  
  ;; normalize the variances to 1: (used for cross_diff_var1sim_compare_hist.pro)
  IF KEYWORD_SET(normalize) THEN BEGIN
    ACube = ACube / SQRT(ACube_sigma2)
    wh_sig0 = WHERE(ACube_sigma2 EQ 0, count_sig0)
    IF count_sig0 GT 0 THEN ACube[wh_sig0] = 0
    BCube = BCube / SQRT(BCube_sigma2)
    wh_sig0 = WHERE(BCube_sigma2 EQ 0, count_sig0)
    IF count_sig0 GT 0 THEN BCube[wh_sig0] = 0
  ENDIF
  
  ;; calculate cross difference
  cube_cross = ACube * CONJ(BCube)
  
  IF calculate_var EQ 1 THEN BEGIN
    sigma2_approx = (ACube_sigma2 + BCube_sigma2)^2/2
    IF use_var_correct EQ 1 THEN BEGIN
      correction_factor = VAR_PROP_CORRECT(ACube_sigma2, BCube_sigma2)
      sigma2_correct = sigma2_approx / correction_factor
      wh_sig0 = WHERE(correction_factor EQ 0 OR ~FINITE(correction_factor), count_sig0)
      IF count_sig0 GT 0 THEN sigma2_correct[wh_sig0] = 0
    ENDIF ELSE sigma2_correct = sigma2_approx
  ENDIF
  
  ;; for historical reasons (variable name changes):
  Aug23_diff = ACube
  Aug23_sigma2 = ACube_sigma2
  Aug27_diff = BCube
  Aug27_sigma2 = BCube_sigma2
  diff_cross = cube_cross
  
END