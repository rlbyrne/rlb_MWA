;; old and pointless

PRO cross_diff_var1sim, kperp_wavelength_max = kperp_wavelength_max, choose_term = choose_term, $
    ;;OUTPUTS:
    kx_mpc = kx_mpc, ky_mpc = ky_mpc, kz_mpc = kz_mpc, kperp_lambda_conv = kperp_lambda_conv, $
    diff_cross = diff_cross, sigma2_correct = sigma2_correct
  
  ;; get constants from the Aug23 run:
  kperp_lambda_conv = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
    'KPERP_LAMBDA_CONV')
  kx_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
    'KX_MPC')
  ky_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
    'KY_MPC')
  kz_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
    'KZ_MPC')
    
    
  IF SIZE(kperp_wavelength_max, /N_ELEMENTS) GT 0 THEN BEGIN
    IF kperp_wavelength_max NE 0 THEN BEGIN
      kx_indices = WHERE(ABS(kx_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = WHERE(ABS(ky_mpc) LE kperp_wavelength_max/kperp_lambda_conv, /NULL)
      ky_indices = ky_indices[1:*] ;remove first bin with zero values
      kx_mpc = kx_mpc[kx_indices]
      ky_mpc = ky_mpc[ky_indices]
    ENDIF
  ENDIF
  
  ;;simulate difference cubes
  Aug23_size = SIZE(Aug23_sigma2)
  Aug23_sim_real = RANDOMN(seed, N_ELEMENTS(kx_mpc), N_ELEMENTS(ky_mpc), N_ELEMENTS(kz_mpc))
  Aug23_sim_imag = RANDOMN(seed, N_ELEMENTS(kx_mpc), N_ELEMENTS(ky_mpc), N_ELEMENTS(kz_mpc))
  Aug23_sim = COMPLEX(Aug23_sim_real, Aug23_sim_imag)
  Aug27_size = SIZE(Aug27_sigma2)
  Aug27_sim_real = RANDOMN(seed, N_ELEMENTS(kx_mpc), N_ELEMENTS(ky_mpc), N_ELEMENTS(kz_mpc))
  Aug27_sim_imag = RANDOMN(seed, N_ELEMENTS(kx_mpc), N_ELEMENTS(ky_mpc), N_ELEMENTS(kz_mpc))
  Aug27_sim = COMPLEX(Aug27_sim_real, Aug27_sim_imag)
  
  ;; calculate cross difference
  diff_cross = Aug23_sim * CONJ(Aug27_sim)
  
END