;; old version of cross_diff_sim

PRO cross_difference_debug, Aug23_sigma2_const, Aug27_sigma2_const, percent_error, kx_mpc, ky_mpc, kz_mpc

  IF SIZE(kx_mpc, /N_ELEMENTS) EQ 0 THEN BEGIN
    kx_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KX_MPC')
  ENDIF
  IF SIZE(ky_mpc, /N_ELEMENTS) EQ 0 THEN BEGIN
    ky_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KY_MPC')
  ENDIF
  IF SIZE(kz_mpc, /N_ELEMENTS) EQ 0 THEN BEGIN
    kz_mpc = getvar_savefile('/nfs/eor-03/r1/EoR2013/fhd_apb_EoR0_high_sem1_1/ps/Combined_obs_jd2456528_pointing0_wedge_cut_cubeXX__even_odd_joint_res_xx_bh_kcube.idlsave', $
      'KZ_MPC')
  ENDIF
  
  
  Aug23_diff_real = RANDOMN(seed, SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS)) * SQRT(Aug23_sigma2_const)
  Aug23_diff_imag = RANDOMN(seed, SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS)) * SQRT(Aug23_sigma2_const)
  Aug23_diff = COMPLEX(Aug23_diff_real, Aug23_diff_imag)
  Aug23_sigma2 = MAKE_ARRAY(SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS), /FLOAT, VALUE = Aug23_sigma2_const)
  
  Aug27_diff_real = RANDOMN(seed, SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS)) * SQRT(Aug27_sigma2_const)
  Aug27_diff_imag = RANDOMN(seed, SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS)) * SQRT(Aug27_sigma2_const)
  Aug27_diff = COMPLEX(Aug27_diff_real, Aug27_diff_imag)
  Aug27_sigma2 = MAKE_ARRAY(SIZE(kx_mpc, /N_ELEMENTS), SIZE(ky_mpc, /N_ELEMENTS), SIZE(kz_mpc, /N_ELEMENTS), /FLOAT, VALUE = Aug27_sigma2_const)
  
  
  ;; calculate cross difference and weights
  diff_cross = Aug23_diff * CONJ(Aug27_diff)
  sigma2 = (Aug23_sigma2 + Aug27_sigma2)^2/2
  weights = 1/sigma2
  wh_sig0 = WHERE(sigma2 EQ 0, count_sig0)
  IF count_sig0 GT 0 THEN weights[wh_sig0] = 0
  
  variance_measured = stddev(real_part(diff_cross))^2
  variance_calc = (Aug23_sigma2_const+Aug27_sigma2_const)^2 / 2
  percent_error = variance_calc / variance_measured
  
  
;; rebin to 2D
;diff_cross_2d = kspace_rebinning_2d(diff_cross, kx_mpc, ky_mpc, kz_mpc, $
;     kperp_edges_mpc, kpar_edges_mpc, weights = weights, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
;     binned_weights = weights_2d)
  

END