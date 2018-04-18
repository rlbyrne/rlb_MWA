pro plot_bl_dist

  ;cal_file = '/home/rlbyrne/cal_simulation/fhd_rlb_phaseI_Barry_effect_sim_Mar2018/calibration/1061316296_cal.sav'
  ;cal_file = '/home/rlbyrne/cal_simulation/fhd_rlb_phaseII_Barry_effect_sim_Mar2018/calibration/1163765304_cal.sav'

  cal_file = '/home/rlbyrne/cal_simulation/fhd_rlb_phaseI_Barry_effect_sim_max_bl_100_Apr2018/calibration/1061316296_cal.sav'
  ;cal_file = '/home/rlbyrne/cal_simulation/fhd_rlb_phaseII_Barry_effect_sim_max_bl_100_Apr2018/calibration/1163765304_cal.sav'


  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  u_coords = cal.uu
  v_coords = cal.vv
  bl_lengths_m = (cal.uu^2. + cal.vv^2.)^0.5 * 3.e8
  
  cgps_open, '/home/rlbyrne/cal_simulation/bl_density_phaseI.png'
  ;cgWindow_SetDefs, IM_WIDTH=600
  cghistoplot, bl_lengths_m, binsize=5., xtitle='Baseline Length (m)', ytitle='Number of Baselines', $
    maxinput=cal.max_cal_baseline, mininput=cal.min_cal_baseline, xrange=[0.,100.], /fill, charsize=1.
  cgps_close, /png, /delete_ps

end
