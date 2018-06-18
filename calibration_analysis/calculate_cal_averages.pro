pro calculate_cal_averages

  ;cal_file = '/Users/rubybyrne/array_simulation/fhd_rlb_hex_array_Barry_effect_sim_May2018/calibration/hex_array_sim_cal.sav'
  ;output_name = 'hex_array_cal_amp_errors'
  
  cal_file = '/Users/rubybyrne/array_simulation/fhd_rlb_random_array_Barry_effect_sim_May2018/calibration/random1_array_sim_cal.sav'
  output_name = 'random1_array_cal_amp_errors'
  
  output_path = '/Users/rubybyrne/array_simulation/plots'
  
  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=0.)
  gains[*,*,0] = (*cal.gain[0]) + (*cal.gain_residual[0])
  gains[*,*,1] = (*cal.gain[1]) + (*cal.gain_residual[1])
  gains_amp = abs(gains)
  mean_amp = mean(gains_amp, dimension=3, /nan)  ;average over polarizations
  mean_amp = mean(mean_amp, dimension=2, /nan)  ;average over tiles
  mean_amp_dev = mean_amp-1.
  print, stddev(mean_amp_dev)
  freq = cal.freq/1.e6
  cgWindow_SetDefs, IM_WIDTH=600
  cgplot, freq, mean_amp_dev, PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Frac. Error in Average Gain Amp.', xrange=[min(freq)-2., max(freq)+2.], yrange=[-.016,-.003], $
    output=output_path+'/'+output_name+'.png', charsize=1.
    
  ;fourier, y_time=func_u, time=u_vals, y_freq=mean_amp_dev-mean(mean_amp_dev), frequencies=freq*1.e6, /inverse
    
  ;fourier, func_x=mean_amp_dev-mean(mean_amp_dev), x_vals=freq*1.e6, func_u=func_u, u_vals=u_vals
  ;power_dist = abs(func_u)^2. / total(abs(func_u)^2.) 
  ;cgplot, u_vals*3e8, power_dist, PSym=-16, color='blue',symsize=.8, aspect=.8, xtitle='Calibration Light Travel Dist. (m)', $
  ;  ytitle='Average Gain Amp. FT', xrange=[0, 4000], yrange=[0,max(power_dist)], $
  ;  output=output_path+'/'+output_name+'_ft.png', charsize=1.
    
end