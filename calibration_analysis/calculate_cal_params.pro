pro calculate_cal_params

  run_path = '/Users/rubybyrne/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018'
  obsids = ['hex_array_sim_10m','hex_array_sim_15m','hex_array_sim_20m', $
    'random1_array_sim_10m','random1_array_sim_15m','random1_array_sim_20m', $
    'random2_array_sim_10m','random2_array_sim_15m','random2_array_sim_20m', $
    'random3_array_sim_10m','random3_array_sim_15m','random3_array_sim_20m']
  output_path = '/Users/rubybyrne/array_simulation/cal_params_plots'
  
  for obs_index = 0,n_elements(obsids)-1 do begin
    
    obsid = obsids[obs_index]
    cal_file = run_path+'/calibration/'+obsid+'_cal.sav'
    obs_file = run_path+'/metadata/'+obsid+'_obs.sav'
  
    cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
    
    gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=0.)
    gains[*,*,0] = (*cal.gain[0]) + (*cal.gain_residual[0])
    gains[*,*,1] = (*cal.gain[1]) + (*cal.gain_residual[1])
    gains_amp = abs(gains)
    mean_amp = mean(gains_amp, dimension=3, /nan)  ;average over polarizations
    mean_amp = mean(mean_amp, dimension=2, /nan)  ;average over tiles
    mean_amp_dev = mean_amp-1.
    freq = cal.freq/1.e6
    cgWindow_SetDefs, IM_WIDTH=600
    cgplot, freq, mean_amp_dev, PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Frac. Error in Average Gain Amp.', xrange=[min(freq)-2., max(freq)+2.], yrange=[-.022,0.01], $
      output=output_path+'/'+obsid+'_cal_amp_errors.png', charsize=1.
  
    ;fourier, y_time=func_u, time=u_vals, y_freq=mean_amp_dev-mean(mean_amp_dev), frequencies=freq*1.e6, /inverse
  
    ;fourier, func_x=mean_amp_dev-mean(mean_amp_dev), x_vals=freq*1.e6, func_u=func_u, u_vals=u_vals
    ;power_dist = abs(func_u)^2. / total(abs(func_u)^2.)
    ;cgplot, u_vals*3e8, power_dist, PSym=-16, color='blue',symsize=.8, aspect=.8, xtitle='Calibration Light Travel Dist. (m)', $
    ;  ytitle='Average Gain Amp. FT', xrange=[0, 4000], yrange=[0,max(power_dist)], $
    ;  output=output_path+'/'+output_name+'_ft.png', charsize=1.
  
    gains_phase = atan(gains,/phase)
  
    obs_struct = getvar_savefile(obs_file, 'obs', /compatibility_mode)
    ant_info = *obs_struct.meta_data
    ant_pos = make_array(128, 2, /float, value=0.)
    for ant_ind=0, n_elements(ant_info)-1 do begin
      if ant_info[ant_ind].pol eq 'X' then begin  ; use x pol positions only
        ant_pos[ant_info[ant_ind].antenna, 0] = ant_info[ant_ind].north
        ant_pos[ant_info[ant_ind].antenna, 1] = ant_info[ant_ind].east
      endif
    endfor
    weights = make_array(128, 2, /float, value=1.)  ; use equal weighting for now
  
    fit_mat = make_array(3, 3, /float, value=0.)  ; initialize array
    fit_mat[0,0] = total(weights, /nan)
    fit_mat[0,1] = total([weights[*,0]*ant_pos[*,0], weights[*,1]*ant_pos[*,0]], /nan)
    fit_mat[0,2] = total([weights[*,0]*ant_pos[*,1], weights[*,1]*ant_pos[*,1]], /nan)
    fit_mat[1,0] = total([weights[*,0]*ant_pos[*,0], weights[*,1]*ant_pos[*,0]], /nan)
    fit_mat[1,1] = total([weights[*,0]*(ant_pos[*,0]^2), weights[*,1]*(ant_pos[*,0]^2)], /nan)
    fit_mat[1,2] = total([weights[*,0]*ant_pos[*,0]*ant_pos[*,1], weights[*,1]*ant_pos[*,0]*ant_pos[*,1]], /nan)
    fit_mat[2,0] = total([weights[*,0]*ant_pos[*,1], weights[*,1]*ant_pos[*,1]], /nan)
    fit_mat[2,1] = total([weights[*,0]*ant_pos[*,0]*ant_pos[*,1], weights[*,1]*ant_pos[*,0]*ant_pos[*,1]], /nan)
    fit_mat[2,2] = total([weights[*,0]*(ant_pos[*,1]^2), weights[*,1]*(ant_pos[*,1]^2)], /nan)
  
    fit_mat_inv = invert(fit_mat)
  
    fit_params = make_array(384, 3, /float, value=0.)
    for freq_ind=0,383 do begin
      fit_vec = make_array(3, /float, value=0.)
      fit_vec[0] = total(weights[*,*]*reform(gains_phase[freq_ind,*,*]), /nan)
      fit_vec[1] = total([weights[*,0]*reform(gains_phase[freq_ind,*,0])*ant_pos[*,0], weights[*,1]*reform(gains_phase[freq_ind,*,1])*ant_pos[*,0]], /nan)
      fit_vec[2] = total([weights[*,0]*reform(gains_phase[freq_ind,*,0])*ant_pos[*,1], weights[*,1]*reform(gains_phase[freq_ind,*,1])*ant_pos[*,1]], /nan)
  
      fit_params[freq_ind,*] = fit_mat_inv#fit_vec
  
    endfor
  
    freq = cal.freq/1.e6
  
    cgps_open, output_path+'/'+obsid+'_cal_phase_errors.png'
    cgWindow_SetDefs, IM_WIDTH=600
    cgplot, freq, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Average Gain Phase Error', xrange=[min(freq)-2., max(freq)+2.], yrange=[-.015,.025], charsize=1.
    cgps_close, /png, /delete_ps
  
    ;  fourier, func_x=fit_params[*,0]-mean(fit_params[*,0]), x_vals=freq*1.e6, func_u=func_u, u_vals=u_vals
    ;  power_dist = abs(func_u)^2. / total(abs(func_u)^2.)
    ;  cgplot, u_vals*3e8, power_dist, PSym=-16, color='blue',symsize=.8, aspect=.8, xtitle='Calibration Light Travel Dist. (m)', $
    ;    ytitle='Average Gain Phase. FT', xrange=[0,4000], yrange=[0,max(power_dist)], $
    ;    output=output_path+'/'+output_name+'_ft.png', charsize=1.
  
    cgps_open, output_path+'/'+obsid+'_cal_phase_grad_errors.png'
    cgplot, freq, fit_params[*,1], PSym=16, symsize=.5, color='red', aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Average Gain Phase Gradient Error', xrange=[min(freq)-2., max(freq)+2.], yrange=[-5e-6,6e-6], charsize=1.
    cgplot, freq, fit_params[*,2], PSym=16, symsize=.5, color='blue', /overplot
    cglegend, title=['X', 'Y'], psym=16, color=['red','blue'], length=0., /center_sym, location=[.85,.9], charsize=1., /box
    cgps_close, /png, /delete_ps
  
    ;  cgps_open, output_path+'/'+output_name+'_grad_ft.png'
    ;  fourier, func_x=fit_params[*,1]-mean(fit_params[*,1]), x_vals=freq*1.e6, func_u=func_u, u_vals=u_vals
    ;  power_dist = abs(func_u)^2. / total(abs(func_u)^2.)
    ;  cgplot, u_vals*3e8, power_dist, PSym=-16, color='red',symsize=.8, aspect=.8, xtitle='Calibration Light Travel Dist. (m)', $
    ;    ytitle='Average Gain Phase Grad. FT', xrange=[0, 4000], yrange=[0,max(power_dist)*1.1], charsize=1.
    ;  fourier, func_x=fit_params[*,2]-mean(fit_params[*,2]), x_vals=freq*1.e6, func_u=func_u, u_vals=u_vals
    ;  power_dist = abs(func_u)^2. / total(abs(func_u)^2.)
    ;  cgplot, u_vals*3e8, power_dist, PSym=-16, color='blue',symsize=.8, /overplot
    ;  cglegend, title=['X', 'Y'], psym=16, color=['red','blue'], length=0., /center_sym, location=[.85,.9], charsize=1., /box
    ;  cgps_close, /png, /delete_ps
    
  endfor  

end