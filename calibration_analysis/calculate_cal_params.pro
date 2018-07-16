pro calculate_cal_params

  run_path = '/Users/rubybyrne/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018'
  obsids = ['hex_array_sim_10m','hex_array_sim_15m','hex_array_sim_20m', $
    'random1_array_sim_10m','random1_array_sim_15m','random1_array_sim_20m', $
    'random2_array_sim_10m','random2_array_sim_15m','random2_array_sim_20m', $
    'random3_array_sim_10m','random3_array_sim_15m','random3_array_sim_20m']
  output_path = '/Users/rubybyrne/array_simulation/cal_params_plots'
  pol = 0  ; polarization, 0=XX, 1=YY
  
  obsids = obsids[0]  ; for debugging
  
  for obs_index = 0,n_elements(obsids)-1 do begin
    
    obsid = obsids[obs_index]
    cal_file = run_path+'/calibration/'+obsid+'_cal.sav'
    obs_file = run_path+'/metadata/'+obsid+'_obs.sav'
    
    if pol eq 0 then pol_name = 'XX' else pol_name='YY'
  
    cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
    
    legend_loc = [.5,.2]
    ft_legend_loc = [.5,.8]
    ft_x_range = [0,2e-6]
    
    ; Plot average gain amplitude
    gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], /complex, value=0.)
    gains[*,*] = (*cal.gain[pol]) + (*cal.gain_residual[pol])
    gains_amp = abs(gains)
    mean_amp = mean(gains_amp, dimension=2, /nan)  ;average over antennas
    mean_amp_dev = mean_amp-1.
    freq = cal.freq/1.e6
    cgWindow_SetDefs, IM_WIDTH=600
    cgplot, freq, mean_amp, PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Average Gain Amp.', title='Average Gain Amplitude Errors', xrange=[min(freq), max(freq)], yrange=[1-.022,1.01], $
      output=output_path+'/'+obsid+'_cal_amp_errors_'+pol_name+'.png', charsize=1.
      
    ; Plot average gain amplitude with per-antenna gain amplitudes
    cgps_open, output_path+'/'+obsid+'_cal_ant_amp_errors_'+pol_name+'.png'
    cgplot, [freq[0], freq[-1]], [1,1], PSym=-0, symsize=.5, color='black', aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Amplitude', title='Overall Gain Amplitude A',$
      xrange=[min(freq), max(freq)], yrange=[1-.022,1.01], charsize=1.
    for antenna = 0, (size(gains_amp))[2]-1 do begin
      cgplot, freq, gains_amp[*, antenna], PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    cgplot, freq, mean_amp, PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['A (Overall Gain Amplitude)', 'Per-Antenna Gain Amplitudes'], psym=-0, color=['blue','gray'], length=0.03, /center_sym, location=legend_loc, charsize=1.
    cgps_close, /png, /delete_ps
    
    ; Plot Fourier Transform of average gain amplitude with per-antenna gain amplitudes
    cgps_open, output_path+'/'+obsid+'_cal_ant_amp_errors_'+pol_name+'_ft.png'
    cgplot, [freq[0], freq[-1]], [1,1], PSym=-0, symsize=.5, color='blue', aspect=.8, xtitle='Calibration Delay (s)', $
      ytitle='Gain Amplitude Power', title='Overall Gain Amplitude A - Power Spectrum',$
      xrange=ft_x_range, yrange=[0,4e8], charsize=1., /nodata
    for antenna = 0, (size(gains_amp))[2]-1 do begin
      fourier, y_time=gains_amp_ft, time=x_axis_ft, y_freq=gains_amp[*, antenna], frequencies=freq*1.e6, /inverse
      cgplot, x_axis_ft, abs(gains_amp_ft)^2., PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    fourier, y_time=mean_amp_ft, time=x_axis_ft, y_freq=mean_amp, frequencies=freq*1.e6, /inverse
    cgplot, x_axis_ft, abs(mean_amp_ft)^2, PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['A (Overall Gain Amplitude)', 'Per-Antenna Gain Amplitudes'], psym=-0, color=['blue','gray'], length=0.03, /center_sym, location=ft_legend_loc, charsize=1.
    cgps_close, /png, /delete_ps
    
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
    fit_mat[0,0] = total(weights[*,pol], /nan)
    fit_mat[0,1] = total(weights[*,pol]*ant_pos[*,0], /nan)
    fit_mat[0,2] = total(weights[*,pol]*ant_pos[*,1], /nan)
    fit_mat[1,0] = total(weights[*,pol]*ant_pos[*,0], /nan)
    fit_mat[1,1] = total(weights[*,pol]*(ant_pos[*,0]^2), /nan)
    fit_mat[1,2] = total(weights[*,pol]*ant_pos[*,0]*ant_pos[*,1], /nan)
    fit_mat[2,0] = total(weights[*,pol]*ant_pos[*,1], /nan)
    fit_mat[2,1] = total(weights[*,pol]*ant_pos[*,0]*ant_pos[*,1], /nan)
    fit_mat[2,2] = total(weights[*,pol]*(ant_pos[*,1]^2), /nan)
  
    fit_mat_inv = invert(fit_mat)
  
    fit_params = make_array(384, 3, /float, value=0.)
    for freq_ind=0,383 do begin
      fit_vec = make_array(3, /float, value=0.)
      fit_vec[0] = total(weights[*,pol]*reform(gains_phase[freq_ind,*]), /nan)
      fit_vec[1] = total(weights[*,pol]*reform(gains_phase[freq_ind,*])*ant_pos[*,0], /nan)
      fit_vec[2] = total(weights[*,pol]*reform(gains_phase[freq_ind,*])*ant_pos[*,1], /nan)
  
      fit_params[freq_ind,*] = fit_mat_inv#fit_vec
  
    endfor
  
    phase_plots_yrange = [-.045,.045] 
    
    ; Plot the overall gain phase, as fit by the plane fitting above
    cgps_open, output_path+'/'+obsid+'_cal_phase_errors_'+pol_name+'.png'
    cgWindow_SetDefs, IM_WIDTH=600
    cgplot, freq, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Overall Gain Phase', title='Overall Gain Phase $\Delta$', xrange=[min(freq), max(freq)], yrange=phase_plots_yrange, charsize=1.
    cgps_close, /png, /delete_ps
    
    ; Plot the average gain phase with the per-antenna gain phases
    ave_phase = mean(gains_phase, dimension=2., /nan)
    cgps_open, output_path+'/'+obsid+'_cal_ant_phase_errors_'+pol_name+'.png'
    cgplot, freq, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Phase', title='Antenna Gain Phases', xrange=[min(freq), max(freq)], yrange=phase_plots_yrange, charsize=1., /nodata
    for antenna = 0, (size(gains_phase))[2]-1 do begin
      cgplot, freq, gains_phase[*,antenna], PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    cgplot, freq, ave_phase, PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['Average Gain Phase', 'Per-Antenna Gain Phases'], psym=-0, color=['blue','gray'], length=0.03, /center_sym, location=legend_loc, charsize=1.
    cgps_close, /png, /delete_ps
    
    ; Plot Fourier Transform of average gain phase with per-antenna gain phases
    cgps_open, output_path+'/'+obsid+'_cal_ant_phase_errors_'+pol_name+'_ft.png'
    cgplot, [freq[0], freq[-1]], [1,1], PSym=-0, symsize=.5, color='blue', aspect=.8, xtitle='Calibration Delay (s)', $
      ytitle='Gain Phase Power', title='Antenna Gain Phases - Power Spectrum',$
      xrange=ft_x_range, yrange=[0,2e9], charsize=1., /nodata
    for antenna = 0, (size(gains_amp))[2]-1 do begin
      fourier, y_time=gains_phase_ft, time=x_axis_ft, y_freq=gains_phase[*, antenna], frequencies=freq*1.e6, /inverse
      cgplot, x_axis_ft, abs(gains_phase_ft)^2., PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    fourier, y_time=mean_phase_ft, time=x_axis_ft, y_freq=ave_phase, frequencies=freq*1.e6, /inverse
    cgplot, x_axis_ft, abs(mean_phase_ft)^2, PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['Average Gain Phase', 'Per-Antenna Gain Phases'], psym=-0, color=['blue','gray'], length=0.03, /center_sym, location=ft_legend_loc, charsize=1.
    cgps_close, /png, /delete_ps
        
    ; Plot the per-antenna gain phases with the average subtracted
    gains_phase_ave_remove = make_array((size(gains_phase))[1], (size(gains_phase))[2], /float, value=0.)
    for antenna = 0, (size(gains_phase))[2]-1 do begin
      gains_phase_ave_remove[*,antenna]=gains_phase[*,antenna]-ave_phase
    endfor
    gains_phase_ave_remove_x2 = total(gains_phase_ave_remove^2., 2)
    gains_phase_ave_remove_stddev = stddev(gains_phase_ave_remove, dimension=2)
    cgps_open, output_path+'/'+obsid+'_cal_ant_phase_errors_ave_remove_'+pol_name+'.png'
    cgplot, freq, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Phase Error', title='Per-Antenna Gain Phase Errors, Ave. Removed', xrange=[min(freq), max(freq)], yrange=phase_plots_yrange, charsize=1., /nodata
    for antenna = 0, (size(gains_phase))[2]-1 do begin
      cgplot, freq, gains_phase_ave_remove[*,antenna], PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    cgps_close, /png, /delete_ps
    
    ; Plot the per-antenna gain phases with the plane fit subtracted
    gains_phase_residual = make_array((size(gains_phase))[1], (size(gains_phase))[2], /float, value=0.)
    for antenna = 0, (size(gains_phase))[2]-1 do begin
      gains_phase_residual[*,antenna]=gains_phase[*,antenna]-fit_params[*,0]-fit_params[*,1]*ant_pos[antenna,0]-fit_params[*,2]*ant_pos[antenna,1]
    endfor
    gains_phase_residual_x2 = total(gains_phase_residual^2., 2)
    gains_phase_residual_stddev = stddev(gains_phase_residual, dimension=2)
    cgps_open, output_path+'/'+obsid+'_cal_ant_phase_errors_no_grad_'+pol_name+'.png'
    cgplot, freq, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Phase Error', title='Per-Antenna Gain Phase Errors, Surface Fit Subtracted', xrange=[min(freq), max(freq)], yrange=phase_plots_yrange, charsize=1., /nodata
    for antenna = 0, (size(gains_phase))[2]-1 do begin
      cgplot, freq, gains_phase_residual[*,antenna], PSym=-0, symsize=.5, color='gray', /overplot
    endfor
    cgps_close, /png, /delete_ps
  
    ; Plot the x and y components of the gain phase gradient
    cgps_open, output_path+'/'+obsid+'_cal_phase_grad_errors_'+pol_name+'.png'
    cgplot, freq, fit_params[*,1], PSym=-0, symsize=.5, color='red', aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Phase Gradient (1/m)', title='Gain Phase Gradient Parameters $\Delta$x and $\Delta$y', xrange=[min(freq), max(freq)], yrange=[-6e-6,6e-6], charsize=1.
    cgplot, freq, fit_params[*,2], PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['$\Delta$x (Gain Phase North-South Gradient)', '$\Delta$y (Gain Phase North-South Gradient)'], psym=-0, color=['red','blue'], length=0.03, /center_sym, location=[.4, legend_loc[1]], charsize=1.
    cgps_close, /png, /delete_ps
    
    ; Plot Fourier Transform of the x and y components of the gain phase gradient
    cgps_open, output_path+'/'+obsid+'_cal_ant_phase_grad_errors_'+pol_name+'_ft.png'
    fourier, y_time=gains_phase_grad_x_ft, time=x_axis_ft, y_freq=fit_params[*,1], frequencies=freq*1.e6, /inverse
    cgplot, x_axis_ft, abs(gains_phase_grad_x_ft)^2, PSym=-0, symsize=.5, color='red', aspect=.8, xtitle='Calibration Delay (s)', $
      ytitle='Gain Phase Gradient Power', title='Gain Phase Gradient Parameters $\Delta$x and $\Delta$y - Power Spectrum', xrange=ft_x_range, yrange=[0,15], charsize=1.
    fourier, y_time=gains_phase_grad_y_ft, time=x_axis_ft, y_freq=fit_params[*,2], frequencies=freq*1.e6, /inverse
    cgplot, x_axis_ft, abs(gains_phase_grad_y_ft)^2, PSym=-0, symsize=.5, color='blue', /overplot
    cglegend, title=['$\Delta$x (Gain Phase North-South Gradient)', '$\Delta$y (Gain Phase North-South Gradient)'], psym=-0, color=['red','blue'], length=0.03, /center_sym, location=ft_legend_loc, charsize=1.
    cgps_close, /png, /delete_ps
    
  endfor  

end