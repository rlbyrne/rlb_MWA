pro calculate_cal_params

  run_paths = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018'
  ;run_paths = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_GLEAM_catalog_testing_negative_sources_Nov2018'
  obsids = ['hex_array_sim_331']
  ;obsids = ['hex_array_sim_331','split_hex_array_sim_331',$
  ;  'random1_array_sim_331','random2_array_sim_331','random3_array_sim_331','random4_array_sim_331',$
  ;  'random5_array_sim_331','random6_array_sim_331','random7_array_sim_331','random8_array_sim_331',$
  ;  'random9_array_sim_331','random10_array_sim_331']
  output_path = '/Users/rubybyrne/array_simulation_331/cal_params_plots'
  ;output_path = '/Users/rubybyrne/test_barry_effect_nsources/three_quarters_GLEAM_model/negative_sources'
  create_cal_files = 0
  cal_file_save_path = '/Users/rubybyrne/array_simulation_331/cal_files'
  
  ;array_colors = ['navy', 'purple', 'dark green']
  array_colors = ['blue', 'red', 'orange']
  
  for obs_index = 0,n_elements(obsids)-1 do begin
    obsid = obsids[obs_index]
    if obs_index le n_elements(array_colors)-1 then array_color=array_colors[obs_index] else array_color=array_colors[-1]
    if n_elements(run_paths) gt 1 then run_path = run_paths[obs_index] else run_path = run_paths
    per_ant_array_color=array_color
    plot_resolution = 600
    cal_file = run_path+'/calibration/'+obsid+'_cal.sav'
    ;obs_file = run_path+'/metadata/'+obsid+'_obs.sav'  ; Obs structure has the wrong antenna locations
    cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
    abs_cal_solutions = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=complex(0.,0.))
    amp_cal_solutions = make_array((size(*cal.gain[0]))[1], 2, /complex, value=complex(0.,0.))
    for pol = 0,1 do begin

      if pol eq 0 then pol_name = 'XX' else pol_name='YY'
     
      legend_loc = [.3,.2]
      ft_legend_loc = [.2,.35]
      ft_x_range = [1e-1,2]
      kpar_conversion_factor = 5.51e5
      
      ; Get antenna positions
      csv_data = read_csv('/Users/rubybyrne/array_simulation_331/'+obsid+'_antenna_locs.csv', header=csv_header)
      ant_pos = make_array(n_elements(csv_data.field1), 2, /float, value=0.)
      for ant_ind=0, n_elements(csv_data.field1)-1 do begin
        ant_pos[fix(csv_data.field1[ant_ind]), 0] = double(csv_data.field2[ant_ind])
        ant_pos[fix(csv_data.field1[ant_ind]), 1] = double(csv_data.field3[ant_ind])
      endfor
   
      max_bl_len = cal.max_cal_baseline
              
      gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], /complex, value=0.)
      gains[*,*] = (*cal.gain[pol]) + (*cal.gain_residual[pol])
      gains_amp = abs(gains)
      mean_amp = mean(gains_amp, dimension=2, /nan)  ;average over antennas
      mean_amp_dev = mean_amp-1.
      freq = cal.freq/1.e6
      
      if obs_index eq 0 and pol eq 0 then begin
        mean_amp_all_obs = make_array(n_elements(mean_amp), 2, n_elements(obsids), /complex, value=complex(0.,0.))
        freq_reference = freq
        mean_amp_all_obs_ft = make_array(n_elements(mean_amp), 2, n_elements(obsids), /complex, value=complex(0.,0.))
      endif
      mean_amp_all_obs[*, pol, obs_index] = mean_amp

      ; Plot average amplitude with per-antenna solutions
      cgps_open, output_path+'/'+obsid+'_cal_ant_amp_errors_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
      cgplot, freq, mean_amp, linestyle=0, thick=4, color=array_color, aspect=.8, xtitle='Frequency (MHz)', $
        ytitle='Gain Amplitude', $
        xrange=[min(freq), max(freq)], yrange=[.998,1.002], charsize=1., /nodata
      for antenna = 0, (size(gains_amp))[2]-1 do begin
        cgplot, freq, gains_amp[*, antenna], linestyle=0, thick=1, color='medium gray', /overplot
      endfor
      if array_color eq 'orange' then cgplot, freq, mean_amp, linestyle=0, thick=5, color='black', /overplot
      cgplot, freq, mean_amp, linestyle=0, thick=4, color=per_ant_array_color, /overplot
      cgplot, [freq[0], freq[-1]], [1,1], linestyle=0, color='black', thick=1, /overplot
      ;cglegend, title=['Average Gain Amplitude', 'Per-Antenna Gain Amplitudes'], $
      ;  psym=-0, color=[array_color,'medium gray'], thick=4, length=0.03, /center_sym, location=legend_loc, charsize=1.
      cgps_close, /png, /delete_ps, density=plot_resolution
      
      average_gain_var = strtrim(variance(mean_amp), 1)
      average_per_antenna_gain_var = strtrim(mean(variance(gains_amp, dimension=1)), 1)
      
      ; Plot average amplitude ft with per_antenna_solutions
      cgps_open, output_path+'/'+obsid+'_cal_ant_amp_errors_'+pol_name+'_ft.png', XSIZE = 5, YSIZE = 4
      fourier, y_time=mean_amp_ft, time=x_axis_ft, y_freq=mean_amp, frequencies=freq*1.e6, /inverse
      cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=4, color=array_color, $
        aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
        ytitle='Gain Amplitude Power (MHz$\up2$)', $
        xrange=ft_x_range, yrange=[1e2,1e9], charsize=1., /xlog, /ylog, /nodata
      for antenna = 0, (size(gains_amp))[2]-1 do begin
        fourier, y_time=gains_amp_ft, time=x_axis_ft, y_freq=gains_amp[*, antenna], frequencies=freq*1.e6, /inverse
        cgplot, x_axis_ft*kpar_conversion_factor, abs(gains_amp_ft)^2., linestyle=0, thick=1, color='medium gray', /overplot, /xlog, /ylog
      endfor
      if array_color eq 'orange' then cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=5, color='black', /xlog, /ylog, /overplot
      cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=4, color=per_ant_array_color, /xlog, /ylog, /overplot
      ;cglegend, title=['Average Gain Amplitude', 'Per-Antenna Gain Amplitudes'], psym=-0, thick=4, color=[array_color,'medium gray'], length=0.03, /center_sym, location=ft_legend_loc, charsize=1.
      cgps_close, /png, /delete_ps, density=plot_resolution
      
      if obs_index eq 0 and pol eq 0 then begin
        freq_reference_ft = x_axis_ft*kpar_conversion_factor
        max_bl_len_reference = max_bl_len/3.e8*kpar_conversion_factor
      endif
      mean_amp_all_obs_ft[*, pol, obs_index] = abs(mean_amp_ft)^2.
      
      ; Plot average amplitude
      cgps_open, output_path+'/'+obsid+'_cal_amp_errors_'+pol_name+'.png', XSIZE = 7, YSIZE = 7
      cgplot, freq, mean_amp, linestyle=0, thick=4, color=array_color, aspect=.8, xtitle='Frequency (MHz)', $
        ytitle='Gain Amplitude', $
        xrange=[min(freq), max(freq)], yrange=[.999,1.001], charsize=1.
      cgplot, [freq[0], freq[-1]], [1,1], linestyle=0, color='black', thick=1, /overplot
      cgps_close, /png, /delete_ps, density=plot_resolution
      
      ; Plot average amplitude fourier transform
      cgps_open, output_path+'/'+obsid+'_cal_amp_errors_'+pol_name+'_ft.png', XSIZE = 7, YSIZE = 7
      fourier, y_time=mean_amp_ft, time=x_axis_ft, y_freq=mean_amp, frequencies=freq*1.e6, /inverse
      cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=4, color=array_color, $
        aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
        ytitle='Gain Amplitude Power (MHz$\up2$)', $
        xrange=ft_x_range, yrange=[1e2,1e8], charsize=1., /xlog, /ylog
      cgplot, [max_bl_len/3.e8*kpar_conversion_factor, max_bl_len/3.e8*kpar_conversion_factor], [1,1e8], linestyle=0, color=array_color, thick=1, /overplot, /xlog, /ylog
      cgps_close, /png, /delete_ps, density=plot_resolution
      
      gains_phase = atan(gains,/phase)
      weights = make_array((size(gains_phase))[2], 2, /float, value=1.)  ; use equal weighting for now
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
    
      fit_params = make_array((size(gains_phase))[1], 3, /float, value=0.)
      for freq_ind=0, (size(gains_phase))[1]-1 do begin
        fit_vec = make_array(3, /float, value=0.)
        fit_vec[0] = total(weights[*,pol]*reform(gains_phase[freq_ind,*]), /nan)
        fit_vec[1] = total(weights[*,pol]*reform(gains_phase[freq_ind,*])*ant_pos[*,0], /nan)
        fit_vec[2] = total(weights[*,pol]*reform(gains_phase[freq_ind,*])*ant_pos[*,1], /nan)
    
        fit_params[freq_ind,*] = fit_mat_inv#fit_vec
    
      endfor
    
      phase_plots_yrange = [-.045,.045] 
    
      ; Plot the x and y components of the gain phase gradient
      cgps_open, output_path+'/'+obsid+'_cal_phase_grad_errors_'+pol_name+'.png', XSIZE = 7, YSIZE = 7
      cgplot, freq, fit_params[*,1], linestyle=2, thick=4, color=array_color, aspect=.8, xtitle='Frequency (MHz)', $
        ytitle='Gain Phase Gradient (1/m)', xrange=[min(freq), max(freq)], yrange=[-8e-5,6e-5], charsize=1.
      cgplot, freq, fit_params[*,2], linestyle=1, thick=4, color=array_color, /overplot
      cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='black', thick=1, /overplot
      cglegend, title=['Gain Phase E-W Gradient', 'Gain Phase N-S Gradient'], linestyle=[2,1], thick=4, $
        color=[array_color,array_color], length=0.03, /center_sym, location=[.18,.32], charsize=.8, /box, background='white'
      cgps_close, /png, /delete_ps, density=plot_resolution
            
      ; Plot Fourier Transform of the x and y components of the gain phase gradient
      cgps_open, output_path+'/'+obsid+'_cal_phase_grad_errors_'+pol_name+'_ft.png', XSIZE = 7, YSIZE = 7
      fourier, y_time=gains_phase_grad_x_ft, time=x_axis_ft, y_freq=fit_params[*,1], frequencies=freq*1.e6, /inverse
      cgplot, [max_bl_len/3.e8*kpar_conversion_factor, max_bl_len/3.e8*kpar_conversion_factor], [1e-5,1e5], PSym=-0, symsize=.5, color='black', aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
        ytitle='Gain Phase Gradient Power (MHz$\up2$/m$\up2$)', linestyle=2,$
        xrange=ft_x_range, yrange=[1e-3,1e5], charsize=1., /xlog, /ylog, /nodata
      cgplot, x_axis_ft*kpar_conversion_factor, abs(gains_phase_grad_x_ft)^2, linestyle=2, thick=4, color=array_color, /overplot, /xlog, /ylog
      fourier, y_time=gains_phase_grad_y_ft, time=x_axis_ft, y_freq=fit_params[*,2], frequencies=freq*1.e6, /inverse
      cgplot, x_axis_ft*kpar_conversion_factor, abs(gains_phase_grad_y_ft)^2, linestyle=1, thick=4, color=array_color, /overplot, /xlog, /ylog
      cgplot, [max_bl_len/3.e8*kpar_conversion_factor, max_bl_len/3.e8*kpar_conversion_factor], [1e-5,1e5], linestyle=0, thick=1, color=array_color, /overplot
      cglegend, title=['Gain Phase E-W Gradient', 'Gain Phase N-S Gradient'], linestyle=[2,1], thick=4, $
        color=[array_color,array_color], length=0.03, /center_sym, location=[.18,.32], charsize=.8, /box, background='white'
      cgps_close, /png, /delete_ps, density=plot_resolution
      
      for antenna = 0, (size(gains_phase))[2]-1 do begin
        phases = fit_params[*,1]*ant_pos[antenna,0]+fit_params[*,2]*ant_pos[antenna,1]  ; include phase gradient terms only
        abs_cal_solutions[*,antenna,pol] = complex(mean_amp*cos(phases), mean_amp*sin(phases))
      endfor
      amp_cal_solutions[*, pol] = complex(mean_amp, 0.)
      
    endfor
    
    if keyword_set(create_cal_files) then begin
      
      ; create cal file with errors from the average amplitude and phase gradients
      *cal.gain[0] = abs_cal_solutions[*,*,0]
      *cal.gain[1] = abs_cal_solutions[*,*,1]
      *cal.gain_residual[0] = make_array((size(*cal.gain_residual[0]))[1], (size(*cal.gain_residual[0]))[2], /complex, value=complex(0.,0.))
      *cal.gain_residual[1] = make_array((size(*cal.gain_residual[1]))[1], (size(*cal.gain_residual[1]))[2], /complex, value=complex(0.,0.))
      print, 'Saving file to '+cal_file_save_path+'/'+obsid+'_cal_abs_errors_only.sav'
      save, cal, filename=cal_file_save_path+'/'+obsid+'_cal_abs_errors_only.sav'
      
      ; create cal file with errors from the average amplitude only
      *cal.gain[0] = complex(rebin(real_part(reform(amp_cal_solutions[*,0])), (size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2]), $
        rebin(imaginary(reform(amp_cal_solutions[*,0])), (size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2]))
      *cal.gain[1] = complex(rebin(real_part(reform(amp_cal_solutions[*,1])), (size(*cal.gain[1]))[1], (size(*cal.gain[1]))[2]), $
        rebin(imaginary(reform(amp_cal_solutions[*,1])), (size(*cal.gain[1]))[1], (size(*cal.gain[1]))[2]))
      print, 'Saving file to '+cal_file_save_path+'/'+obsid+'_cal_amp_errors_only.sav'
      save, cal, filename=cal_file_save_path+'/'+obsid+'_cal_amp_errors_only.sav'
    endif
            
  endfor  
  
  for pol = 0,1 do begin

    if pol eq 0 then pol_name = 'XX' else pol_name='YY'
     
    cgps_open, output_path+'/all_obs_cal_amp_errors_'+pol_name+'.png', XSIZE = 7, YSIZE = 7
    cgplot, freq_reference, reform(mean_amp_all_obs[*, pol, 0]), linestyle=0, thick=4, color=array_color, aspect=.8, xtitle='Frequency (MHz)', $
      ytitle='Gain Amplitude', $
      xrange=[min(freq_reference), max(freq_reference)], yrange=[.992,1.001], charsize=1., /nodata
    ;cgplot, freq_reference,  reform(mean_amp_all_obs[*, pol, 0]), linestyle=0, thick=4, color='blue', /overplot
    ;cgplot, freq_reference,  reform(mean_amp_all_obs[*, pol, 0]), linestyle=0, thick=4, color='blue', /overplot
    for obs_index = 0,n_elements(obsids)-1 do begin
      if obs_index le n_elements(array_colors)-1 then array_color=array_colors[obs_index] else array_color=array_colors[-1]
      if obs_index le 1 then thick=4 else if obs_index le 2 then thick=5 else thick=.5
      cgplot, freq_reference,  reform(mean_amp_all_obs[*, pol, obs_index]), linestyle=0, thick=thick, color=array_color, /overplot
    endfor
    cgplot, [freq_reference[0], freq_reference[-1]], [1,1], linestyle=0, color='black', thick=1, /overplot
    cglegend, title=['Hexagonal Array', 'Offset Hexagonal Array', 'Random Arrays'], $
      psym=-0, color=array_colors, thick=4, length=0.03, /center_sym, location=[.18,.35], charsize=.8, /box, background='white'
    cgps_close, /png, /delete_ps, density=plot_resolution
    
    cgps_open, output_path+'/all_obs_cal_amp_errors_'+pol_name+'_ft.png', XSIZE = 7, YSIZE = 7
    cgplot, freq_reference_ft, reform(mean_amp_all_obs_ft[*, pol, 0]), linestyle=0, thick=4, color=array_color, $
      aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
      ytitle='Gain Amplitude Power (MHz$\up2$)', $
      xrange=ft_x_range, yrange=[1,1e8], charsize=1., /xlog, /ylog, /nodata
    for obs_index = 0,n_elements(obsids)-1 do begin
      if obs_index le n_elements(array_colors)-1 then array_color=array_colors[obs_index] else array_color=array_colors[-1]
      if obs_index le 1 then thick=4 else if obs_index le 2 then thick=5 else thick=.5
      cgplot, freq_reference_ft,  reform(mean_amp_all_obs_ft[*, pol, obs_index]), linestyle=0, thick=thick, color=array_color, /overplot
    endfor
    cgplot, [max_bl_len_reference, max_bl_len_reference], [1,1e8], linestyle=0, color=array_colors[0], thick=1, /overplot, /xlog, /ylog
    cglegend, title=['Hexagonal Array', 'Offset Hexagonal Array', 'Random Arrays'], $
      psym=-0, color=array_colors, thick=4, length=0.03, /center_sym, location=[.18,.35], charsize=.8, /box, background='white'
    cgps_close, /png, /delete_ps, density=plot_resolution
    
  endfor

end