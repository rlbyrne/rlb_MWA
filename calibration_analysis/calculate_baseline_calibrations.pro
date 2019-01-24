pro calculate_baseline_calibrations

  run_path = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_traditional_cal_large_window_Jul2018'
  obsid = 'hex_array_sim_331'
  output_path = '/Users/rubybyrne/array_simulation_331/cal_params_plots'
  pol = 0
  plot_resolution = 600
  kpar_conversion_factor = 16500.
  ft_x_range = [3e-3,1]

  cal_file = run_path+'/calibration/'+obsid+'_cal.sav'
  ;obs_file = run_path+'/metadata/'+obsid+'_obs.sav'
  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  abs_cal_solutions = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=complex(0.,0.))
  amp_cal_solutions = make_array((size(*cal.gain[0]))[1], 2, /complex, value=complex(0.,0.))

  ; Get antenna positions
  csv_data = read_csv('/Users/rubybyrne/array_simulation_331/'+obsid+'_antenna_locs.csv', header=csv_header)
  ant_pos = make_array(n_elements(csv_data.field1), 2, /float, value=0.)
  for ant_ind=0, n_elements(csv_data.field1)-1 do begin
    ant_pos[fix(csv_data.field1[ant_ind]), 0] = double(csv_data.field2[ant_ind])
    ant_pos[fix(csv_data.field1[ant_ind]), 1] = double(csv_data.field3[ant_ind])
  endfor

  max_bl_len = cal.max_cal_baseline
  if pol eq 0 then pol_name = 'XX' else pol_name='YY'

  gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], /complex, value=0.)
  abs_cal_solutions = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], /complex, value=0.)
  gains[*,*] = (*cal.gain[pol]) + (*cal.gain_residual[pol])
  gains_amp = abs(gains)
  mean_amp = mean(gains_amp, dimension=2, /nan)  ;average over antennas
  freq = cal.freq/1.e6
  
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
  
  for antenna = 0, (size(gains_phase))[2]-1 do begin
    phases = fit_params[*,1]*ant_pos[antenna,0]+fit_params[*,2]*ant_pos[antenna,1]  ; include phase gradient terms only
    abs_cal_solutions[*,antenna] = complex(mean_amp*cos(phases), mean_amp*sin(phases))
  endfor
  
  ;Find baselines:
  baseline_x = []
  baseline_y = []
  baseline_ant1 = []
  baseline_ant2 = []
  which_baseline = []
  for ant1 = 0, (size(gains_phase))[2]-2 do begin
    for ant2 = ant1+1, (size(gains_phase))[2]-1 do begin
      xdiff = ant_pos[ant2,0]-ant_pos[ant1,0] 
      ydiff = ant_pos[ant2,1]-ant_pos[ant1,1]
      if ydiff lt 0 then begin
        xdiff = -xdiff
        ydiff = -ydiff
      endif
      if xdiff^2.+ydiff^2. lt 16.^2. then begin
        baseline_x = [baseline_x, xdiff]
        baseline_y = [baseline_y, ydiff]
        baseline_ant1 = [baseline_ant1, ant1]
        baseline_ant2 = [baseline_ant2, ant1]
        if round(xdiff*10.)/10. eq -7.5 then which_baseline = [which_baseline, 1] else begin
          if round(xdiff*10.)/10. eq 7.5 then which_baseline = [which_baseline, 2] else begin
            if round(xdiff*10.)/10. eq 15 then which_baseline = [which_baseline, 3] else print, 'ERROR: INVALID BASELINE'
          endelse
        endelse
      endif
    endfor
  endfor
  
  baseline_gains = gains[*,baseline_ant1]*gains[*,baseline_ant2]
  baseline_gains_abs = abs_cal_solutions[*,baseline_ant1]*abs_cal_solutions[*,baseline_ant2]
    
  ;PLot the per-baseline gain amplitudes for traditional cal and absolute cal
  cgps_open, output_path+'/'+obsid+'_baseline_gain_amp_'+pol_name+'.png'
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Gain Amplitude', $
    xrange=[min(freq), max(freq)], yrange=[.97,1.016], charsize=1., /nodata
  for bl = 0, n_elements(baseline_x)-1 do begin
    cgplot, freq, abs(baseline_gains[*, bl]), linestyle=0, thick=1, color='blue', /overplot
  endfor
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='black', /overplot
  cgplot, [freq[0], freq[-1]], [1,1], linestyle=0, color='grey', thick=1, /overplot
  cglegend, title=['15m baselines, per-antenna cal', 'abs cal'], linestyle=[0,0], thick=1, color=['blue','black'], length=0.03, /center_sym, location=[.3,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=plot_resolution

  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_amp_'+pol_name+'.png'
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Gain Amplitude', $
    xrange=[min(freq), max(freq)], yrange=[.97,1.016], charsize=1., /nodata
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 1)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 1)]), dimension=2), linestyle=0, thick=1, color='blue', /overplot
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 2)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 2)]), dimension=2), linestyle=0, thick=1, color='red', /overplot
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 3)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 3)]), dimension=2), linestyle=0, thick=1, color='orange', /overplot
  cgplot, freq, mean_amp^2., linestyle=0, thick=1, color='black', /overplot
  cgplot, [freq[0], freq[-1]], [1,1], linestyle=0, color='grey', thick=1, /overplot
  cglegend, title=['(-7.5m,13m) baselines', '(7.5m,13m) baselines', '(15m,0m) baselines', 'all baselines'],$
    linestyle=[0,0,0,0], thick=1, $
    color=['blue','red','orange', 'black'], length=0.03, /center_sym, location=[.3,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=plot_resolution

  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_amp_errors_'+pol_name+'.png'
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Gain Amplitude Error^2', $
    xrange=[min(freq), max(freq)], yrange=[-3e-6,2.5e-4], charsize=1., /nodata
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 1)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 1)]), dimension=2)-1)^2., linestyle=0, thick=1, color='blue', /overplot
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 2)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 2)]), dimension=2)-1)^2., linestyle=0, thick=1, color='red', /overplot
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 3)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 3)]), dimension=2)-1)^2., linestyle=0, thick=1, color='orange', /overplot
  cgplot, freq, (mean_amp^2.-1)^2., linestyle=0, thick=1, color='black', /overplot
  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
  cglegend, title=['(-7.5m,13m) baselines', '(7.5m,13m) baselines', '(15m,0m) baselines', 'all baselines'],$
    linestyle=[0,0,0,0], thick=1, $
    color=['blue','red','orange', 'black'], length=0.03, /center_sym, location=[.3,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=plot_resolution
  
  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_amp_errors_diff_'+pol_name+'.png'
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Gain Amplitude Error^2', $
    xrange=[min(freq), max(freq)], yrange=[-1e-5,4e-5], charsize=1., /nodata
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 1)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 1)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2., linestyle=0, thick=1, color='blue', /overplot
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 2)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 2)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2., linestyle=0, thick=1, color='red', /overplot
  cgplot, freq, (mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 3)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 3)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2., linestyle=0, thick=1, color='orange', /overplot
  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
  cglegend, title=['(-7.5m,13m) baselines', '(7.5m,13m) baselines', '(15m,0m) baselines'],$
    linestyle=[0,0,0], thick=1, $
    color=['blue','red','orange'], length=0.03, /center_sym, location=[.3,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=plot_resolution
  
  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_amp_'+pol_name+'.png'
  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
    ytitle='Gain Amplitude', $
    xrange=[min(freq), max(freq)], yrange=[.97,1.016], charsize=1., /nodata
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 1)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 1)]), dimension=2), linestyle=0, thick=1, color='blue', /overplot
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 2)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 2)]), dimension=2), linestyle=0, thick=1, color='red', /overplot
  cgplot, freq, mean(abs((gains[*,baseline_ant1])[*,where(which_baseline eq 3)]), dimension=2)*mean(abs((gains[*,baseline_ant2])[*,where(which_baseline eq 3)]), dimension=2), linestyle=0, thick=1, color='orange', /overplot
  cgplot, freq, mean_amp^2., linestyle=0, thick=1, color='black', /overplot
  cgplot, [freq[0], freq[-1]], [1,1], linestyle=0, color='grey', thick=1, /overplot
  cglegend, title=['(-7.5m,13m) baselines', '(7.5m,13m) baselines', '(15m,0m) baselines', 'all baselines'],$
    linestyle=[0,0,0,0], thick=1, $
    color=['blue','red','orange', 'black'], length=0.03, /center_sym, location=[.3,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the per-baseline gain phases for traditional cal and absolute cal
  ;  cgps_open, output_path+'/'+obsid+'_baseline_gain_phase_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
  ;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
  ;    ytitle='Gain Amplitude', $
  ;    xrange=[min(freq), max(freq)], yrange=[-.05,.05], charsize=1., /nodata
  ;  for bl = 0, n_elements(baseline_x)-1 do begin
  ;    cgplot, freq, atan(baseline_gains[*, bl], /phase), linestyle=0, thick=1, color='blue', /overplot
  ;  endfor
  ;  for bl = 0, n_elements(baseline_x)-1 do begin
  ;    cgplot, freq, atan(baseline_gains_abs[*, bl], /phase), linestyle=0, thick=1, color='black', /overplot
  ;  endfor
  ;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
  ;  cgps_close, /png, /delete_ps, density=plot_resolution
  
  ;Plot the average per-baseline gain amplitudes for each baseline type for traditional cal - absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_amp_diff_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[-1e-5,5e-5], charsize=1., /nodata
;  cgplot, freq, (mean(abs(baseline_gains[*, where(which_baseline eq 1)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2, linestyle=0, thick=1, color='blue', /overplot
;  cgplot, freq, (mean(abs(baseline_gains[*, where(which_baseline eq 2)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2, linestyle=0, thick=1, color='red', /overplot
;  cgplot, freq, (mean(abs(baseline_gains[*, where(which_baseline eq 3)]), dimension=2)-1)^2.-(mean_amp^2.-1)^2, linestyle=0, thick=1, color='orange', /overplot
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the per-baseline gain amplitudes FT for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_amp_'+pol_name+'_ft.png', XSIZE = 5, YSIZE = 4
;  fourier, y_time=mean_amp_ft, time=x_axis_ft, y_freq=mean_amp^2., frequencies=freq*1.e6, /inverse
;  cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=4, color=array_color, $
;    aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
;    ytitle='Gain Amplitude Power (MHz$\up2$)', $
;    xrange=ft_x_range, yrange=[1e-1,1e10], charsize=1., /xlog, /ylog, /nodata
;  baseline_gain_amp_ft = make_array((size(baseline_gains))[1], (size(baseline_gains))[2], /float, value=0.)
;  for bl = 0, n_elements(baseline_x)-1 do begin
;    fourier, y_time=gains_amp_ft, time=x_axis_ft, y_freq=reform(abs(baseline_gains[*, bl])), frequencies=freq*1.e6, /inverse
;    baseline_gain_amp_ft[*,bl] = gains_amp_ft
;    cgplot, x_axis_ft*kpar_conversion_factor, abs(baseline_gain_amp_ft[*,bl])^2., linestyle=0, thick=1, color='black', /overplot, /xlog, /ylog
;  endfor
;  cgplot, x_axis_ft*kpar_conversion_factor, mean(abs(baseline_gain_amp_ft)^2., dimension=2, /nan), linestyle=0, thick=4, color='red', /overplot
;  cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2., linestyle=0, thick=4, color='blue', /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the square of the difference between the errors in per-baseline gain amplitudes for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_amp_avg_error_diff_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[-.0001,.0007], charsize=1., /nodata
;  for bl=0, n_elements(baseline_x)-1 do begin
;    cgplot, freq, reform(((abs(baseline_gains)-1.)^2.)[*, bl])-(mean_amp-1)^2., linestyle=0, thick=1, color='black', /overplot
;  endfor
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution
  
  ;PLot the per-baseline gain phases for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_phase_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[-.05,.05], charsize=1., /nodata
;  for bl = 0, n_elements(baseline_x)-1 do begin
;    cgplot, freq, atan(baseline_gains[*, bl], /phase), linestyle=0, thick=1, color='black', /overplot
;  endfor
;  for bl = 0, n_elements(baseline_x)-1 do begin
;    cgplot, freq, atan(baseline_gains_abs[*, bl], /phase), linestyle=0, thick=1, color='blue', /overplot
;  endfor
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the per-baseline gain phases FT for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_phase_'+pol_name+'_ft.png', XSIZE = 5, YSIZE = 4
;  fourier, y_time=mean_amp_ft, time=x_axis_ft, y_freq=mean_amp^2., frequencies=freq*1.e6, /inverse
;  cgplot, x_axis_ft*kpar_conversion_factor, abs(mean_amp_ft)^2, linestyle=0, thick=4, color=array_color, $
;    aspect=.8, xtitle='Calibration Delay (h/Mpc)', $
;    ytitle='Gain Amplitude Power (MHz$\up2$)', $
;    xrange=ft_x_range, yrange=[1e-1,1e10], charsize=1., /xlog, /ylog, /nodata
;  baseline_gain_phase_ft = make_array((size(baseline_gains))[1], (size(baseline_gains))[2], /float, value=0.)
;  for bl = 0, n_elements(baseline_x)-1 do begin
;    fourier, y_time=gains_phase_ft, time=x_axis_ft, y_freq=reform(atan(baseline_gains[*, bl], /phase)), frequencies=freq*1.e6, /inverse
;    baseline_gain_phase_ft[*,bl] = gains_phase_ft
;    cgplot, x_axis_ft*kpar_conversion_factor, abs(gains_phase_ft)^2., linestyle=0, thick=1, color='grey', /overplot, /xlog, /ylog
;  endfor
;  baseline_gain_phase_abs_ft = make_array((size(baseline_gains))[1], (size(baseline_gains))[2], /float, value=0.)
;  for bl = 0, n_elements(baseline_x)-1 do begin
;    fourier, y_time=gains_phase_ft, time=x_axis_ft, y_freq=reform(atan(baseline_gains_abs[*, bl], /phase)), frequencies=freq*1.e6, /inverse
;    baseline_gain_phase_abs_ft[*,bl] = gains_phase_ft
;    cgplot, x_axis_ft*kpar_conversion_factor, abs(gains_phase_ft)^2., linestyle=0, thick=1, color='black', /overplot, /xlog, /ylog
;  endfor
;  cgplot, x_axis_ft*kpar_conversion_factor, mean(abs(baseline_gain_phase_ft)^2., dimension=2, /nan), linestyle=0, thick=4, color='red', /overplot
;  cgplot, x_axis_ft*kpar_conversion_factor, mean(abs(baseline_gain_phase_abs_ft)^2., dimension=2, /nan), linestyle=0, thick=4, color='blue', /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the average of the square of the errors in per-baseline gain phases for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_phase_avg_error_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[0,.0004], charsize=1., /nodata
;  cgplot, freq, mean(atan(baseline_gains[*, *], /phase)^2., dimension=2), linestyle=0, thick=1, color='black', /overplot
;  cgplot, freq, mean(atan(baseline_gains_abs[*, *], /phase)^2., dimension=2), linestyle=0, thick=1, color='blue', /overplot
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

;PLot the average of the square of the errors in per-baseline gain phases for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_type_gain_phase_avg_error_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[0,.0004], charsize=1., /nodata
;  cgplot, freq, mean(atan(baseline_gains[*, where(which_baseline eq 1)], /phase)^2., dimension=2), linestyle=0, thick=1, color='blue', /overplot
;  cgplot, freq, mean(atan(baseline_gains[*, where(which_baseline eq 2)], /phase)^2., dimension=2), linestyle=0, thick=1, color='red', /overplot
;  cgplot, freq, mean(atan(baseline_gains[*, where(which_baseline eq 3)], /phase)^2., dimension=2), linestyle=0, thick=1, color='orange', /overplot
;  cgplot, freq, mean(atan(baseline_gains_abs[*, *], /phase)^2., dimension=2), linestyle=0, thick=1, color='black', /overplot
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  ;PLot the square of the difference between the errors in per-baseline gain phases for traditional cal and absolute cal
;  cgps_open, output_path+'/'+obsid+'_baseline_gain_phase_error_diff_'+pol_name+'.png', XSIZE = 5, YSIZE = 4
;  cgplot, freq, mean_amp^2., linestyle=0, thick=4, color='blue', aspect=.8, xtitle='Frequency (MHz)', $
;    ytitle='Gain Amplitude', $
;    xrange=[min(freq), max(freq)], yrange=[-.0005,.0013], charsize=1., /nodata
;  avg_deviation = mean((abs(baseline_gains)-1.)^2., dimension=2)
;  for bl=0, n_elements(baseline_x)-1 do begin
;    cgplot, freq, reform(atan(baseline_gains[*, bl], /phase)^2.)-reform(atan(baseline_gains_abs[*, bl], /phase)^2.), linestyle=0, thick=1, color='black', /overplot
;  endfor
;  cgplot, [freq[0], freq[-1]], [0,0], linestyle=0, color='grey', thick=1, /overplot
;  cgps_close, /png, /delete_ps, density=plot_resolution

  stop

end