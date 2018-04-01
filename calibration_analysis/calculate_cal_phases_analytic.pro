pro calculate_cal_phases_analytic

  cal_file = '/home/rlbyrne/cal_simulation/1163765304_cal.sav'
  obs_file = '/home/rlbyrne/cal_simulation/1163765304_obs.sav'
  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=0.)
  gains[*,*,0] = (*cal.gain[0]) + (*cal.gain_residual[0])
  gains[*,*,1] = (*cal.gain[1]) + (*cal.gain_residual[1])
  ;gains_amp = abs(gains)
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
  
  ;stop
    
  cgps_open, '/home/rlbyrne/cal_simulation/cal_phase_errors_phaseII.png'
  ;cgWindow_SetDefs, IM_WIDTH=600
  cgplot, fit_params[*,0], PSym=16, symsize=.5, aspect=.8, xtitle='Frequency Channel', $
    ytitle='Average Gain Phase Error'
  ;cgplot, fit_params[*,1], PSym=16, symsize=.5, fit_params[*,1], color='red', aspect=.8, xtitle='Frequency Channel', $
  ;  ytitle='Average Gain Phase Gradient Error'
  ;cgplot, fit_params[*,2], PSym=16, symsize=.5, fit_params[*,1], color='blue', /overplot
  cgps_close, /png, /delete_ps
end
