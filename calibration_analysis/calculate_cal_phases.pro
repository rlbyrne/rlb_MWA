pro phase_calc, ant_pos, fit_params, phase_calculated, phase_derivs

  phase_calculated = fit_params[0]+fit_params[1]*ant_pos[*,0]+fit_params[2]*ant_pos[*,1]
  IF N_PARAMS() GE 4 THEN $
    phase_derivs = [[make_array(128, /float, value=1.)], [ant_pos[*,0]], [ant_pos[*,1]]]
    
end

pro fit_phases

  cal_file = '/home/rlbyrne/cal_simulation/1061316296_cal.sav'
  obs_file = '/home/rlbyrne/cal_simulation/1061316296_obs.sav'
  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=0.)
  gains[*,*,0] = (*cal.gain[0]) + (*cal.gain_residual[0])
  gains[*,*,1] = (*cal.gain[1]) + (*cal.gain_residual[1])
  ;gains_amp = abs(gains)
  gains = mean(gains, dimension=3, /nan)  ; average over polarizations
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
  
  weights = make_array(128, /float, value=1.)  ; use equal weighting for now
  fit_params = [0.,0.,0.]
  
  stop
  
  for freq_ind=0,383 do begin
    phases_fitted = curvefit(ant_pos, reform(gains_phase[freq_ind,*]), weights, fit_params, function_name='PHASE_CALC')
    print, fit_params
  endfor
  
  
  stop
  
;cgWindow_SetDefs, IM_WIDTH=600
;cgplot, mean_amp_dev, PSym=16, symsize=.5, aspect=.8, xtitle='Frequency Channel', $
;  ytitle='Frac. Error in Average Gain Amp.', output='/home/rlbyrne/cal_simulation/cal_amp_errors_phaseI.png'
  
end

pro calculate_cal_phases
  fit_phases
end