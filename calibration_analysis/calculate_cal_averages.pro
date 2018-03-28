pro calculate_cal_averages

  cal_file = '/home/rlbyrne/cal_simulation/1061316296_cal.sav'
  cal = getvar_savefile(cal_file, 'cal', /compatibility_mode)
  gains = make_array((size(*cal.gain[0]))[1], (size(*cal.gain[0]))[2], 2, /complex, value=0.)
  gains[*,*,0] = (*cal.gain[0]) + (*cal.gain_residual[0])
  gains[*,*,1] = (*cal.gain[1]) + (*cal.gain_residual[1])
  gains_amp = abs(gains)
  gains_phase = atan(U_Q_mix_model[0],/phase)
  mean_amp = mean(gains_amp, dimension=3, /nan)  ;average over polarizations
  mean_amp = mean(mean_amp, dimension=2, /nan)  ;average over tiles
  mean_amp_dev = mean_amp-1.
  print, mean_amp_dev
  cgWindow_SetDefs, IM_WIDTH=600
  cgplot, mean_amp_dev, PSym=16, symsize=.5, aspect=.8, xtitle='Frequency Channel', $
    ytitle='Frac. Error in Average Gain Amp.', output='/home/rlbyrne/cal_simulation/cal_amp_errors_phaseI.png'
    
end