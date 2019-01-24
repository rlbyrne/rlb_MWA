pro plot_abscal_compare

  optimal_amp = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/optimal_abscal.sav', 'amp')
  skycal_amp = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/from_skycal_abscal.sav', 'amp')
  optimal_delta_x = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/optimal_abscal.sav', 'delta_x')
  skycal_delta_x = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/from_skycal_abscal.sav', 'delta_x')
  optimal_delta_y = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/optimal_abscal.sav', 'delta_y')
  skycal_delta_y = getvar_savefile('/Users/rubybyrne/abscal_method_comparison/from_skycal_abscal.sav', 'delta_y')
  
  ; Get frequency axis
  run_path = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018'
  cal = getvar_savefile(run_path+'/calibration/hex_array_sim_331_cal.sav', 'cal', /compatibility_mode)
  freq = cal.freq/1.e6
  
  plot_resolution = 600
  cgps_open, '/Users/rubybyrne/abscal_method_comparison/amplitude.png', XSIZE = 7, YSIZE = 7
  !P.Multi = [0, 1, 2]
  !X.OMargin = [2, 0]
  cgplot, freq, skycal_amp, linestyle=0, thick=1, color='blue', xrange=[min(freq), max(freq)], yrange=[.992,1.001], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Gain Amplitude'
  cgplot, freq, optimal_amp, linestyle=2, thick=4, color='blue', /overplot
  cglegend, title=['Abs. Cal. from Sky Cal.', 'Optimal Abs. Cal.'], linestyle=[0,2], thick=4, $
    color=['blue','blue'], length=0.03, /center_sym, location=[.18,.70], charsize=.8, /box, background='white'
  cgplot, freq, skycal_amp-optimal_amp, linestyle=0, thick=1, color='black', xrange=[min(freq), max(freq)], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Residual (from sky cal. - optimal)'
  cgps_close, /png, /delete_ps, density=plot_resolution

  cgps_open, '/Users/rubybyrne/abscal_method_comparison/phase_grad_x.png', XSIZE = 7, YSIZE = 7
  !P.Multi = [0, 1, 2]
  !X.OMargin = [2, 0]
  cgplot, freq, skycal_delta_x, linestyle=0, thick=1, color='blue', xrange=[min(freq), max(freq)], yrange=[-10e-5,6e-5], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Gain Phase Gradient (1/m)'
  cgplot, freq, optimal_delta_x, linestyle=2, thick=4, color='blue', /overplot
  cglegend, title=['E-W Phase Grad. from Sky Cal.', 'Optimal E-W Phase Grad.'], linestyle=[0,2], thick=4, $
    color=['blue','blue'], length=0.03, /center_sym, location=[.18,.70], charsize=.8, /box, background='white'
  cgplot, freq, skycal_delta_x-optimal_delta_x, linestyle=0, thick=1, color='black', xrange=[min(freq), max(freq)], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Residual (from sky cal. - optimal)', yrange=[-1e-6,1e-6]
  cgps_close, /png, /delete_ps, density=plot_resolution
  
  cgps_open, '/Users/rubybyrne/abscal_method_comparison/phase_grad_y.png', XSIZE = 7, YSIZE = 7
  !P.Multi = [0, 1, 2]
  !X.OMargin = [2, 0]
  cgplot, freq, skycal_delta_y, linestyle=0, thick=1, color='blue', xrange=[min(freq), max(freq)], yrange=[-10e-5,6e-5], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Gain Phase Gradient (1/m)'
  cgplot, freq, optimal_delta_y, linestyle=2, thick=4, color='blue', /overplot
  cglegend, title=['N-S Phase Grad. from Sky Cal.', 'Optimal N-S Phase Grad.'], linestyle=[0,2], thick=4, $
    color=['blue','blue'], length=0.03, /center_sym, location=[.18,.70], charsize=.8, /box, background='white'
  cgplot, freq, skycal_delta_y-optimal_delta_y, linestyle=0, thick=1, color='black', xrange=[min(freq), max(freq)], $
    charsize=1., xtitle='Frequency (MHz)', ytitle='Residual (from sky cal. - optimal)', yrange=[-1e-6,1e-6]
  cgps_close, /png, /delete_ps, density=plot_resolution
  
  stop

end