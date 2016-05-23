pro fourier_wrapper

  n = 100
  t_range = [-10,10]
  t_vals = findgen(n)*float(t_range[1]-t_range[0])/(n-1)+t_range[0]
  y_time = sin(t_vals/5)
  
  fourier, y_time=y_time, t_range=t_range, time=time, y_freq=y_freq, frequencies=frequencies
  fourier, y_time=y_time_new, t_range=t_range, time=time_new, y_freq=y_freq, frequencies=frequencies, /inverse
  cgplot, time, y_time, thick=2
  cgplot, time_new, real_part(y_time_new), /overplot, color='red', thick=2
  
  clip_vals = [6,2,3]
  colors = ['blue','green','lavender','turquoise']
  for clip = 0, n_elements(clip_vals)-1 do begin
    y_freq_clip = y_freq[n/2-n/(2*clip_vals[clip]):n/2+n/(2*clip_vals[clip])]
    frequencies_clip = frequencies[n/2-n/(2*clip_vals[clip]):n/2+n/(2*clip_vals[clip])]
    
    fourier, y_time=y_time_clip, t_range=t_range, time=time_clip, y_freq=y_freq_clip, frequencies=frequencies_clip, /inverse
    cgplot, time_clip, real_part(y_time_clip), /overplot, color=colors[clip], psym=7, symsize=1, thick=2
    
    undefine, y_time_clip
    undefine, time_clip
  endfor
  stop
  
end