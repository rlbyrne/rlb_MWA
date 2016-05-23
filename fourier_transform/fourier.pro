pro fourier, y_time=y_time, time=time, t_range=t_range, $
    y_freq=y_freq, frequencies=frequencies, freq_range=freq_range, inverse=inverse
  
  if ~keyword_set(inverse) then begin
    ; forward FT:
    ; inputs: y_time, time OR t_range
    ; outputs: y_freq, frequencies, freq_range
  
    n = n_elements(y_time)
    if n_elements(time) gt 0 then t_range = minmax(time)
    t_step = float(t_range[1]-t_range[0])/(n-1)
    if n_elements(time) lt 2 then time = (findgen(n)*t_step)+t_range[0]
    y_time_shift = shift(y_time, floor(t_range[0]/t_step))
    
    y_freq_unnorm = fft(y_time_shift,/center)
    y_freq = n*t_step*y_freq_unnorm
    is_n_even = (n mod 2) eq 0
    freq_step = 2*!Pi/(n*t_step)
    if is_n_even then frequencies = (findgen(n)-n/2)*freq_step else frequiencies = (findgen(n)-(n-1)/2)*freq_step
    freq_range = minmax(frequencies)
        
  endif else begin
    ; inverse FT:
    ; inputs: y_freq, freq OR freq_range, t_range
    ; outputs: y_time, time
  
    n = n_elements(y_freq)
    if n_elements(frequencies) gt 0 then freq_range = minmax(frequencies)
    freq_step = float(freq_range[1]-freq_range[0])/(n-1)
    y_time_unnorm = fft(y_freq, /center, /inverse)
    y_time_shift = freq_step/(2*!Pi)*y_time_unnorm
    t_step = 2*!Pi/(n*freq_step)
    time = (findgen(n)*t_step) + t_range[0]
    y_time = shift(y_time_shift, -floor(t_range[0]/t_step))
  endelse
  
end