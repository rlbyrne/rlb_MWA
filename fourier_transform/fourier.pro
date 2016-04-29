pro fourier

trange = [0,100]
tstep = .3
t = (findgen(round((trange[1]-trange[0])/tstep)+1)*tstep)+trange[0]
n = n_elements(t)
;y_time = make_array(n,value=5)
y_time = exp(-(t-50)^2/100)

y_freq_unnorm = fft(y_time,/center)
y_freq = n*tstep*y_freq_unnorm

stop
is_n_even = (n mod 2) eq 0
freq_step = 2*!Pi/(n*tstep)
if is_n_even then freq = [findgen(n)-n/2]*freq_step $
  else freq = [findgen(n)-(n-1)/2]*freq_step
  
y_time_new = fft(y_freq, /center, /inverse)
tstep_new = 2*!Pi/(n*freq_step)
t_new = findgen(n)*tstep_new
  
stop
end