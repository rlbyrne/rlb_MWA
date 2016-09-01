PRO FFT1D, t, func, freq, func2

n = SIZE(t, /N_ELEMENTS)
tstep = (MAX(t)-MIN(t))/(n-1)

func1 = func
func1 = SHIFT(func1, ROUND(t[0]/tstep)) ;shift func1
func2 = FFT(func1, /CENTER) ;take Fourier Transform, returns func2 centered around zero
func2 = func2*n*tstep ;scale amplitudes
freq = (FINDGEN(n) - FLOOR(n/2)) * 2*!PI / (n*tstep) ;define "frequency" domain

END