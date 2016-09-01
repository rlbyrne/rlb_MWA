PRO invFFT1D, freq, func, t, func1

n = SIZE(freq, /N_ELEMENTS)
freqstep = (MAX(freq)-MIN(freq))/(n-1)

func2 = func
func2 = SHIFT(func2, ROUND(freq[0]/freqstep)) ;shift func2
func1 = FFT(func2, /INVERSE) ;take Fourier Transform
func1 = func1 * freqstep / (2*!PI) ;scale amplitudes
func1 = SHIFT(func1, FLOOR(n/2)) ;shift func1 such that it is centered around zero
t = (FINDGEN(n) - FLOOR(n/2)) * (2*!PI) / (n * freqstep) ;define "time" domain



END