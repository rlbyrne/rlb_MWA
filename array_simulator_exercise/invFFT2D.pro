PRO invFFT2D, freqX, freqY, func, tX, tY, func1

nX = SIZE(freqX, /N_ELEMENTS)
nY = SIZE(freqY, /N_ELEMENTS)
freqXstep = (MAX(freqX)-MIN(freqX))/(nX-1)
freqYstep = (MAX(freqY)-MIN(freqY))/(nY-1)

func2 = func
func2 = SHIFT(func2, ROUND(freqX[0]/freqXstep), ROUND(freqY[0]/freqYstep)) ;shift func2
func1 = FFT(func2, /INVERSE) ;take Fourier Transform
func1 = func1 * freqXstep * freqYstep / (2*!PI)^2 ;scale amplitudes
func1 = SHIFT(func1, FLOOR(nX/2), FLOOR(nY/2)) ;shift func1 such that it is centered around zero
tX = (FINDGEN(nX) - FLOOR(nX/2)) * (2*!PI) / (nX * freqXstep) ;define "time" domain
tY = (FINDGEN(nY) - FLOOR(nY/2)) * (2*!PI) / (nY * freqYstep)


END