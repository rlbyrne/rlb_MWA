PRO FFT2D, tX, tY, func, freqX, freqY, func2

nX = SIZE(tX, /N_ELEMENTS)
nY = SIZE(tY, /N_ELEMENTS)
tXstep = (MAX(tX)-MIN(tX))/(nX-1)
tYstep = (MAX(tY)-MIN(tY))/(nY-1)

func1 = func
func1 = SHIFT(func1, ROUND(tX[0]/tXstep), ROUND(tY[0]/tYstep)) ;shift func1
func2 = FFT(func1, /CENTER) ;take Fourier Transform, returns func2 centered around zero
func2 = func2*nX*nY*tXstep*tYstep ;scale amplitudes
freqX = (FINDGEN(nX) - FLOOR(nX/2)) * 2*!PI / (nX*tXstep) ;define "frequency" domain
freqY = (FINDGEN(nY) - FLOOR(nY/2)) * 2*!PI / (nY*tYstep)

END