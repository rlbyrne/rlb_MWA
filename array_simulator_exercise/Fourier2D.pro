PRO Fourier2D


;domain parameters:
timeX_min = -1
timeX_max = 1
timeXstep = .1
nX = FLOOR(timeX_max/timeXstep) - CEIL(timeX_min/timeXstep) + 1
timeX = (FINDGEN(nX) + CEIL(timeX_min/timeXstep))*timeXstep

timeY_min = -1
timeY_max = 1
timeYstep = .1
nY = FLOOR(timeY_max/timeYstep) - CEIL(timeY_min/timeYstep) + 1
timeY = (FINDGEN(nY) + CEIL(timeY_min/timeYstep))*timeYstep

;establish input function:
func1 = MAKE_ARRAY(nX, nY)
FOR x = 0, nX-1 DO BEGIN
   FOR y = 0, nY-1 DO BEGIN
      func1[x,y] = EXP((-timeX[x]^2-timeY[y]^2))
   ENDFOR
ENDFOR

;WINDOW, 1
;CGCONTOUR, func1, timeX, timeY, /ISOTROPIC
;CGSURFACE, func1, timeX, timeY

func1shifted = func1
FOR x = 0, nX-1 DO BEGIN
   func1shifted[x,*] = SHIFT(func1shifted[x,*], CEIL(timeX_min/timeXstep))
ENDFOR
FOR y = 0, nY-1 DO BEGIN
   func1shifted[*,y] = SHIFT(func1shifted[*,y], CEIL(timeY_min/timeYstep))
ENDFOR

;WINDOW, 2   
;CGCONTOUR, func1shifted, timeX, timeY, /ISOTROPIC
;CGSURFACE, func1shifted, timeX, timeY  

;forward fourier transform
func2unscaled = FFT(func1shifted, /CENTER) ;take Fourier Transform
func2 = func2unscaled*nX*nY*timeXstep*timeYstep ;scale amplitudes
freqX = (FINDGEN(nX) - FLOOR(nX/2)) *2*!PI / (nX*timeXstep)
freqY = (FINDGEN(nY) - FLOOR(nY/2)) *2*!PI / (nY*timeYstep)


;WINDOW, 3   
;CGCONTOUR, func2, freqX, freqY, /ISOTROPIC
;CGSURFACE, func2, freqX, freqY



;inverse fourier transform:
invfunc2unscaled = func2/(nX*nY*timeXstep*timeYstep) ;remove amplitude scaling (invfunc2unscaled = func2unscaled)
invfunc1shifted = FFT(invfunc2unscaled, /INVERSE, /CENTER)

invfunc1 = invfunc1shifted
FOR x = 0, nX-1 DO BEGIN
   invfunc1[x,*] = SHIFT(invfunc1[x,*], -CEIL(timeX_min/timeXstep))
ENDFOR
FOR y = 0, nY-1 DO BEGIN
   invfunc1[*,y] = SHIFT(invfunc1[*,y], -CEIL(timeY_min/timeYstep))
ENDFOR

;CGSURFACE, invfunc1-func1


END