PRO FastFourier

;domain parameters:
time_min = -5
time_max = 5
timestep = .001
n = FLOOR(time_max/timestep) - CEIL(time_min/timestep) + 1
time = (FINDGEN(n) + CEIL(time_min/timestep))*timestep

;establish input function:
func1 = EXP(-time^2)

;forward fourier transform
func1shifted = SHIFT(func1, CEIL(time_min/timestep))
func2unscaled = FFT(func1shifted) ;take Fourier Transform
func2 = SHIFT(func2unscaled*n*timestep, FLOOR(n/2)) ;scale amplitudes, shift domain
freq = (FINDGEN(n)-FLOOR(n/2))*2*!PI/(n*timestep) ;produce frequency domain

CGPLOT, time, func1, XRANGE = [time_min, time_max], YRANGE = [-.1, 2], COLOR = "blue"
;CGOPLOT, freq, func2, COLOR = "black"


;inverse fourier transform:
invfunc2 = func2
invfunc1shifted = FFT(invfunc2, /INVERSE, /CENTER)
invfunc1unscaled = SHIFT(invfunc1shifted, -CEIL(time_min/timestep))
invfunc1 = invfunc1unscaled/(n*timestep)
CGOPLOT, time, invfunc1, COLOR = "red", LINESTYLE = 2
CGOPLOT, time, func1-invfunc1, COLOR = "cyan", LINESTYLE = 2


END