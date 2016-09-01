PRO AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured

IF SIZE(X, /N_ELEMENTS) EQ 0 THEN BEGIN
   X = FINDGEN(100) - 50
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)

   ;place antennas
   array_X[25] = 1
   ;array_X[27] = 1
   ;array_X[53] = 1
   ;CGPLOT, X, array_X, TITLE = "Antenna Placement"

ENDIF

;PRINT, array_X

invFFT1D, X, array_X, thetaX, psf
WINDOW, 3
;CGPLOT, thetaX, ABS(psf), TITLE = "Array PSF"
PRINT, MAX(psf)

psf_sq = ABS(psf)^2

;CGOPLOT, thetaX, psf_sq


FFT1D, thetaX, psf_sq, U, array_U
CGPLOT, U, array_U, TITLE = "Array Baselines", XRANGE = [-2, 10], PSYM = -4


brightness_sky = MAKE_ARRAY(SIZE(thetaX, /N_ELEMENTS), /FLOAT, VALUE = 0)

;make stars:
brightness_sky[SIZE(thetaX, /N_ELEMENTS)/8] = 1
brightness_sky[SIZE(thetaX, /N_ELEMENTS)/6] = .8
brightness_sky[SIZE(thetaX, /N_ELEMENTS)*3/4] = 2


;CGPLOT, thetaX, brightness_sky, TITLE = "Brightness, Sky"

FFT1D, thetaX, brightness_sky, U, brightness_U
;WINDOW, 3
;CGPLOT, U, brightness_U, TITLE = "Intensity, Ground"
brightness_U_measured = brightness_U * array_U
invFFT1D, U, brightness_U_measured, thetaX, brightness_sky_measured
;CGPLOT, thetaX, brightness_sky_measured, TITLE = "Measured Sky Intensity"

END