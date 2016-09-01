PRO SkyBrightness

;;;;ARRAY (GROUND);;;;

X = FINDGEN(80) - 40
Y = FINDGEN(80) - 40
array_XY = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), SIZE(Y, /N_ELEMENTS), /FLOAT, VALUE = 0)

;place antennas
;array_XY[10, 10] = 1
;array_XY[10, 30] = 1
;array_XY[60, 30] = 1
;array_XY[60, 70] = 1

;FOR i = 0, 79 DO BEGIN
;   FOR j = 0, 79 DO BEGIN
;      IF X[i]^2+Y[j]^2 LE 9 THEN BEGIN
;         array_XY[i,j] = 1
;      ENDIF
;      IF X[i-10]^2+Y[j-10]^2 LE 9 THEN BEGIN
;         array_XY[i,j] = 1
;      ENDIF
;   ENDFOR
;ENDFOR


;CGSURFACE, array_XY, X_ground, Y_ground


;;;;ARRAY POINT SPREAD FUNCTION;;;;

invFFT2D, X, Y, array_XY, thetaX, thetaY, psf
;CGCONTOUR, psf, thetaX, thetaY, /FILL
;CGSURFACE, psf, thetaX, thetaY


;;;;ARRAY POINT SPREAD FUNCTION, SQUARED;;;;

psf_sq = ABS(psf)^2
;CGCONTOUR, psf_sq, thetaX, thetaY, /FILL


;;;;BASELINES (ARRAY IN UV SPACE);;;;

FFT2D, thetaX, thetaY, psf_sq, U, V, array_UV
;CGSURFACE, array_UV, U, V


;;;;SKY INTENSITY;;;;

brightness_sky = MAKE_ARRAY(SIZE(thetaX, /N_ELEMENTS), SIZE(thetaY, /N_ELEMENTS), /FLOAT, VALUE = 0)

;make stars:
brightness_sky[SIZE(thetaX, /N_ELEMENTS)/2, SIZE(thetaY, /N_ELEMENTS)*3/4] = 1
brightness_sky[SIZE(thetaX, /N_ELEMENTS)/3, SIZE(thetaY, /N_ELEMENTS)*5/8] = 1
brightness_sky[SIZE(thetaX, /N_ELEMENTS)*3/4, SIZE(thetaY, /N_ELEMENTS)/3] = 2

;CGCONTOUR, brightness_sky, thetaX, thetaY, /FILL
CGSURFACE, brightness_sky, thetaX, thetaY


;;;;GROUND INTENSITY;;;;

FFT2D, thetaX, thetaY, brightness_sky, U, V, brightness_UV
;CGCONTOUR, brightness_UV, U, V, /FILL


;;;;MEASURED VALUES;;;;

brightness_UV_measured = brightness_UV * array_UV
;CGCONTOUR, brightness_UV_measured, U, V, /FILL


;;;;MEASURED SKY INTENSITY;;;;

invFFT2D, U, V, brightness_UV_measured, thetaX, thetaY, brightness_sky_measured
;CGCONTOUR, brightness_sky_measured, thetaX, thetaY, /FILL
CGSURFACE, brightness_sky_measured, thetaX, thetaY


END