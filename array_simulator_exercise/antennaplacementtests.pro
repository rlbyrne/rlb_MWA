PRO AntennaPlacementTests

X = FINDGEN(100)-50

choose = 4

;one small antenna:
IF choose EQ -1 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   array_X[60] = 1
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 1
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, 2*MAX(brightness_sky_measured)], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF

;total coverage:
IF choose EQ -2 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 1)
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 1
   CGPLOT, thetaX, brightness_sky_measured, TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF


;;;;arrange 10 (delta function) antennas;;;;
antennas = 10

;antennas clumped together
IF choose EQ 1 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   array_X[45:54] = 4
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, X, array_X, YRANGE = [0,4.5], TITLE = 'ANTENNA PLACEMENT'
   WINDOW, 1
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 2
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, MAX([MAX(brightness_sky_measured), MAX(brightness_sky)])*1.1], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF

;antennas less clumped together
IF choose EQ 2 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   array_X[40] = 4
   array_X[42] = 4
   array_X[44] = 4
   array_X[46] = 4
   array_X[48] = 4
   array_X[50] = 4
   array_X[52] = 4
   array_X[54] = 4
   array_X[56] = 4
   array_X[58] = 4
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, X, array_X, YRANGE = [0,4.5], TITLE = 'ANTENNA PLACEMENT'
   WINDOW, 1
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 2
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, MAX([MAX(brightness_sky_measured), MAX(brightness_sky)])*1.1], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF

;antennas spread evenly
IF choose EQ 3 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   array_X[5] = 4
   array_X[5*3] = 4
   array_X[5*5] = 4
   array_X[5*7] = 4
   array_X[5*9] = 4
   array_X[5*11] = 4
   array_X[5*13] = 4
   array_X[5*15] = 4
   array_X[5*17] = 4
   array_X[5*19] = 4
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, X, array_X, YRANGE = [0,4.5], TITLE = 'ANTENNA PLACEMENT'
   WINDOW, 1
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 2
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, MAX([MAX(brightness_sky_measured), MAX(brightness_sky)])*1.1], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF

;antennas spread randomly
IF choose EQ 4 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   FOR i = 1, antennas DO BEGIN
      array_X[FIX(29*RANDOMU(seed))] = 4
   ENDFOR
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, X, array_X, YRANGE = [0,4.5], TITLE = 'ANTENNA PLACEMENT'
   WINDOW, 1
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 2
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, MAX([MAX(brightness_sky_measured), MAX(brightness_sky)])*1.1], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF


IF choose EQ 5 THEN BEGIN
   array_X = MAKE_ARRAY(SIZE(X, /N_ELEMENTS), /FLOAT, VALUE = 0)
   array_X[42] = 4
   array_X[43] = 4
   array_X[45] = 4
   array_X[46] = 4
   array_X[49] = 4
   array_X[50] = 4
   array_X[53] = 4
   array_X[54] = 4
   array_X[56] = 4
   array_X[57] = 4
   AntennaPlacement, X, array_X, U, array_U, thetaX, brightness_sky, brightness_sky_measured
   WINDOW, 0
   CGPLOT, X, array_X, YRANGE = [0,4.5], TITLE = 'ANTENNA PLACEMENT'
   WINDOW, 1
   CGPLOT, U, array_U, TITLE = "UV COVERAGE"
   WINDOW, 2
   CGPLOT, thetaX, brightness_sky_measured, YRANGE = [0, MAX([MAX(brightness_sky_measured), MAX(brightness_sky)])*1.1], TITLE = "MEASURED SKY BRIGHTNESS"
   CGOPLOT, thetaX, brightness_sky, COLOR = "blue"
ENDIF

END