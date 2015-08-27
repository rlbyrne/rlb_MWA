FUNCTION hist_contour_levels, histogram_vals, input_contours, LABELS = LABELS

  contours = FLOAT(input_contours)
  
  IF SIZE(WHERE(contours GT 100 OR contours LT 0, /NULL), /N_ELEMENTS) NE 0 THEN BEGIN
    PRINT, 'ERROR: contour levels must be valid percentages between 0 and 100'
    RETURN, []
  ENDIF
  
  histvals = ROUND(histogram_vals)
  histvals = histvals[REVERSE(SORT(histvals))]
  contour_levels = MAKE_ARRAY(SIZE(contours, /N_ELEMENTS), VALUE = 0.)
  
  FOR i = 0, SIZE(contour_levels, /N_ELEMENTS) - 1 DO BEGIN
  
    num_above = ROUND(TOTAL(histvals) * contours[i] / 100)
    minindex = 0
    maxindex = SIZE(histvals, /N_ELEMENTS) - 1
    test = 0
    
    IF histvals[0] GT num_above THEN BEGIN
      test = 1
      contour_levels[i] = histvals[0]
    ENDIF
    IF TOTAL(histvals) LE num_above THEN BEGIN
      test = 1
      contour_levels[i] = histvals[SIZE(histvals, /N_ELEMENTS) - 1]
    ENDIF
    
    WHILE test EQ 0 DO BEGIN
      stepsize = CEIL((maxindex - minindex) / 10.)
      maxindex = minindex
      WHILE TOTAL( histvals[0 : MIN([maxindex,SIZE(histvals, /N_ELEMENTS) - 1])] ) LE num_above DO BEGIN
        minindex = maxindex
        maxindex = maxindex + stepsize
      ENDWHILE
      IF TOTAL( histvals[0 : MIN([maxindex - 1,SIZE(histvals, /N_ELEMENTS) - 1])] ) LE num_above THEN BEGIN
        test = 1
        contour_levels[i] = MEAN( histvals[MIN([maxindex - 1, SIZE(histvals, /N_ELEMENTS) - 1]) : MIN([maxindex, SIZE(histvals, /N_ELEMENTS) - 1])])
      ENDIF
    ENDWHILE
    
  ENDFOR
  
  contour_levels_sorted = contour_levels[UNIQ(contour_levels, SORT(contour_levels))]
  contours_sorted = ' ' + STRTRIM(STRING(FIX(contours[UNIQ(contour_levels, SORT(contour_levels))])), 2) + ' % '
  
  contour_levels = contour_levels_sorted
  contours = contours_sorted
 
 IF KEYWORD_SET(LABELS) THEN RETURN, contours ELSE RETURN, contour_levels 

  
END