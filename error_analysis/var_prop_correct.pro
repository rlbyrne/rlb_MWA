
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For two values A and B sampled from Gaussian distributions with variances varA and varB, ;;
  ;; the variance of A*CONJ(B) is (varA + varB)^2 / 2 when varA EQ varB. When varA NE varB    ;;
  ;; the propagated variance differs from the above by a factor that depends on the ratio of  ;;
  ;; varA to varB. This function gives the correction factor for each input pair of variances ;;
  ;; where correction_factors = calculated variances / actual variances and the calculated    ;;
  ;; variances are calculated as above.                                                       ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION var_prop_correct, var1, var2

  var1 = FLOAT(var1)
  var2 = FLOAT(var2)
  
  var_proportion = var2 / (var1 + var2)
  RESTORE, '~/MWA/IDL_code/rlb_MWA/error_analysis/var_prop_corrections.sav'
  correction_factors = MAKE_ARRAY(N_ELEMENTS(var_proportion), /FLOAT)
  FOR i = 0, N_ELEMENTS(var_proportion) - 1 DO BEGIN
  if var_proportion[i] lt 0.02 or var_proportion[i] gt 0.98 then correction_factors[i] = 0 else $
    correction_factors[i] = INTERPOL(correction_factor, var_ratio, var_proportion[i], /NAN)
  ENDFOR
  
  RETURN, correction_factors
  
  
END