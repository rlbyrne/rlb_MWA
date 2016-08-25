PRO vissig_freq, restore_all=restore_all, $
    ;if restore_all is not set:
    obs = obs, vis_arr = vis_arr, cal = cal, flag_arr = flag_arr
    
  ; Script to plot visibility sigmas given various time interleaving or 2 second time visibilities
    
  if keyword_set(restore_all) then begin
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav' ;;restore obs structure
    vis_XX_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_XX.sav', 'vis_ptr') ;;restore array of calibrated visibilities
    vis_YY_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_YY.sav', 'vis_ptr')
    vis_arr = [vis_XX_restored, vis_YY_restored]
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_flags.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/1061316296_cal.sav'
  endif
  
  uv=sqrt(cal.UU^2.+cal.VV^2.) ;;baseline lengths
  pol_i=0
  flag_diff = ((*flag_arr[pol_i])>0)
  vis_arr_use = IMAGINARY(*vis_arr[pol_i])
  wh0 = WHERE(flag_diff EQ 0, count0)
  IF count0 GT 0 THEN vis_arr_use[wh0] = !VALUES.F_NAN
  
  sigma_arr=FLTARR(2,6)
  
  bi_use_even = (INDGEN(CEIL((SIZE(vis_arr_use))[1] / 2.)) * 2)[1:*]
  bi_use_odd = ((INDGEN(CEIL((SIZE(vis_arr_use))[1] / 2.)) * 2) + 1)[0:N_ELEMENTS(bi_use_even) - 1]
  data_diff = vis_arr_use[bi_use_even,*] - vis_arr_use[bi_use_odd,*] ; only use imaginary part
  
  
  ;FOR i = 1, (SIZE(*vis_arr[pol_i]))[1] - 1 DO BEGIN
  ;  data_diff =Imaginary( (*vis_arr[pol_i])[i,*])-Imaginary((*vis_arr[pol_i])[i-1,*])
  ;ENDFOR
  
  uv=sqrt(cal.UU^2.+cal.VV^2.)
  ;uv = uv[bi_use_even]
  
  binsize=10.
  result=histogram(uv*3.*10^8.,binsize=binsize,reverse_indices=ri, /NAN, locations=locations, omax=omax)
  vis_sigma=FLTARR(N_elements(result))
  num_vis_sigma = FLTARR(N_ELEMENTS(result))
  for i=0, N_elements(result)-1 do begin
    if result[i] GT 0 then begin
      data_diff_values = data_diff[*,ri[ri[i]:ri[i+1]-1]]
      vis_sigma[i]=stddev(data_diff_values, /NAN)/sqrt(2.)
      num_vis_sigma[i] = FLOAT(TOTAL(FINITE(data_diff_values))) / ((SIZE(data_diff))[1]-1)
    endif
  endfor
  y_arr=vis_sigma
  x_arr=locations+binsize/2
  
  cgPS_Open,'/nfs/eor-00/h1/rbyrne/vis_sig_analysis/vissig_freq_diff.png',/quiet,/nomatch
  cgplot, x_arr, y_arr, psym=10, Ystyle=8, xrange=[0,2000], title='Visibility Sigma vs Wavelength, Interleaved Frequencies ', ytitle='visibility sigma', charsize=1, position=[.15,.15,.85,.85]
  cgaxis, yaxis=1, yrange=[0,9000], /save, title='vis # in bin', charsize=1, color='blue'
  cgtext, .65,.75,' $\tex\sigma$ = '+strtrim(stddev(data_diff[*], /NAN)/sqrt(2.),2), /normal, charsize=1
  cgoplot, x_arr,num_vis_sigma, psym = 10, color='blue', xrange=[0,2000]
  cgPS_Close,/png,Density=300,Resize=100.,/allow_transparent,/nomessage
    
END