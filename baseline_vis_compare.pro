PRO baseline_vis_compare, restore_all=restore_all, $
    ;if restore_all is not set:
    obs = obs, vis_arr = vis_arr, cal = cal, flag_arr = flag_arr
    
  if keyword_set(restore_all) then begin
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/metadata/1061316296_obs.sav' ;;restore obs structure
    vis_XX_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_XX.sav', 'vis_ptr') ;;restore array of calibrated visibilities
    vis_YY_restored = GETVAR_SAVEFILE('/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_vis_YY.sav', 'vis_ptr')
    vis_arr = [vis_XX_restored, vis_YY_restored]
    RESTORE, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/vis_data/1061316296_flags.sav'
    restore, '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_devel_June2015/calibration/1061316296_cal.sav'
  endif
  
  bin_start=(*obs.baseline_info).bin_offset
  baselines_uu = (cal.UU)[0:bin_start[1]-1]
  baselines_vv = (cal.VV)[0:bin_start[1]-1]
  
  pol_i=0
  freq = 5
  timestep = 28
  flag_diff = ((*flag_arr[pol_i])>0)
  vis_arr_use = IMAGINARY(*vis_arr[pol_i])
  wh0 = WHERE(flag_diff EQ 0, count0)
  IF count0 GT 0 THEN vis_arr_use[wh0] = !VALUES.F_NAN
  
  histmin = min([[baselines_uu],[baselines_vv]])
  histmax = max([[baselines_uu],[baselines_vv]])
  binsize = (histmax-histmin)/10.
  result_uu = histogram(baselines_uu, binsize=binsize, reverse_indices=ri_uu, /NAN, min = histmin, max = histmax)

  histvals = make_array(n_elements(result_uu), n_elements(result_uu), /integer, value = 0)
  vis_mean = histvals
  vis_stddev = histvals
  for i = 0, n_elements(result_uu) - 1 do begin
    if result_uu[i] GT 1 then begin
      result_vv = histogram(baselines_vv[ri_uu[ri_uu[i]:ri_uu[i+1]-1]], binsize = binsize, reverse_indices = ri_vv, /NAN, min = histmin, max = histmax)
      places = where(result_vv gt 1, count)
      histvals[i,*] = result_vv
    endif
    
    for j = 0, n_elements(result_vv)-1 do begin
      if result_uu[i] gt 1 and result_vv[j] gt 1 then indices = ([ri_uu[ri_uu[i]:ri_uu[i+1]-1]])[ri_vv[ri_vv[j]:ri_vv[j+1]-1]] else indices = -1
      print, 'i, j, equal' + string(i) + ', ' + string(j)
      print, indices[0: min([10,n_elements(indices)-1])]
      if n_elements(indices) gt 1 then vis_mean[i,j] = mean(vis_arr_use[freq, indices + bin_start[timestep]], /NAN) else vis_mean[i,j] = 0
      if n_elements(indices) gt 1 then vis_stddev[i,j] = stddev(vis_arr_use[freq, indices + bin_start[timestep]], /NAN) else vis_mean[i,j] = 0
      stop
    endfor
    
  endfor
  
  stop
  
END