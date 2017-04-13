pro tukey_test_kzslice_plot, tukey=tukey

  choose_term = 2
  pol = 'yy'
  if keyword_set(tukey) then begin
    note_part = 'tukey filter, pol ' + pol + ', term ' + strtrim(string(choose_term),2)
  endif else begin
    note_part = 'no tukey filter, pol ' + pol + ', term ' + strtrim(string(choose_term),2)
  endelse
  
  filepath = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/ps/'
  if keyword_set(tukey) then begin
    filename1 = filepath + '1061316176_cubeXX__even_odd_joint_tk_dirty_' + pol + '_averemove_bh_kcube.idlsave'
    filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_tk_dirty_' + pol + '_averemove_bh_kcube.idlsave'
  endif else begin
    filename1 = filepath + 'Combined_obs_1061316176_cubeXX__even_odd_joint_dirty_' + pol + '_averemove_bh_kcube.idlsave'
    filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_dirty_' + pol + '_averemove_bh_kcube.idlsave'
  endelse
  
  check_file1 = file_test(filename1)
  check_file2 = file_test(filename2)
  if check_file1 eq 0 or check_file2 eq 0 then begin
    print, '***ERROR*** : files not found'
  endif
  
  ACube = getvar_savefile(filename1, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  ACube = real_part(ACube)
  BCube = getvar_savefile(filename2, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  BCube = real_part(BCube)
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  kz_mpc = getvar_savefile(filename1, 'KZ_MPC')
  delay_params = getvar_savefile(filename1, 'DELAY_PARAMS')
  
  use_slice = 48
  ACube_slice = ACube[*,*,use_slice]
  BCube_slice = BCube[*,*,use_slice]
  
  if keyword_set(tukey) then begin
    saveloc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/kzslice_tukey_dirty_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
  endif else begin
    saveloc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/kzslice_tophat_dirty_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
  endelse
  
  CGPS_OPEN, saveloc, /FONT, XSIZE = 30, YSIZE = 10
    
  QUICK_IMAGE, ACube_slice, kx_mpc, ky_mpc, $
    DATA_ASPECT = .5, TITLE = 'Even/Odd Diff Slice, 1061316176', $
    XTITLE = 'kx', YTITLE = 'ky', window = 1,$
    start_multi_params = {nrow:1, ncol:2},$
    note = 'kz slice ' + STRTRIM(STRING(kz_mpc[use_slice]),2) + ', ' + pol + ', term ' + STRTRIM(STRING(choose_term), 2),$
    MULTI_POS = multi_pos, /NOERASE, data_range = [-1e5,1e5]
    
  QUICK_IMAGE, BCube_slice, kx_mpc, ky_mpc, $
    DATA_ASPECT = .5, TITLE = 'Even/Odd Diff Slice, 1061316296', $
    XTITLE = 'kx', YTITLE = 'ky', window = 1,$
    MULTI_POS = multi_pos[*,1], /NOERASE, data_range = [-1e5,1e5]
    
  CGPS_CLOSE, /PNG, /DELETE_PS
  undefine, multi_pos
  
end