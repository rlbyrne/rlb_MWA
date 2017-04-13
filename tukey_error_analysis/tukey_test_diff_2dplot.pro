pro tukey_test_diff_2dplot, tukey = tukey

  obsid1 = '1061316176'
  obsid2 = '1061316296'
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
  
  if keyword_set(tukey) then begin
    savefile_name_A = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/2d_diff_tukey_dirty_1061316176_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
    savefile_name_B = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/2d_diff_tukey_dirty_1061316296_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
  endif else begin
    savefile_name_A = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/2d_diff_tophat_dirty_1061316176_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
    savefile_name_B = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/2d_diff_tophat_dirty_1061316296_'+pol+'_term' + STRTRIM(STRING(choose_term), 2)
  endelse
  
  ACube = getvar_savefile(filename1, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  ACube = real_part(ACube)
  BCube = getvar_savefile(filename2, 'DATA_DIFF_' + STRTRIM(STRING(choose_term),2))
  BCube = real_part(BCube)
  kx_mpc = getvar_savefile(filename1, 'KX_MPC')
  ky_mpc = getvar_savefile(filename1, 'KY_MPC')
  kz_mpc = getvar_savefile(filename1, 'KZ_MPC')
  delay_params = getvar_savefile(filename1, 'DELAY_PARAMS')
  
  ACube_2d = kspace_rebinning_2d(ACube, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc_A, kpar_edges_mpc_A, kperp_bin = kperp_bin_A, kpar_bin = kpar_bin_A)
  BCube_2d = kspace_rebinning_2d(BCube, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc_B, kpar_edges_mpc_B, kperp_bin = kperp_bin_B, kpar_bin = kpar_bin_B)
  
  kpower_2d_plots, power = ACube_2d, delay_params = delay_params, $
    kperp_edges = kperp_edges_mpc_A, kpar_edges = kpar_edges_mpc_A, kperp_bin = kperp_bin_A, kpar_bin = kpar_bin_A, $
    plotfile = savefile_name_A, /PNG, FULL_TITLE = 'Even/Odd Diff', $
    NOTE = 'obsid 1061316176, '+ pol +', term' + STRTRIM(STRING(choose_term), 2), DATA_RANGE = [-1e7,1e7], color_profile = 'sym_log'

  kpower_2d_plots, power = BCube_2d, delay_params = delay_params, $
    kperp_edges = kperp_edges_mpc_B, kpar_edges = kpar_edges_mpc_B, kperp_bin = kperp_bin_B, kpar_bin = kpar_bin_B, $
    plotfile = savefile_name_B, /PNG, FULL_TITLE = 'Even/Odd Diff', $
    NOTE = 'obsid 1061316296, '+ pol +', term' + STRTRIM(STRING(choose_term), 2), DATA_RANGE = [-1e7,1e7], color_profile = 'sym_log'
  
end