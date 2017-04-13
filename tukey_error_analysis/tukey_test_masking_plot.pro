pro tukey_test_masking_plot, tukey=tukey

  if n_elements(xplotrange) lt 1 then xplotrange = [1e5, 1e30]
  if n_elements(yplotrange) lt 1 then yplotrange = [1e5, 1e30]
  
  obsid1 = '1061316176'
  obsid2 = '1061316296'
  use_cube = 'res'
  choose_terms = [1]
  polarizations = ['xx']
  data_range = [1e-6,1e-2]
  cut_val = 0.1
  
  if keyword_set(tukey) then begin
    note_part_1 = use_cube+' cubes, obsid '+obsid1+', tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
    note_part_2 = use_cube+' cubes, obsid '+obsid2+', tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endif else begin
    note_part_1 = use_cube+' cubes, obsid '+obsid1+', no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
    note_part_2 = use_cube+' cubes, obsid '+obsid2+', no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endelse
  
  for term = 0, n_elements(choose_terms) - 1 do begin
    for pol = 0, n_elements(polarizations) - 1 do begin
      
      filepath = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/ps/'
      if keyword_set(tukey) then begin
        filename1 = filepath + '1061316176_cubeXX__even_odd_joint_tk_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_tk_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endif else begin
        filename1 = filepath + 'Combined_obs_1061316176_cubeXX__even_odd_joint_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
        filename2 = filepath + 'Combined_obs_1061316296_cubeXX__even_odd_joint_'+use_cube+'_' + polarizations[pol] + '_averemove_bh_kcube.idlsave'
      endelse
      
      check_file1 = file_test(filename1)
      check_file2 = file_test(filename2)
      if check_file1 eq 0 or check_file2 eq 0 then begin
        print, '***ERROR*** : files not found'
      endif
      
      wt_meas_min_1 = getvar_savefile(filename1, "wt_meas_min")
      wt_meas_min_2 = getvar_savefile(filename2, "wt_meas_min")
      
      kx_mpc = getvar_savefile(filename1, 'KX_MPC')
      ky_mpc = getvar_savefile(filename1, 'KY_MPC')
      
      mask_1 = wt_meas_min_1
      mask_1[*,*] = 1
      mask_2 = wt_meas_min_2
      mask_2[*,*] = 1
      cut_1 = WHERE(wt_meas_min_1 LT cut_val, count)
      IF count GT 0 THEN mask_1[cut_1] = 0
      cut_2 = WHERE(wt_meas_min_2 LT cut_val, count)
      IF count GT 0 THEN mask_2[cut_2] = 0
      
      if keyword_set(tukey) then begin
        output_loc_1 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_'+obsid1+'_tukey_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
        output_loc_2 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_'+obsid2+'_tukey_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endif else begin
        output_loc_1 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_'+obsid1+'_tophat_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
        output_loc_2 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_'+obsid2+'_tophat_'+use_cube+'_'+polarizations[pol]+'_term' + STRTRIM(STRING(choose_terms[term]), 2)
      endelse
      
      QUICK_IMAGE, mask_1, kx_mpc, ky_mpc, $
        DATA_ASPECT = .5,$
        TITLE = 'Masking', $
        XTITLE = 'kx', YTITLE = 'ky',$
        note = note_part_1, $
        data_range = [0,1], savefile = output_loc_1
        
      QUICK_IMAGE, mask_2, kx_mpc, ky_mpc, $
        DATA_ASPECT = .5,$
        TITLE = 'Masking', $
        XTITLE = 'kx', YTITLE = 'ky',$
        note = note_part_2, $
        data_range = [0,1], savefile = output_loc_2
        
    endfor
  endfor
end