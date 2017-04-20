pro masking_gradient_plot, tukey=tukey

  obsid1 = '1061316176'
  obsid2 = '1061316296'
  use_cube = 'res'
  choose_terms = [1]
  polarizations = ['xx']
  data_range = [1e-6,1e-2]
  cut_val = 0.1
  
  if keyword_set(tukey) then begin
    note_part_1 = use_cube+' cube, obsid '+obsid1+', tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
    note_part_2 = use_cube+' cube, obsid '+obsid2+', tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
  endif else begin
    note_part_1 = use_cube+' cube, obsid '+obsid1+', no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
    note_part_2 = use_cube+' cube, obsid '+obsid2+', no tukey filter, pol ' + polarizations + ', term ' + strtrim(string(choose_terms),2)
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
      
      wt_meas_min = getvar_savefile(filename1, "wt_meas_min")
      wt_meas_ave = getvar_savefile(filename1, "wt_meas_ave")
      
      kx_mpc = getvar_savefile(filename1, 'KX_MPC')
      ky_mpc = getvar_savefile(filename1, 'KY_MPC')
      kx_step = kx_mpc[1]-kx_mpc[0]
      ky_step = ky_mpc[1]-ky_mpc[0]
      
      ;wt_meas_min_grad = make_array((size(wt_meas_min))[1]-1, (size(wt_meas_min))[2]-1, /float, value=0)
      ;denominator = make_array((size(wt_meas_min))[1]-1, (size(wt_meas_min))[2]-1, /float, value=0)
      ;for i = 0, (size(wt_meas_min))[1]-2 do begin
      ;  for j = 0, (size(wt_meas_min))[2]-2 do begin
      ;    xgrad = (wt_meas_min[i+1,j]-wt_meas_min[i,j])/kx_step
      ;    ygrad = (wt_meas_min[i,j+1]-wt_meas_min[i,j])/ky_step
      ;    grad_amp = sqrt(xgrad^2+ygrad^2)
      ;    wt_meas_min_grad[i,j] = grad_amp
      ;    denominator[i,j] = (2*wt_meas_ave[i,j]+wt_meas_ave[i+1,j]+wt_meas_ave[i,j+1])/4
      ;  endfor
      ;endfor
      
      ;wt_meas_ratio = wt_meas_min_grad/denominator
      
      wt_meas_ave_grad = make_array((size(wt_meas_ave))[1]-1, (size(wt_meas_ave))[2]-1, /float, value=0)
      denominator = make_array((size(wt_meas_ave))[1]-1, (size(wt_meas_ave))[2]-1, /float, value=0)
      for i = 0, (size(wt_meas_ave))[1]-2 do begin
        for j = 0, (size(wt_meas_ave))[2]-2 do begin
          xgrad = (wt_meas_ave[i+1,j]-wt_meas_ave[i,j])/kx_step
          ygrad = (wt_meas_ave[i,j+1]-wt_meas_ave[i,j])/ky_step
          grad_amp = sqrt(xgrad^2+ygrad^2)
          wt_meas_ave_grad[i,j] = grad_amp
          denominator[i,j] = (2*wt_meas_ave[i,j]+wt_meas_ave[i+1,j]+wt_meas_ave[i,j+1])/4
        endfor
      endfor
      
      wt_meas_ratio = wt_meas_ave_grad/denominator
      
      kx_mpc_grad = kx_mpc[0:n_elements(kx_mpc)-2]+(kx_step/2)
      ky_mpc_grad = ky_mpc[0:n_elements(ky_mpc)-2]+(ky_step/2)
      
      if keyword_set(tukey) then begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_gradient_plot_tukey
      endif else begin
        output_loc = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/masking_gradient_plot_tophat
      endelse
      
      QUICK_IMAGE, wt_meas_ratio, kx_mpc_grad, ky_mpc_grad, $
        DATA_ASPECT = .5,$
        TITLE = 'Gradient of Ave Weights/Ave Weights', $
        XTITLE = 'kx', YTITLE = 'ky',$
        note = note_part_1, $
        data_range = [0,1e3], $
        savefile = output_loc
        
        
    endfor
  endfor
end