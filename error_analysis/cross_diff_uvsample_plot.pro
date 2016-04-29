pro cross_diff_uvsample_plot

  uv_img_clip = 6
  uv_img_clip = number_formatter(round(uv_img_clip))
  sample_factors = [.0002,.0005,.001,.005,.01,.05,.1,.5,1,5]
  sample_factors = number_formatter(sample_factors)
  
  for sample = 0, n_elements(sample_factors)-1 do begin
  
    pos = strpos(sample_factors[sample],'.')
    if pos ne -1 then filename_sample = strmid(sample_factors[sample],0,pos) + 'p' + $
      strmid(sample_factors[sample], pos+1, strlen(sample_factors[sample])-pos) else $
      filename_sample = sample_factors[sample]
      
    kx_mpc = getvar_savefile('/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/cross_diff_noise_sim_flatUV'+ $
      filename_sample+'_uvsample'+uv_img_clip+'_xx_term1.sav', 'kx_mpc')
    ky_mpc = getvar_savefile('/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/cross_diff_noise_sim_flatUV'+ $
      filename_sample+'_uvsample'+uv_img_clip+'_xx_term1.sav', 'ky_mpc')
    cubecross_sample = getvar_savefile('/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/cross_diff_noise_sim_flatUV'+ $
      filename_sample+'_uvsample'+uv_img_clip+'_xx_term1.sav', 'cube_cross')
    cubecross_clip = getvar_savefile('/data3/MWA/FHD_Aug23/fhd_rlb_cross_diff_noise_sim_flatUV/cross_diff_noise_sim_flatUV'+ $
      filename_sample+'_uvimgclip'+uv_img_clip+'_xx_term1.sav', 'cube_cross')
    varmeas_sample = variance(real_part(cubecross_sample), dimension=3)
    varmeas_clip = variance(real_part(cubecross_clip), dimension=3)
    varmeas_ratio = varmeas_sample / varmeas_clip
    varmeas_ratio2 = sqrt(varmeas_sample) / varmeas_clip
    where0 = where(varmeas_clip eq 0, count)
    if count gt 0 then begin
      varmeas_ratio[where0]=0
      varmeas_ratio2[where0]=0
      endif
    varmeas_diff = (varmeas_clip - varmeas_sample)/(varmeas_clip + varmeas_sample)
    where0 = where((varmeas_clip+varmeas_sample) eq 0, count)
    if count gt 0 then varmeas_diff[where0]=0
    note = 'crossed sim noise with ' + sample_factors[sample] + ' UV coverage, single obs, UV img clip ' + uv_img_clip
    output_loc = '/home/rlbyrne/error_analysis_plots/uvsample_plot_UVsim'+filename_sample+'_uvimgclip'+uv_img_clip+'_xx_term1'
    ;CGPS_OPEN, output_loc, /FONT, XSIZE = 30, YSIZE = 15
    CGPS_OPEN, output_loc, /FONT, XSIZE = 13, YSIZE = 6
        
    data_range = [0,1e26]
    ;QUICK_IMAGE, varmeas_sample, kx_mpc, ky_mpc, DATA_RANGE = data_range, $
    ;  TITLE = 'Meas. Var. Sampled', XTITLE = 'kx', YTITLE = 'ky', $
    ;  /log, start_multi_params = {nrow:2, ncol:2}, multi_pos = multi_pos, $
    ;  /noerase, no_ps_close=savefile
    ;QUICK_IMAGE, varmeas_clip, kx_mpc, ky_mpc, DATA_RANGE = data_range, $
    ;  TITLE = 'Meas. Var. Img. Clipped', XTITLE = 'kx', YTITLE = 'ky', $
    ;  /log, multi_pos = multi_pos[*,1], /noerase
    ;QUICK_IMAGE, varmeas_ratio, kx_mpc, ky_mpc, DATA_RANGE = [0,2], $
    ;  TITLE = 'Sampled/Clipped', XTITLE = 'kx', YTITLE = 'ky', $
    ;  multi_pos = multi_pos[*,2], /noerase
    ;QUICK_IMAGE, varmeas_diff, kx_mpc, ky_mpc, DATA_RANGE = [0,1], $
    ;  TITLE = '(C-S)/(C+S)', XTITLE = 'kx', YTITLE = 'ky', $
    ;  NOTE = note, $
    ;  multi_pos = multi_pos[*,3], no_ps_close=savefile, /noerase
    quick_image, varmeas_ratio2, kx_mpc, ky_mpc, data_range = sqrt(data_range), $
      TITLE = 'Sampled^2/Clipped', XTITLE = 'kx', YTITLE = 'ky', note=note, no_ps_close=savefile, $
      data_aspect = 1, /log
      
    undefine, multi_pos
    CGPS_CLOSE, /PNG, /DELETE_PS
    
  endfor
  
end