pro weights_plotting, after_eppsilon=after_eppsilon, pol=pol
  
  if keyword_set(after_eppsilon) then begin
    
    fhd_path = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/'
    fhd_run = 'fhd_rlb_array_sim_Barry_effect_Jun2018'
    ;obs_name = 'sidelobe_survey_obsIDs_firstobs'
    obsid = 'hex_array_sim_331'

    if ~keyword_set(pol) then pol = 'XX'
    if pol ne 'XX' and pol ne 'YY' then pol = 'XX'
  
    ;weights_odd = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uvf.idlsave', 'weights_cube')
    ;weights_even = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_even_cube' + pol + '_weights_uvf.idlsave', 'weights_cube')
    weights_odd = getvar_savefile(fhd_path + fhd_run + '/ps/data/uvf_cubes/'+obsid+'_odd_cube'+pol+'_weights_uvf.idlsave', 'weights_cube')
    weights_even = getvar_savefile(fhd_path + fhd_run + '/ps/data/uvf_cubes/'+obsid+'_even_cube'+pol+'_weights_uvf.idlsave', 'weights_cube')
    weights_evenodd = weights_odd + weights_even
    weights_cube = total(abs(weights_evenodd), 3)
    weights_cube_squared = total(abs(weights_evenodd^2.), 3)
    
    ;xarr = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'xarr')
    ;yarr = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'yarr')
    ;kperp_lambda_conv = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'kperp_lambda_conv')
    ;xarr = kperp_lambda_conv * xarr
    ;yarr = kperp_lambda_conv * yarr
    xarr = getvar_savefile(fhd_path + fhd_run + '/ps/data/uvf_cubes/'+obsid+'_odd_cube'+pol+'_weights_uvf.idlsave', 'kx_rad_vals')
    yarr = getvar_savefile(fhd_path + fhd_run + '/ps/data/uvf_cubes/'+obsid+'_odd_cube'+pol+'_weights_uvf.idlsave', 'ky_rad_vals')
    
  endif else begin
  
    weights_cube = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018/grid_data/hex_array_sim_331_uv_weights_XX.sav', 'weights_uv')
    xarr = (findgen(2048)-1024)*.5
    yarr = (findgen(2048)-1024)*.5
    stop
  endelse
  
  quick_image, abs(weights_cube), xarr, yarr, xtitle = 'U (wavelengths)', ytitle = 'V (wavelengths)', cb_title = 'Weights (Jy/beam)', data_range = [0,2e-5], savefile='/Users/rubybyrne/weights_plots_for_cal_paper/hex_array.png'
  
  
end