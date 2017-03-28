pro weights_plotting, after_eppsilon=after_eppsilon, pol=pol

  fhd_path = '/nfs/mwa-08/d1/DiffuseSurvey2015/'
  fhd_run = 'fhd_rlb_diffuse_survey_oneobs_nodiffuse'
  obs_name = 'sidelobe_survey_obsIDs_firstobs'
  obsid = '1131454296'
  
  if ~keyword_set(pol) then pol = 'XX'
  if pol ne 'XX' and pol ne 'YY' then pol = 'XX'
  
  if keyword_set(after_eppsilon) then begin
  
    weights_odd = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uvf.idlsave', 'weights_cube')
    weights_even = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_even_cube' + pol + '_weights_uvf.idlsave', 'weights_cube')
    weights_evenodd = weights_odd + weights_even
    weights_cube = total(abs(weights_evenodd), 3)
    weights_cube_squared = total(abs(weights_evenodd^2.), 3)
    
    xarr = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'xarr')
    yarr = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'yarr')
    kperp_lambda_conv = getvar_savefile(fhd_path + fhd_run + '/ps/Combined_obs_' + obs_name + '_odd_cube' + pol + '_weights_uv_plane.idlsave', 'kperp_lambda_conv')
    xarr = kperp_lambda_conv * xarr
    yarr = kperp_lambda_conv * yarr
    
  endif else begin
  
    weights_cube = getvar_savefile(fhd_path + fhd_run + '/grid_data/' + obsid + '_uv_weights_' + pol + '.sav', 'weights_uv')
    
  endelse
  
  quick_image, weights_cube, xarr, yarr, xtitle = 'U (wavelengths)', ytitle = 'V (wavelengths)', title = 'Weights'
  
  
end