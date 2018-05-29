pro pol_leakage_calc

  obsid = '1131478776'
  deconvolution_catalog = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018/1131478776_decon_catalog.sav'
  fit_sources_number = 500
  fit_radius = 15  ; use only sources within 15 degrees of the pointing center
  output_path = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018/'
  
  catalog = getvar_savefile(deconvolution_catalog, 'catalog', /compatibility_mode)
  obs = getvar_savefile(output_path+'metadata/'+obsid+'_obs.sav', 'obs', /compatibility_mode)
  apparent_fluxes = make_array(n_elements(catalog), /float, value=0.)
  for source_ind = 0, n_elements(catalog)-1 do begin
    if (catalog[source_ind].ra-obs.obsra)^2. + (catalog[source_ind].dec-obs.obsdec)^2. lt fit_radius^2. then begin
   ; if catalog[source_ind].extend eq !null then begin  ; use point sources only (no extended)
      apparent_fluxes[source_ind] = catalog[source_ind].flux.I
   ; endif
    endif
  endfor
  brightest_indices = (reverse(sort(apparent_fluxes)))[0:fit_sources_number-1]
  fit_sources = catalog[brightest_indices]
  
  residual_I = readfits(output_path + 'output_data/' + obsid + '_uniform_Residual_I.fits', header)  ; assume all headers are the same
  residual_Q = readfits(output_path + 'output_data/' + obsid + '_uniform_Residual_Q.fits')
  residual_U = readfits(output_path + 'output_data/' + obsid + '_uniform_Residual_U.fits')
  residual_V = readfits(output_path + 'output_data/' + obsid + '_uniform_Residual_V.fits')
  
  pol_fluxes = make_array(fit_sources_number, 4, /float, value=0.)
  source_ras = make_array(fit_sources_number, 1, /float, value=0.)
  source_decs = make_array(fit_sources_number, 1, /float, value=0.)
  source_xvals = make_array(fit_sources_number, 1, /float, value=0.)
  source_yvals = make_array(fit_sources_number, 1, /float, value=0.)
  for source_ind = 0, fit_sources_number-1 do begin
    pol_fluxes[source_ind, 0] = fit_sources[source_ind].flux.I + residual_I[round(fit_sources[source_ind].X), round(fit_sources[source_ind].Y)]
    pol_fluxes[source_ind, 1] = residual_Q[round(fit_sources[source_ind].X), round(fit_sources[source_ind].Y)]
    pol_fluxes[source_ind, 2] = residual_U[round(fit_sources[source_ind].X), round(fit_sources[source_ind].Y)]
    pol_fluxes[source_ind, 3] = residual_V[round(fit_sources[source_ind].X), round(fit_sources[source_ind].Y)]
    source_ras[source_ind] = fit_sources[source_ind].RA
    source_decs[source_ind] = fit_sources[source_ind].DEC
    source_xvals[source_ind] = fit_sources[source_ind].X
    source_yvals[source_ind] = fit_sources[source_ind].Y
  endfor
  
  frac_pol_leakage = make_array(fit_sources_number, 3, /float, value=0.)
  frac_pol_leakage[*, 0] = pol_fluxes[*, 1] / pol_fluxes[*, 0]
  frac_pol_leakage[*, 1] = pol_fluxes[*, 2] / pol_fluxes[*, 0]
  frac_pol_leakage[*, 2] = pol_fluxes[*, 3] / pol_fluxes[*, 0]
    
  fit_matrix = make_array(6, 6, /float, value=0.)
  fit_matrix[0,*] = [total(source_xvals^4.), total(source_xvals^2.*source_yvals^2.), total(source_xvals^3.*source_yvals), $
    total(source_xvals^3.), total(source_xvals^2.*source_yvals), total(source_xvals^2.)]
  fit_matrix[1,*] = [total(source_xvals^2.*source_yvals^2.), total(source_yvals^4.), total(source_xvals*source_yvals^3.), $
    total(source_xvals*source_yvals^2.), total(source_yvals^3.), total(source_yvals^2.)]
  fit_matrix[2,*] = [total(source_xvals^3.*source_yvals), total(source_xvals*source_yvals^3.), total(source_xvals^2.*source_yvals^2.), $
    total(source_xvals^2.*source_yvals), total(source_xvals*source_yvals^2.), total(source_xvals*source_yvals)]
  fit_matrix[3,*] = [total(source_xvals^3.), total(source_xvals*source_yvals^2.), total(source_xvals^2.*source_yvals), $
    total(source_xvals^2.), total(source_xvals*source_yvals), total(source_xvals)]
  fit_matrix[4,*] = [total(source_xvals^2.*source_yvals), total(source_yvals^3.), total(source_xvals*source_yvals^2.), $
    total(source_xvals*source_yvals), total(source_yvals^2.), total(source_yvals)]
  fit_matrix[5,*] = [total(source_xvals^2.), total(source_yvals^2.), total(source_xvals*source_yvals), $
    total(source_xvals), total(source_yvals), fit_sources_number]
    
  fit_vector = [total(source_xvals^2.*frac_pol_leakage[*,1]), total(source_yvals^2.*frac_pol_leakage[*,1]), total(source_xvals*source_yvals*frac_pol_leakage[*,1]), $
    total(source_xvals*frac_pol_leakage[*,1]), total(source_yvals*frac_pol_leakage[*,1]), total(frac_pol_leakage[*,1])]
  
  fit_result = matrix_multiply(invert(fit_matrix), fit_vector)
  fit_result = reform(fit_result)
  print, fit_result
  
  polarizations = ['Q', 'U', 'V']
  for pol = 0, 2 do begin
    cgPS_Open, obsid+'_stokes_'+polarizations[pol]+'_source_leakage.png'
    colorbar_scale_factor = 127./max(abs(frac_pol_leakage[*,pol]))
    colors = reform(frac_pol_leakage[*,pol])*colorbar_scale_factor+127.
    cgLoadCT, 70
    cgplot, source_ras, source_decs, /nodata, xtitle='RA (deg)', ytitle='Dec (deg)', title='Stokes '+polarizations[pol]+' Leakage', $
      aspect=1.0, xrange=[obs.obsra-fit_radius-5, obs.obsra+fit_radius+5], yrange=[obs.obsdec-fit_radius-5, obs.obsdec+fit_radius+5]
    for j=0,fit_sources_number-1 do cgplots, source_ras[j], source_decs[j], psym=16, color=fix(colors[j])
    cgcolorbar, range=[-max(abs(frac_pol_leakage[*,pol])), max(abs(frac_pol_leakage[*,pol]))], /vertical, $
      title = 'Fractional Polarization Leakage'
    cgps_close, /delete_ps
  endfor
end