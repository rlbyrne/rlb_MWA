pro pol_leakage_calc

  obsid = '1131478776'
  deconvolution_catalog = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018/1131478776_decon_catalog.sav'
  fit_sources_number = 400
  fit_radius = 12  ; use only sources within 15 degrees of the pointing center
  source_size = 1.5  ; stddev of the Gaussian sources fit in pixels
  output_path = '/Users/rubybyrne/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2018/'
  
  catalog = getvar_savefile(deconvolution_catalog, 'catalog', /compatibility_mode)
  obs = getvar_savefile(output_path+'metadata/'+obsid+'_obs.sav', 'obs', /compatibility_mode)
  apparent_fluxes = make_array(n_elements(catalog), /float, value=0.)
  for source_ind = 0, n_elements(catalog)-1 do begin
    if (catalog[source_ind].ra-obs.obsra)^2. + (catalog[source_ind].dec-obs.obsdec)^2. lt fit_radius^2. then begin
      if catalog[source_ind].extend eq !null then begin  ; use point sources only (no extended)
        apparent_fluxes[source_ind] = catalog[source_ind].flux.I
      endif
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
    xval = fit_sources[source_ind].X
    yval = fit_sources[source_ind].Y
    source_flux_I = fit_sources[source_ind].flux.I
    source_flux_Q = 0.
    source_flux_U = 0.
    source_flux_V = 0.
    for pixel_x = floor(-source_size*2.)-1, ceil(source_size*2.)+1 do begin
      for pixel_y = floor(-source_size*2.)-1, ceil(source_size*2.)+1 do begin
        if (round(xval)+pixel_x-xval)^2.+(round(yval)+pixel_y-yval)^2. lt (source_size*3.)^2. then begin  ; use pixels within 3*source_size of the source center
          weighting = exp(-((round(xval)+pixel_x-xval)^2.+(round(yval)+pixel_y-yval)^2.)/(source_size*2.))
          source_flux_I = source_flux_I + weighting*residual_I[round(xval)+pixel_x, round(yval)+pixel_y]
          source_flux_Q = source_flux_Q + weighting*residual_Q[round(xval)+pixel_x, round(yval)+pixel_y]
          source_flux_U = source_flux_U + weighting*residual_U[round(xval)+pixel_x, round(yval)+pixel_y]
          source_flux_V = source_flux_V + weighting*residual_V[round(xval)+pixel_x, round(yval)+pixel_y]
        endif
      endfor
    endfor
    pol_fluxes[source_ind, 0] = source_flux_I
    pol_fluxes[source_ind, 1] = source_flux_Q
    pol_fluxes[source_ind, 2] = source_flux_U
    pol_fluxes[source_ind, 3] = source_flux_V
    source_ras[source_ind] = fit_sources[source_ind].RA
    source_decs[source_ind] = fit_sources[source_ind].DEC
    source_xvals[source_ind] = xval
    source_yvals[source_ind] = yval
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

  polarizations = ['Q', 'U', 'V']
  fit_params = make_array(6, 2, /float, value=0.)
  for pol = 0, 2 do begin
    
    cgPS_Open, obsid+'_stokes_'+polarizations[pol]+'_source_leakage.png'
    colorbar_scale_factor = 127./max(abs(frac_pol_leakage[*,pol]))
    colors = reform(frac_pol_leakage[*,pol])*colorbar_scale_factor+127.
    cgLoadCT, 70
    cgplot, source_ras, source_decs, /nodata, xtitle='RA (deg)', ytitle='Dec (deg)', title='Stokes '+polarizations[pol]+' Leakage', $
      aspect=1.0, xrange=[obs.obsra-fit_radius-5, obs.obsra+fit_radius+5], yrange=[obs.obsdec-fit_radius-5, obs.obsdec+fit_radius+5]
    for j=0,fit_sources_number-1 do cgplots, source_ras[j], source_decs[j], psym=16, color=fix(colors[j])
    if pol ne 2 then begin
      fit_vector = [total(source_xvals^2.*frac_pol_leakage[*,pol]), total(source_yvals^2.*frac_pol_leakage[*,pol]), total(source_xvals*source_yvals*frac_pol_leakage[*,pol]), $
        total(source_xvals*frac_pol_leakage[*,pol]), total(source_yvals*frac_pol_leakage[*,pol]), total(frac_pol_leakage[*,pol])]

      fit_result = matrix_multiply(invert(fit_matrix), fit_vector)
      fit_result = reform(fit_result)
      fit_params[*,pol] = fit_result
      fitted_vals = fit_result[0]*source_xvals^2. + fit_result[1]*source_yvals^2. + fit_result[2]*source_xvals*source_yvals $
        + fit_result[3]*source_xvals + fit_result[4]*source_yvals + fit_result[5]
      errors = frac_pol_leakage[*,pol]-fitted_vals
      
      colors_fit = fitted_vals*colorbar_scale_factor+127.
      colors_fit[where(colors_fit lt 0, /null)] = 0  ; this is buggy
      colors_fit[where(colors_fit gt 254, /null)] = 254  ; this is buggy
      for j=0,fit_sources_number-1 do cgplots, source_ras[j], source_decs[j], psym=7, color=fix(colors_fit[j])
      cgcolorbar, range=[-max(abs(frac_pol_leakage[*,pol])), max(abs(frac_pol_leakage[*,pol]))], /vertical, $
        title = 'Fractional Polarization Leakage'
      cgps_close, /delete_ps
      
      cgPS_Open, obsid+'_stokes_'+polarizations[pol]+'_source_leakage_residual.png'
      cgplot, source_ras, source_decs, /nodata, xtitle='RA (deg)', ytitle='Dec (deg)', title='Stokes '+polarizations[pol]+' Resiudal Leakage', $
        aspect=1.0, xrange=[obs.obsra-fit_radius-5, obs.obsra+fit_radius+5], yrange=[obs.obsdec-fit_radius-5, obs.obsdec+fit_radius+5]
      colors_residual = reform(errors)*colorbar_scale_factor+127.
      for j=0,fit_sources_number-1 do cgplots, source_ras[j], source_decs[j], psym=16, color=fix(colors_residual[j])
      cgcolorbar, range=[-max(abs(frac_pol_leakage[*,pol])), max(abs(frac_pol_leakage[*,pol]))], /vertical, $
        title = 'Polarization Residual'
      cgps_close, /delete_ps
    endif else begin
    cgcolorbar, range=[-max(abs(frac_pol_leakage[*,pol])), max(abs(frac_pol_leakage[*,pol]))], /vertical, $
      title = 'Fractional Polarization Leakage'
    cgps_close, /delete_ps
    endelse
  endfor
  
  stop
  
  for source_ind = 0, n_elements(catalog)-1 do begin
    q_leakage = fit_params[0,0]*(catalog[source_ind].X)^2. + fit_params[1,0]*(catalog[source_ind].Y)^2. + fit_params[2,0]*(catalog[source_ind].X)*(catalog[source_ind].Y) $
      + fit_params[3,0]*(catalog[source_ind].X) + fit_params[4,0]*(catalog[source_ind].Y) + fit_params[5,0]
    u_leakage = fit_params[0,1]*(catalog[source_ind].X)^2. + fit_params[1,1]*(catalog[source_ind].Y)^2. + fit_params[2,1]*(catalog[source_ind].X)*(catalog[source_ind].Y) $
      + fit_params[3,1]*(catalog[source_ind].X) + fit_params[4,1]*(catalog[source_ind].Y) + fit_params[5,1]
    catalog[source_ind].flux.Q = catalog[source_ind].flux.I * q_leakage
    catalog[source_ind].flux.U = catalog[source_ind].flux.I * u_leakage
    if catalog[source_ind].extend eq !null then begin
      for comp_ind = 0, n_elements(*catalog[source_ind].extend)-1 do begin
        (*catalog[source_ind].extend)[comp_ind].flux.Q = (*catalog[source_ind].extend)[comp_ind].flux.I * q_leakage
        (*catalog[source_ind].extend)[comp_ind].flux.U = (*catalog[source_ind].extend)[comp_ind].flux.I * u_leakage
      endfor
    endif
  endfor
  
end