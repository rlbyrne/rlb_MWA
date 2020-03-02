pro pol_leakage_calc, $
  plot=plot, $ ; produce plots of the polarization leakage fits
  make_catalog=make_catalog, $ ; produce catalogs that correct for polarization leakage
  write_fit_params=write_fit_params, $ ; record the leakage fit parameters in a csv
  create_decon_catalogs=create_decon_catalogs, $ ; produce catalogs from the FHD deconvolution outputs (must be set if those catalogs do not exist)
  image_units_jy_per_sr=image_units_jy_per_sr, $ ; set this keyword if the FHD outputs are in units Jy/sr (i.e. created after Feb 2020)
  obs_list=obs_list ; optional list of obsids or path to a list of obsids. If unset, code will use all available observations.
  
  fhd_path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_Aug2019/'
  save_path = '/Users/rubybyrne/polarization_leakage/pol_leakage_Feb2020_new/'
  
  ; Find obsids where all files are present
  if ~keyword_set(obs_list) then begin
    decon_filenames = file_search(fhd_path+'deconvolution/*') 
    obsids_decon = []
    for file=0,n_elements(decon_filenames)-1 do begin
      obsids_decon = [obsids_decon, strmid(decon_filenames, strlen(fhd_path+'deconvolution/'), 10)]
    endfor
  endif else begin
    if isa(obs_list, /array) then begin ;obs_list is an array
      obsids_decon = obs_list
    endif else begin ;else assume obs_list is a path to a .txt file of obsids
      nlines = FILE_LINES(obs_list)
      obsids_decon = STRARR(nlines)
      OPENR, file, obs_list, /GET_LUN
      READF, file, obsids_decon
      FREE_LUN, file
    endelse
  endelse
  
  residual_I_filenames = file_search(fhd_path+'output_data/*_uniform_Residual_I.fits')
  obsids_residual_I = []
  for file=0,n_elements(residual_I_filenames)-1 do begin
    obsids_residual_I = [obsids_residual_I, strmid(residual_I_filenames, strlen(fhd_path+'output_data/'), 10)]
  endfor
  residual_Q_filenames = file_search(fhd_path+'output_data/*_uniform_Residual_Q.fits')
  obsids_residual_Q = []
  for file=0,n_elements(residual_Q_filenames)-1 do begin
    obsids_residual_Q = [obsids_residual_Q, strmid(residual_Q_filenames, strlen(fhd_path+'output_data/'), 10)]
  endfor
  residual_U_filenames = file_search(fhd_path+'output_data/*_uniform_Residual_U.fits')
  obsids_residual_U = []
  for file=0,n_elements(residual_U_filenames)-1 do begin
    obsids_residual_U = [obsids_residual_U, strmid(residual_U_filenames, strlen(fhd_path+'output_data/'), 10)]
  endfor
  residual_V_filenames = file_search(fhd_path+'output_data/*_uniform_Residual_V.fits')
  obsids_residual_V = []
  for file=0,n_elements(residual_V_filenames)-1 do begin
    obsids_residual_V = [obsids_residual_V, strmid(residual_V_filenames, strlen(fhd_path+'output_data/'), 10)]
  endfor
  obs_filenames = file_search(fhd_path+'metadata/*_obs.sav')
  obsids_obsfile = []
  for file=0,n_elements(obs_filenames)-1 do begin
    obsids_obsfile = [obsids_obsfile, strmid(obs_filenames, strlen(fhd_path+'metadata/'), 10)]
  endfor
  obsids = cgsetintersection(cgsetintersection(cgsetintersection(cgsetintersection(cgsetintersection(long(obsids_decon), long(obsids_residual_I)), long(obsids_residual_Q)), long(obsids_residual_U)), long(obsids_residual_V)), long(obsids_obsfile))
  obsids = strtrim(obsids[uniq(obsids, sort(obsids))],1)
  print, 'Processing '+string(n_elements(obsids))+' observations.'
  
  if keyword_set(create_decon_catalogs) then begin
    for obs_index=0, n_elements(obsids)-1 do begin
      obsid = obsids[obs_index]
      print, 'Generating a catalog from deconvolution outputs, saving to '+fhd_path+'decon_catalogs/'+obsid+'_decon_catalog.sav'
      generate_calibration_catalog, source_list=fhd_path+'deconvolution/'+obsid+'_fhd.sav', file_path=fhd_path+'decon_catalogs/'+obsid+'_decon_catalog.sav'
    endfor
  endif  
  
  if keyword_set(write_fit_params) then begin
    spawn, 'echo "obsid,Q leakage param1,Q leakage param2,Q leakage param3,Q leakage param4,Q leakage param5,Q leakage param6,U leakage param1,U leakage param2,U leakage param3,U leakage param4,U leakage param5,U leakage param6" > '+fhd_path+'pol_leakage_fit_params.csv'
  endif
    
  for obs_index=0, n_elements(obsids)-1 do begin
    obsid = obsids[obs_index]
    deconvolution_catalog = fhd_path+'decon_catalogs/'+obsid+'_decon_catalog.sav'
    fit_sources_number = 2000
    search_radius = 25. ; use sources within this radius of the observation center in degrees
    use_extended = 1  ; if set, use extended sources
    source_size = 1.  ; stddev of the Gaussian sources fit in pixels
    isolated_source_radius = .02  ; distance between sources is at least this in degrees
    isolated_source_flux_fraction = .9  ; if sources aren't isolated, they must contain this fraction of the flux in the region
    output_path = fhd_path+'plots/'
    if n_elements(image_units_jy_per_sr) eq 0 then image_units_jy_per_sr=1
    
    ; grab sources and sort by flux
    catalog = getvar_savefile(deconvolution_catalog, 'catalog', /compatibility_mode)
    obs = getvar_savefile(fhd_path+'metadata/'+obsid+'_obs.sav', 'obs', /compatibility_mode)
    ;source_fluxes = catalog.flux.XX+catalog.flux.YY ;things work better when you use true, not apparent, fluxes
    source_fluxes = catalog.flux.I
    brightest_indices = reverse(sort(source_fluxes))
    
    ; remove extended sources if use_extended is not set
    if ~keyword_set(use_extended) then begin
      extended_indices = []
      for source_ind = 0, n_elements(brightest_indices)-1 do begin
        if catalog[brightest_indices[source_ind]].extend ne !null then extended_indices=[extended_indices, source_ind]
      endfor
      remove, extended_indices, brightest_indices
    endif
    
    ; remove sources far from the obs center
    far_indices = []
    for source_ind = 0, n_elements(brightest_indices)-1 do begin
      if (catalog[brightest_indices[source_ind]].X-obs.obsx)^2.+(catalog[brightest_indices[source_ind]].Y-obs.obsy)^2. gt (search_radius/obs.degpix)^2. then begin
        far_indices = [far_indices, source_ind]
      endif
    endfor
    remove, far_indices, brightest_indices
    
    ; remove sources close to another bright source
    close_sources_indices = []
    for source_ind_1 = 0, n_elements(brightest_indices)-1 do begin
      close_flux = 0.
      for source_ind_2 = 0, n_elements(brightest_indices)-1 do begin
        if source_ind_2 ne source_ind_1 then begin
          if (catalog[brightest_indices[source_ind_1]].X-catalog[brightest_indices[source_ind_2]].X)^2.+(catalog[brightest_indices[source_ind_1]].Y-catalog[brightest_indices[source_ind_2]].Y)^2. lt (isolated_source_radius/obs.degpix)^2. then begin
            close_flux += catalog[brightest_indices[source_ind_2]].flux.I
          endif
        endif
      endfor
      if close_flux gt (1.-isolated_source_flux_fraction)*catalog[brightest_indices[source_ind_1]].flux.I then close_sources_indices = [close_sources_indices, source_ind_1]
    endfor
    remove, close_sources_indices, brightest_indices
    
    if fit_sources_number lt n_elements(brightest_indices) then fit_sources = catalog[brightest_indices[0:fit_sources_number-1]] else begin
      fit_sources = catalog[brightest_indices]
      fit_sources_number = n_elements(brightest_indices)
    endelse
    
    residual_I = readfits(fhd_path + 'output_data/' + obsid + '_uniform_Residual_I.fits', header)  ; assume all headers are the same
    residual_Q = readfits(fhd_path + 'output_data/' + obsid + '_uniform_Residual_Q.fits')
    residual_U = readfits(fhd_path + 'output_data/' + obsid + '_uniform_Residual_U.fits')
    residual_V = readfits(fhd_path + 'output_data/' + obsid + '_uniform_Residual_V.fits')
    if keyword_set(image_units_jy_per_sr) then begin
      image_norm = (obs.degpix*!DtoR)^2.
    endif else image_norm = 1.
    
    pol_fluxes = make_array(fit_sources_number, 4, /float, value=0.)
    source_ras = make_array(fit_sources_number, 1, /float, value=0.)
    source_decs = make_array(fit_sources_number, 1, /float, value=0.)
    source_xvals = make_array(fit_sources_number, 1, /float, value=0.)
    source_yvals = make_array(fit_sources_number, 1, /float, value=0.)
    for source_ind = 0, fit_sources_number-1 do begin
      source_flux_I = fit_sources[source_ind].flux.I
      source_flux_Q = 0.
      source_flux_U = 0.
      source_flux_V = 0.
      image_mask = make_array((size(residual_I))[1], (size(residual_I))[2], /float, value=0.)
      
      if fit_sources[source_ind].extend eq !null then begin
        xvals = [fit_sources[source_ind].X]
        yvals = [fit_sources[source_ind].Y]
        n_comps = 1
      endif else begin
        xvals = (*fit_sources[source_ind].extend).X
        yvals = (*fit_sources[source_ind].extend).Y
        n_comps = n_elements(xvals)
      endelse
      
      for comp_ind = 0, n_comps-1 do begin
        xval = xvals[comp_ind]
        yval = yvals[comp_ind]
        for pixel_x = floor(-source_size*2.)-1, ceil(source_size*2.)+1 do begin
          for pixel_y = floor(-source_size*2.)-1, ceil(source_size*2.)+1 do begin
            if (round(xval)+pixel_x-xval)^2.+(round(yval)+pixel_y-yval)^2. lt (source_size*3.)^2. then begin  ; use pixels within 3*source_size of the source center
              image_mask[round(xval)+pixel_x, round(yval)+pixel_y] += exp(-((round(xval)+pixel_x-xval)^2.+(round(yval)+pixel_y-yval)^2.)/(source_size*2.))
            endif
          endfor
        endfor
      endfor
      image_mask[where(image_mask gt 1.)] = 1.
      ;source_flux_I += total(residual_I*image_mask)*image_norm ;commented out: use the deconvolved flux as reference
      source_flux_Q += total(residual_Q*image_mask)*image_norm
      source_flux_U += total(residual_U*image_mask)*image_norm
      source_flux_V += total(residual_V*image_mask)*image_norm
      
      pol_fluxes[source_ind, 0] = source_flux_I
      pol_fluxes[source_ind, 1] = source_flux_Q
      pol_fluxes[source_ind, 2] = source_flux_U
      pol_fluxes[source_ind, 3] = source_flux_V
      source_ras[source_ind] = fit_sources[source_ind].RA
      source_decs[source_ind] = fit_sources[source_ind].DEC
      source_xvals[source_ind] = fit_sources[source_ind].X
      source_yvals[source_ind] = fit_sources[source_ind].Y
    endfor
    
    frac_pol_leakage = make_array(fit_sources_number, 3, /float, value=0.)
    frac_pol_leakage[*, 0] = pol_fluxes[*, 1] / pol_fluxes[*, 0]
    frac_pol_leakage[*, 1] = pol_fluxes[*, 2] / pol_fluxes[*, 0]
    frac_pol_leakage[*, 2] = pol_fluxes[*, 3] / pol_fluxes[*, 0]
    weights = pol_fluxes[*, 0] ; Weight by the source flux
      
    fit_matrix = make_array(6, 6, /float, value=0.)
    fit_matrix[0,*] = [total(weights*source_xvals^4.), total(weights*source_xvals^2.*source_yvals^2.), total(weights*source_xvals^3.*source_yvals), $
      total(weights*source_xvals^3.), total(weights*source_xvals^2.*source_yvals), total(weights*source_xvals^2.)]
    fit_matrix[1,*] = [total(weights*source_xvals^2.*source_yvals^2.), total(weights*source_yvals^4.), total(weights*source_xvals*source_yvals^3.), $
      total(weights*source_xvals*source_yvals^2.), total(weights*source_yvals^3.), total(weights*source_yvals^2.)]
    fit_matrix[2,*] = [total(weights*source_xvals^3.*source_yvals), total(weights*source_xvals*source_yvals^3.), total(weights*source_xvals^2.*source_yvals^2.), $
      total(weights*source_xvals^2.*source_yvals), total(weights*source_xvals*source_yvals^2.), total(weights*source_xvals*source_yvals)]
    fit_matrix[3,*] = [total(weights*source_xvals^3.), total(weights*source_xvals*source_yvals^2.), total(weights*source_xvals^2.*source_yvals), $
      total(weights*source_xvals^2.), total(weights*source_xvals*source_yvals), total(weights*source_xvals)]
    fit_matrix[4,*] = [total(weights*source_xvals^2.*source_yvals), total(weights*source_yvals^3.), total(weights*source_xvals*source_yvals^2.), $
      total(weights*source_xvals*source_yvals), total(weights*source_yvals^2.), total(weights*source_yvals)]
    fit_matrix[5,*] = [total(weights*source_xvals^2.), total(weights*source_yvals^2.), total(weights*source_xvals*source_yvals), $
      total(weights*source_xvals), total(weights*source_yvals), total(weights)]
  
    polarizations = ['Q', 'U', 'V']
    
    fit_params = make_array(6, 2, /float, value=0.)
    for pol = 0, 1 do begin
      fit_vector = [total(weights*source_xvals^2.*frac_pol_leakage[*,pol]), total(weights*source_yvals^2.*frac_pol_leakage[*,pol]), total(weights*source_xvals*source_yvals*frac_pol_leakage[*,pol]), $
        total(weights*source_xvals*frac_pol_leakage[*,pol]), total(weights*source_yvals*frac_pol_leakage[*,pol]), total(weights*frac_pol_leakage[*,pol])]
      fit_result = matrix_multiply(invert(fit_matrix), fit_vector)
      fit_result = reform(fit_result)
      fit_params[*,pol] = fit_result  
    endfor
    
    if keyword_set(plot) then begin
      for pol = 0, 2 do begin
        print, 'Saving plot to '+save_path+obsid+'_stokes_'+polarizations[pol]+'_source_leakage.png'
        cgPS_Open, save_path+obsid+'_stokes_'+polarizations[pol]+'_source_leakage.png'
        ;colorbar_extent = max(abs(frac_pol_leakage[*,pol]))  ; autoscaling color bar
        colorbar_extent = .3
        plot_range = [1024-400, 1024+400]
        colorbar_scale_factor = 127./colorbar_extent
        colors = reform(frac_pol_leakage[*,pol])*colorbar_scale_factor+127.
        for color_ind = 0,n_elements(colors)-1 do begin
          if colors[color_ind] gt 255 then colors[color_ind]=255 else begin
            if colors[color_ind] lt 0 then colors[color_ind]=0
          endelse
        endfor
        cgLoadCT, 70
        cgplot, source_xvals, source_yvals, /nodata, xtitle='RA (pixel)', ytitle='Dec (pixel)', title='Stokes '+polarizations[pol]+' Leakage', $
          aspect=1.0, xrange=[plot_range[0], plot_range[1]], yrange=[plot_range[0], plot_range[1]], position=pos
        if pol ne 2 then begin
          resolution=200
          pixel_size = (plot_range[1]-plot_range[0])/float(200)
          image_xvals = findgen(resolution, start=plot_range[0]+pixel_size/2., increment=pixel_size)
          image_yvals = findgen(resolution, start=plot_range[0]+pixel_size/2., increment=pixel_size)
          fit_image = make_array(resolution, resolution, /float, value=0.)
          for xind = 0, resolution-1 do begin
            for yind = 0, resolution-1 do begin
              fit_val = fit_params[0,pol]*image_xvals[xind]^2. + fit_params[1,pol]*image_yvals[yind]^2. + fit_params[2,pol]*image_xvals[xind]*image_yvals[yind] $
                + fit_params[3,pol]*image_xvals[xind] + fit_params[4,pol]*image_yvals[yind] + fit_params[5,pol]
              fit_color = fit_val*colorbar_scale_factor+127.
              if fit_color gt 255 then fit_color=255 else begin
                if fit_color lt 0 then fit_color=0
              endelse
              fit_image[xind, yind] = fit_color
            endfor
          endfor
          cgimage, fit_image, position=pos, /noerase
        endif
        for j=0,fit_sources_number-1 do cgplots, source_xvals[j], source_yvals[j], psym=16, color=fix(colors[j])
        cgplot, source_xvals, source_yvals, psym=9, color='black', thick=0.1, /overplot
        cgcolorbar, range=[-colorbar_extent, colorbar_extent], /vertical, $
          title = 'Fractional Polarization Leakage'
        cgps_close, /delete_ps, width=1200  ; width determines the image resolution
      endfor
    endif
      
    if keyword_set(make_catalog) then begin
      leakage_threshold = .4
      for source_ind = 0, n_elements(catalog)-1 do begin
        q_leakage = fit_params[0,0]*(catalog[source_ind].X)^2. + fit_params[1,0]*(catalog[source_ind].Y)^2. + fit_params[2,0]*(catalog[source_ind].X)*(catalog[source_ind].Y) $
          + fit_params[3,0]*(catalog[source_ind].X) + fit_params[4,0]*(catalog[source_ind].Y) + fit_params[5,0]
        if q_leakage lt -leakage_threshold then q_leakage=-leakage_threshold
        if q_leakage gt leakage_threshold then q_leakage=leakage_threshold
        u_leakage = fit_params[0,1]*(catalog[source_ind].X)^2. + fit_params[1,1]*(catalog[source_ind].Y)^2. + fit_params[2,1]*(catalog[source_ind].X)*(catalog[source_ind].Y) $
          + fit_params[3,1]*(catalog[source_ind].X) + fit_params[4,1]*(catalog[source_ind].Y) + fit_params[5,1]
        if u_leakage lt -leakage_threshold then u_leakage=-leakage_threshold
        if u_leakage gt leakage_threshold then u_leakage=leakage_threshold
        catalog[source_ind].flux.Q = catalog[source_ind].flux.I * q_leakage
        catalog[source_ind].flux.U = catalog[source_ind].flux.I * u_leakage
        if catalog[source_ind].extend ne !null then begin
          for comp_ind = 0, n_elements(*catalog[source_ind].extend)-1 do begin
            (*catalog[source_ind].extend)[comp_ind].flux.Q = (*catalog[source_ind].extend)[comp_ind].flux.I * q_leakage
            (*catalog[source_ind].extend)[comp_ind].flux.U = (*catalog[source_ind].extend)[comp_ind].flux.I * u_leakage
          endfor
        endif
      endfor
      print, 'Saving polarization leakage corrected catalog to '+save_path+obsid+'_decon_catalog_pol_leakage_corrected.sav'
      save, catalog, filename=save_path+obsid+'_decon_catalog_pol_leakage_corrected.sav'
    endif
    
    if keyword_set(write_fit_params) then begin
      line = obsid+','+strtrim(fit_params[0,0],2)+','+strtrim(fit_params[1,0],2)+','+strtrim(fit_params[2,0],2)+','+strtrim(fit_params[3,0],2)+','+strtrim(fit_params[4,0],2)+','+strtrim(fit_params[5,0],2) $
        +','+strtrim(fit_params[0,1],2)+','+strtrim(fit_params[1,1],2)+','+strtrim(fit_params[2,1],2)+','+strtrim(fit_params[3,1],2)+','+strtrim(fit_params[4,1],2)+','+strtrim(fit_params[5,1],2)
      print, 'Saving fit parameters to '+save_path+'pol_leakage_fit_params.csv'
      spawn, 'printf "'+line+'\n" >> '+save_path+'pol_leakage_fit_params.csv'
    endif

  endfor
  
end