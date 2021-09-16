pro export_pol_leakage_surface_as_healpix, obs_list=obs_list

  ;fhd_path = '/Users/ruby/Astro/diffuse_survey/fhd_rlb_diffuse_survey_decon_4pol_May2021/'
  ;leakage_params_csv = '/Users/ruby/Astro/diffuse_survey/pol_leakage_fit_params.csv'
  ;save_path = '/Users/ruby/Astro/diffuse_survey/pol_leakage_hpx_surfaces/'
  fhd_path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_diffuse_survey_decon_4pol_May2021/'
  save_path = '/Users/rubybyrne/polarization_leakage/pol_leakage_Jun2021/healpix_leakage_maps/'
  leakage_params_csv = '/Users/rubybyrne/polarization_leakage/pol_leakage_Jun2021/pol_leakage_fit_params.csv'
  
  ; Find obsids where all files are present
  if ~keyword_set(obs_list) then begin
    decon_filenames = file_search(fhd_path+'metadata/*')
    obsids_decon = [strmid(decon_filenames, strlen(fhd_path+'metadata/'), 10)]
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
  obsids_decon = obsids_decon[uniq(obsids_decon)] ;remove duplicates
  
  ; Get surface fit parameters
  pol_leakage_params_file = read_csv(leakage_params_csv, header=csv_header)
  pol_leakage_params = make_array(n_elements(obsids_decon), 12, /double, value=0.)
  
  for obsind=0,n_elements(obsids_decon)-1 do begin
    
    obsid = obsids_decon[obsind]

    obs = getvar_savefile(fhd_path+'metadata/'+obsid+'_obs.sav', 'obs', /compatibility_mode)
    nside = 512
    hpx_ordering = 'ring'
    coord_sys = 'equatorial'
    n_hpx = nside2npix(nside)
    hpx_inds = Lindgen(n_hpx)
    
    IF hpx_ordering EQ 'ring' THEN pix2vec_ring, nside, hpx_inds, pix_coords $
    ELSE pix2vec_nest, nside, hpx_inds, pix_coords
  
    vec2ang,pix_coords,pix_dec,pix_ra,/astro
    IF coord_sys EQ 'galactic' THEN glactc,pix_ra,pix_dec,2000.,pix_ra,pix_dec,2, /degree
    IF coord_sys EQ 'equatorial' THEN Hor2Eq,pix_dec,pix_ra,obs.JD0,pix_ra,pix_dec,lat=obs.lat,lon=obs.lon,alt=obs.alt,precess=1,/nutate
  
    apply_astrometry, obs, ra_arr=pix_ra, dec_arr=pix_dec, x_arr=xv_hpx, y_arr=yv_hpx, /ad2xy
  
    hpx_i_use=where((xv_hpx GT 0) AND (xv_hpx LT (obs.dimension-1)) AND (yv_hpx GT 0) AND (yv_hpx LT (obs.elements-1)),n_hpx_use)
    IF n_hpx_use EQ 0 THEN BEGIN
      print,"Error: Map has no valid Healpix indices"
    ENDIF
    xv_hpx=xv_hpx[hpx_i_use]
    
    params_file_obsind = where(pol_leakage_params_file.field01 eq obsid, count)
    if count ne 1 then begin
      print, 'ERROR: Obsid not present in pol leakage parameter file.'
      return
    endif
    fit_params = make_array(6, 2, /float, value=0.)
    fit_params[0,0] = pol_leakage_params_file.field02[params_file_obsind]
    fit_params[1,0] = pol_leakage_params_file.field03[params_file_obsind]
    fit_params[2,0] = pol_leakage_params_file.field04[params_file_obsind]
    fit_params[3,0] = pol_leakage_params_file.field05[params_file_obsind]
    fit_params[4,0] = pol_leakage_params_file.field06[params_file_obsind]
    fit_params[5,0] = pol_leakage_params_file.field07[params_file_obsind]
    fit_params[0,1] = pol_leakage_params_file.field08[params_file_obsind]
    fit_params[1,1] = pol_leakage_params_file.field09[params_file_obsind]
    fit_params[2,1] = pol_leakage_params_file.field10[params_file_obsind]
    fit_params[3,1] = pol_leakage_params_file.field11[params_file_obsind]
    fit_params[4,1] = pol_leakage_params_file.field12[params_file_obsind]
    fit_params[5,1] = pol_leakage_params_file.field13[params_file_obsind]
    
    q_leakage = make_array(n_elements(hpx_i_use), /float, value=0.)
    u_leakage = make_array(n_elements(hpx_i_use), /float, value=0.)
    leakage_threshold = .4
    for source_ind = 0, n_elements(hpx_i_use)-1 do begin
      q_leakage[source_ind] = fit_params[0,0]*(xv_hpx[source_ind])^2. + fit_params[1,0]*(yv_hpx[source_ind])^2. + fit_params[2,0]*(xv_hpx[source_ind])*(yv_hpx[source_ind]) $
        + fit_params[3,0]*(xv_hpx[source_ind]) + fit_params[4,0]*(yv_hpx[source_ind]) + fit_params[5,0]
      if q_leakage[source_ind] lt -leakage_threshold then q_leakage[source_ind]=-leakage_threshold
      if q_leakage[source_ind] gt leakage_threshold then q_leakage[source_ind]=leakage_threshold
      u_leakage[source_ind] = fit_params[0,1]*(xv_hpx[source_ind])^2. + fit_params[1,1]*(yv_hpx[source_ind])^2. + fit_params[2,1]*(xv_hpx[source_ind])*(yv_hpx[source_ind]) $
        + fit_params[3,1]*(xv_hpx[source_ind]) + fit_params[4,1]*(yv_hpx[source_ind]) + fit_params[5,1]
      if u_leakage[source_ind] lt -leakage_threshold then u_leakage[source_ind]=-leakage_threshold
      if u_leakage[source_ind] gt leakage_threshold then u_leakage[source_ind]=leakage_threshold
    endfor
    
    save, obsid, nside, hpx_ordering, coord_sys, hpx_i_use, q_leakage, u_leakage, filename=save_path+obsid+'_hpx_pol_leakage_surface.sav'
  
  endfor


end