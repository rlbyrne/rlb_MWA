pro plot_pol_leakage, obs_list=obs_list

  fhd_path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_diffuse_survey_decon_4pol_May2021/'
  save_path = '/Users/rubybyrne/polarization_leakage/pol_leakage_Jun2021/'
  
  ; Find obsids where all files are present
  if ~keyword_set(obs_list) then begin
    decon_filenames = file_search(fhd_path+'deconvolution/*')
    obsids_decon = [strmid(decon_filenames, strlen(fhd_path+'deconvolution/'), 10)]
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
  
  for obsind=0,n_elements(obsids_decon)-1 do begin
    
    obsid = obsids_decon[obsind]

    obs = getvar_savefile(fhd_path+'metadata/'+obsid+'_obs.sav', 'obs', /compatibility_mode)
    nside = 512
    hpx_ordering = 'ring'
    coord_sys = 'equatorial'
    
    IF hpx_ordering EQ 'ring' THEN pix2vec_ring, nside, hpx_inds, pix_coords $
    ELSE pix2vec_nest, nside, hpx_inds, pix_coords
  
    vec2ang,pix_coords,pix_dec,pix_ra,/astro
    IF coord_sys EQ 'galactic' THEN glactc,pix_ra,pix_dec,2000.,pix_ra,pix_dec,2, /degree
    IF coord_sys EQ 'equatorial' THEN Hor2Eq,pix_dec,pix_ra,Jdate_use,pix_ra,pix_dec,lat=obs.lat,lon=obs.lon,alt=obs.alt,precess=1,/nutate
  
    apply_astrometry, obs, ra_arr=pix_ra, dec_arr=pix_dec, x_arr=xv_hpx, y_arr=yv_hpx, /ad2xy
  
    hpx_i_use=where((xv_hpx GT 0) AND (xv_hpx LT (dimension-1)) AND (yv_hpx GT 0) AND (yv_hpx LT (elements-1)),n_hpx_use)
    IF n_hpx_use EQ 0 THEN BEGIN
      print,"Error: Map has no valid Healpix indices"
    ENDIF
    xv_hpx=xv_hpx[hpx_i_use]
    yv_hpx=yv_hpx[hpx_i_use]
  
  endfor


end