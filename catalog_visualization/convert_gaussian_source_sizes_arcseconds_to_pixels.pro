; This function replicates part of source_dft_multi

pro convert_gaussian_source_sizes_arcseconds_to_pixels

  restore, '/Users/ruby/EoR/extended_source_models_from_Ben_Fall2018/FornaxA_gaussian_model.sav', /relaxed
  restore, '/Users/ruby/EoR/gaussian_model_debugging_Dec18/1130776864_obs.sav', /relaxed
  source_array = []
  for i=0,n_elements(catalog)-1 do begin
    source_array = [source_array, (*catalog[0].extend)[i]]
  endfor

  gaussian_x=make_array(n_elements(source_array),/float)
  gaussian_y=make_array(n_elements(source_array),/float)
  gaussian_rot=make_array(n_elements(source_array),/float)
  for source_ind=0,n_elements(source_array)-1 do begin
    if tag_exist(source_array[source_ind], 'shape') then begin
      ;Convert from FWHM in arcsec to stddev in deg
      gaussian_x[source_ind]=source_array[source_ind].shape.x/(7200*sqrt(2*alog(2)))
      gaussian_y[source_ind]=source_array[source_ind].shape.y/(7200*sqrt(2*alog(2)))
      ;Convert from deg to rad
      gaussian_rot[source_ind]=source_array[source_ind].shape.angle*!Pi/180.
    endif else begin
      gaussian_x[source_ind]=0
      gaussian_y[source_ind]=0
      gaussian_rot[source_ind]=0
    endelse
  endfor
  null=where(gaussian_x,n_gauss_params)
  if n_gauss_params eq 0 then begin
    print, 'Catalog does not contain Gaussian shape parameters. Unsetting keyword gaussian_source_models.'
    undefine, gaussian_source_models
  endif else begin
    ;Convert from deg to pixels
    gaussian_ra_vals = [[source_array.ra+0.5*gaussian_x*Cos(gaussian_rot)], [source_array.ra-0.5*gaussian_x*Cos(gaussian_rot)], $
      [source_array.ra+0.5*gaussian_y*Sin(gaussian_rot)], [source_array.ra-0.5*gaussian_y*Sin(gaussian_rot)]]
    gaussian_dec_vals = [[source_array.dec+0.5*gaussian_x*Sin(gaussian_rot)], [source_array.dec-0.5*gaussian_x*Sin(gaussian_rot)], $
      [source_array.dec+0.5*gaussian_y*Cos(gaussian_rot)], [source_array.dec-0.5*gaussian_y*Cos(gaussian_rot)]]
    apply_astrometry, obs, x_arr=gaussian_x_vals, y_arr=gaussian_y_vals, ra_arr=gaussian_ra_vals, dec_arr=gaussian_dec_vals, /ad2xy
    gaussian_x = sqrt((gaussian_x_vals[*,0]-gaussian_x_vals[*,1])^2.+(gaussian_y_vals[*,0]-gaussian_y_vals[*,1])^2.)
    gaussian_y = sqrt((gaussian_x_vals[*,2]-gaussian_x_vals[*,3])^2.+(gaussian_y_vals[*,2]-gaussian_y_vals[*,3])^2.)
  endelse
  
  for i=0,n_elements(source_array)-1 do begin
    (*catalog[0].extend)[i].shape.x = gaussian_x[i]
    (*catalog[0].extend)[i].shape.y = gaussian_y[i]
  endfor
  
  save, catalog, filename = '/Users/ruby/EoR/gaussian_model_debugging_Dec18/FornaxA_gaussian_model_pixels_units.sav'

end