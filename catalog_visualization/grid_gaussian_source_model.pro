; Script that converts gaussian source models into a gridded point source component model
; Note that this script does not support gaussian angle parameters
; This script works for a single source with gaussian component models
; Written by R. Byrne, Jan. 2019

pro grid_gaussian_source_model
  
  catalog_path = '/Users/ruby/EoR/source_models_for_supplementing_GLEAM/CenA_gaussian_model.sav'
  source_index = 0
  grid_resolution = 0.01 ;Resolution of the gridded model in degrees
  flux_cut = 0.01
  template_catalog_path = '/Users/ruby/EoR/source_models_for_supplementing_GLEAM/FornaxA_from_FHD.sav'
  
  catalog = getvar_savefile(catalog_path, 'catalog')
  source_ra = catalog[source_index].ra
  source_dec = catalog[source_index].dec
  components = *catalog[source_index].extend
  
  comp_ras = []
  comp_decs = []
  comp_fluxes = []
  comp_sizes_x = []
  comp_sizes_y = []
  for comp_ind = 0,n_elements(components)-1 do begin
    ra = components[comp_ind].ra
    if ra gt source_ra+180. then ra=ra-360.
    if ra lt source_ra-180. then ra=ra+360.
    comp_ras = [comp_ras, ra]
    comp_decs = [comp_decs, components[comp_ind].dec]
    comp_fluxes = [comp_fluxes, components[comp_ind].flux.I]
    ;Convert size parameters to stddev in degrees
    comp_sizes_x = [comp_sizes_x, components[comp_ind].shape.x/(7200.*sqrt(2.*alog(2.)))]
    comp_sizes_y = [comp_sizes_y, components[comp_ind].shape.y/(7200.*sqrt(2.*alog(2.)))]
  endfor
  
  ra_vals = findgen(floor((max(comp_ras)-min(comp_ras)+2)/grid_resolution), start=min(comp_ras)-1., increment=grid_resolution)
  dec_vals = findgen(floor((max(comp_decs)-min(comp_decs)+2)/grid_resolution), start=min(comp_decs)-1., increment=grid_resolution)
  gridded_vals = fltarr(n_elements(ra_vals), n_elements(dec_vals))
  
  for comp_ind = 0,n_elements(components)-1 do begin
    if comp_sizes_x[comp_ind] gt 0 and comp_sizes_y[comp_ind] gt 0 then begin
      for ra_ind = 0,n_elements(ra_vals)-1 do begin
        for dec_ind = 0,n_elements(dec_vals)-1 do begin
          pixel_val = comp_fluxes[comp_ind]*grid_resolution^2./(2.*!Pi*comp_sizes_x[comp_ind]*comp_sizes_y[comp_ind]) $
            * exp(-(ra_vals[ra_ind]-comp_ras[comp_ind])^2./(2.*comp_sizes_x[comp_ind]^2.)) $
            * exp(-(dec_vals[dec_ind]-comp_decs[comp_ind])^2./(2.*comp_sizes_y[comp_ind]^2.))
          gridded_vals[ra_ind, dec_ind] = gridded_vals[ra_ind, dec_ind] + pixel_val
        endfor
      endfor
    endif else begin
      ra_ind = value_locate(ra_vals, comp_ras[comp_ind])
      dec_ind = value_locate(dec_vals, comp_decs[comp_ind])
      gridded_vals[ra_ind, dec_ind] = gridded_vals[ra_ind, dec_ind] + comp_fluxes[comp_ind]
    endelse
  endfor
  
  template_catalog = getvar_savefile(template_catalog_path, 'catalog', /compatibility_mode)
  template_catalog = template_catalog[0]
  template_catalog.ra = catalog[source_index].ra
  template_catalog.dec = catalog[source_index].dec
  template_catalog.freq = catalog[source_index].freq
  template_catalog.flux = catalog[source_index].flux
  
  template_component = (*template_catalog.extend)[0]
  output_components = replicate(template_component, n_elements(gridded_vals))
  source_boundary_warning = 0
  n_comps = 0
  for ra_ind = 0,n_elements(ra_vals)-1 do begin
    for dec_ind = 0,n_elements(dec_vals)-1 do begin
      if gridded_vals[ra_ind, dec_ind] gt flux_cut then begin
        if ra_ind eq 0 or dec_ind eq 0 or ra_ind eq n_elements(ra_vals)-1 or dec_ind eq n_elements(dec_vals)-1 then source_boundary_warning = 1
        output_components[n_comps].ra = ra_vals[ra_ind]
        output_components[n_comps].dec = dec_vals[dec_ind]
        output_components[n_comps].freq = catalog[source_index].freq
        output_components[n_comps].flux.I = gridded_vals[ra_ind, dec_ind]
        n_comps += 1
      endif
    endfor
  endfor
  if source_boundary_warning then print, 'WARNING: Source components extend to the source boundary.'
  output_components = output_components[0:n_comps-1]
  stop
  template_catalog.extend = ptr_new(output_components)
  
  catalog = template_catalog
  save, catalog, filename = '/Users/ruby/EoR/source_models_for_supplementing_GLEAM/CenA_gridded_model_'+strtrim(n_comps,1)+'comp.sav'

end