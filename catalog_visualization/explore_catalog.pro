pro explore_catalog

  source_ra = 50.5
  source_dec = 37.
  restore, '/Users/rubybyrne/FHD/catalog_data/GLEAM_plus_rlb2017.sav', /relaxed
  new_catalog = []
  fluxes = []
  for index = 0, n_elements(catalog)-1 do begin
    if (abs(catalog[index].dec-source_dec) lt 20. and abs(catalog[index].ra - source_ra) lt 20.) then begin
      stop
      new_catalog = [new_catalog, catalog[index]]
      fluxes = [fluxes, catalog[index].flux.I]
    endif
  endfor
  brightest_indices = (reverse(sort(fluxes)))
  new_catalog_sorted = new_catalog[brightest_indices]
  stop
  
end
