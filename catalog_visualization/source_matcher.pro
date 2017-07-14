pro source_matcher

  source_ra = 0.
  source_dec = 0.
  source_flux = 0.
  matching_catalog_filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/GLEAMIDR4_181_consistent.sav'
  restore, matching_catalog_filepath, /relaxed
  
  for index = 0, n_elements(catalog)-1 do begin
    matching_flux = (catalog[index].flux).I
    if matching_flux gt source_flux*0.2 then begin
      matching_ra = catalog[index].ra
      matching_dec = catalog[index].dec
      ra_distance = min([abs(source_ra-matching_ra), abs(source_ra-matching_ra-360.), abs(source_ra-matching_ra+360.)])
      dec_distance = 
    endif
    
  endfor
  
end