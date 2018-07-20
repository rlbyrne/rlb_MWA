;Script that replaces sources in a catalog.
;Limitations: This script only works when the original and replacement sources have a 1-to-1 correspondence
;Written by R. Byrne 09/17

pro catalog_find_and_replace

  delete_only = 1
  add_only = 1
  add_the_weird_screwy_fornax_model = 1
  
  original_catalog_path = '/nfs/eor-00/h1/rbyrne/catalogs/GLEAM_plus_nofornax.sav'
  source_catalog_path = ''
  output_path = '/nfs/eor-00/h1/rbyrne/catalogs/GLEAM_plus_extended_fornax.sav'
  
  target_ra = 50.7
  target_dec = -37.2
  search_radius_deg = .5
  target_flux = 100
  flux_tolerance_jy = 300
  
  original_catalog = getvar_savefile(original_catalog_path, 'catalog')
  
  original_candidates = []
  for i = 0, n_elements(original_catalog)-1 do begin
    distance = sqrt((original_catalog[i].ra-target_ra)^2+(original_catalog[i].dec-target_dec)^2)
    flux_distance = abs(original_catalog[i].flux.I-target_flux)
    if distance lt search_radius_deg and flux_distance lt flux_tolerance_jy then begin
      original_candidates = [original_candidates, i]
    endif
  endfor
  
  if n_elements(original_candidates) gt 1 then begin
    print, 'WARNING: ' + string(n_elements(original_candidates)) + ' original source candidates found.'
    print, 'ORIGINAL SOURCE CANDIDATES:'
    for i = 0, n_elements(original_candidates)-1 do begin
      print, string(i+1)+')'
      print, 'RA: '+string(original_catalog[original_candidates[i]].ra)
      print, 'Dec: '+string(original_catalog[original_candidates[i]].dec)
      print, 'Flux: '+string(original_catalog[original_candidates[i]].flux.I)
      print, '..........'
    endfor
    choose_og_source = ''
    read, choose_og_source, prompt='Choose which source to replace: (<enter> to exit)'
    if choose_og_source eq '' then return
    original_source = original_candidates[fix(choose_og_source)-1]
    print, 'USING SOURCE CANDIDATE NUMBER ' + string(fix(choose_og_source))
    print, 'RA: '+string(original_catalog[original_source].ra)
    print, 'Dec: '+string(original_catalog[original_source].dec)
    print, 'Flux: '+string(original_catalog[original_source].flux.I)
  endif else begin
    if n_elements(original_candidates) eq 0 then begin
      print, 'ERROR: No original source found to replace. Consider relaxing search parameters.'
      if add_only ne 1 then return
    endif else original_source = original_candidates[0]
  endelse
  
  if delete_only eq 0 then begin
    source_catalog = getvar_savefile(source_catalog_path, 'catalog')
    source_candidates = []
    for i = 0, n_elements(source_catalog) do begin
      distance = sqrt((source_catalog[i].ra-original_catalog[original_source].ra)^2+(source_catalog[i].dec-original_catalog[original_source].dec)^2)
      flux_distance = abs(source_catalog[i].flux.I-original_catalog[original_source].flux.I)
      if distance lt search_radius_deg and flux_distance lt flux_tolerance_jy then begin
        source_candidates = [source_candidates, i]
      endif
    endfor
    
    if n_elements(source_candidates) gt 1 then begin
      print, 'WARNING: ' + string(n_elements(source_candidates)) + ' replacement source candidates found.'
      print, 'REPLACEMENT SOURCE CANDIDATES:'
      for i = 0, n_elements(source_candidates)-1 do begin
        print, string(i+1)+')'
        print, 'RA: '+string(source_catalog[source_candidates[i]].ra)
        print, 'Dec: '+string(source_catalog[source_candidates[i]].dec)
        print, 'Flux: '+string(source_catalog[source_candidates[i]].flux.I)
        print, '..........'
      endfor
      choose_replacement_source = ''
      read, choose_replacement_source, prompt='Choose which source to replace with: (<enter> to exit)'
      if choose_replacement_source eq '' then return
      replacement_source = source_candidates[fix(choose_replacement_source)-1]
      print, 'USING SOURCE CANDIDATE NUMBER ' + string(fix(choose_replacement_source))
      print, 'RA: '+string(source_catalog[replacement_source].ra)
      print, 'Dec: '+string(source_catalog[replacement_source].dec)
      print, 'Flux: '+string(source_catalog[replacement_source].flux.I)
    endif else begin
      if n_elements(source_candidates) eq 0 then begin
        print, 'ERROR: No replacement sources found. Consider relaxing search parameters.'
        return
      endif else replacement_source = source_candidates[0]
    endelse
  endif
  
  catalog = original_catalog
  if add_only eq 0 then begin
    remove, original_source, catalog
  endif
  
  if delete_only eq 0 then begin
    catalog = [catalog, source_catalog[replacement_source]]
  endif
  
  if add_the_weird_screwy_fornax_model eq 1 then begin
  
    ra_min = 50.0
    ra_max = 51.5
    dec_min = -37.7
    dec_max = -36.7
    
    source_catalog = getvar_savefile('/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/master_sgal_fornax_cat.sav','catalog',/compatibility_mode)
    component_array = []
    for i = 0, n_elements(source_catalog)-1 do begin
      if source_catalog[i].ra gt ra_min and source_catalog[i].ra lt ra_max and source_catalog[i].dec gt dec_min and source_catalog[i].dec lt dec_max then begin
        component_array = [component_array, source_catalog[i]]
      endif
    endfor
    extend = ptr_new(component_array)
    
    ;sum Stokes I flux of the components to get total flux
    ;centroid the source based on the Stokes I flux of the components
    ra_weighted_sum = 0.
    dec_weighted_sum = 0.
    flux_sum = 0.
    for comp = 0, n_elements(component_array)-1 do begin
      flux_I_comp = (component_array[comp].flux).I
      flux_sum = flux_sum + flux_I_comp
      ra_weighted_sum = ra_weighted_sum + (component_array[comp].ra)*flux_I_comp
      dec_weighted_sum = dec_weighted_sum + (component_array[comp].dec)*flux_I_comp
    endfor
    ra_deg_source = ra_weighted_sum/flux_sum
    dec_deg_source = dec_weighted_sum/flux_sum
    freq_MHz_source = component_array[0].freq ;assume the frequency is equal for all components
    
    flux_source = source_catalog[0].flux
    flux_source.XX = 0.
    flux_source.YY = 0.
    flux_source.XY = 0.
    flux_source.YX = 0.
    flux_source.I = flux_sum
    flux_source.Q = 0.
    flux_source.U = 0.
    flux_source.V = 0.
    
    source = source_catalog[0]
    source.id = 0
    source.x = 0.
    source.y = 0.
    source.ra = ra_deg_source
    source.dec = dec_deg_source
    source.ston = 0.
    source.freq = freq_MHz_source
    source.alpha = -.8
    source.gain = 1
    source.flag = 0
    source.extend = extend
    source.flux = flux_source
    
    catalog = [catalog, source]
    
  endif
  
  print, 'Saving new catalog to ' + output_path
  save, filename = output_path, catalog
  
end