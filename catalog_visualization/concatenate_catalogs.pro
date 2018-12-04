; Script to concatenate catalogs with conflicting data structures
; Adds a catalog that includes shape structure (catalog1) to one that does not (catalog2)
; Written by Ruby Byrne, 12/18

pro concatenate_catalogs

  catalog1_path = '/Users/ruby/EoR/extended_source_models_from_Ben_Fall2018/FornaxA_from_FHD.sav'
  catalog2_path = '/Users/ruby/EoR/GLEAM_EGC_v2_181MHz.sav'
  
  restore, catalog1_path, /relaxed
  catalog1 = catalog
  restore, catalog2_path, /relaxed
  catalog2 = catalog
  
  catalog2_reformatted=Replicate(catalog1[0],n_elements(catalog2)>1)
  for source=0,n_elements(catalog2)-1 do begin
    catalog2_reformatted[source].id = catalog2[source].id
    catalog2_reformatted[source].x = catalog2[source].x
    catalog2_reformatted[source].y = catalog2[source].y
    catalog2_reformatted[source].ra = catalog2[source].ra
    catalog2_reformatted[source].dec = catalog2[source].dec
    catalog2_reformatted[source].ston = catalog2[source].ston
    catalog2_reformatted[source].freq = catalog2[source].freq
    catalog2_reformatted[source].alpha = catalog2[source].alpha
    catalog2_reformatted[source].gain = catalog2[source].gain
    catalog2_reformatted[source].flag = catalog2[source].flag
    catalog2_reformatted[source].extend = catalog2[source].extend
    catalog2_reformatted[source].flux = catalog2[source].flux
  endfor
  
  catalog = [catalog1, catalog2_reformatted]
  save, catalog, filename='/Users/ruby/EoR/extended_source_models_from_Ben_Fall2018/GLEAM_v2_plus_FornaxA_from_FHD.sav'

end