pro make_gaussian_source_catalog

  restore, catalog_path, /relaxed
  new_struct = {ID:catalog[0].ID, X:catalog[0].X,Y:catalog[0].Y,RA:catalog[0].RA,DEC:catalog[0].DEC,STON:catalog[0].STON,FREQ:catalog[0].FREQ,ALPHA:catalog[0].ALPHA,GAIN:catalog[0].GAIN,FLAG:catalog[0].FLAG,EXTEND:catalog[0].EXTEND,FLUX:catalog[0].FLUX,SHAPE:{x:0.,y:0.,angle:0.}}
  struct_arr = replicate(new_struct,n_elements(catalog))
  struct_assign, catalog, struct_arr
  
  ;put in gaussian source parameters
  struct_arr[1].shape.x=3600.
  struct_arr[1].shape.y=7200.
  
  catalog=struct_arr
  save, catalog, filename=save_filepath
end