pro catalog_sav_to_txt

  ;catalog_name = 'GLEAM'
  ;catalog_path = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/GLEAMIDR4_181_consistent.sav'
  catalog_name = '1130781304_run1_catalog'
  catalog_path = '/nfs/mwa-08/d1/DiffuseSurvey2015/1130781304_run1_catalog.sav'
  restore, catalog_path, /relaxed
  
  openw, outfile, '/nfs/eor-00/h1/rbyrne/catalog_data/'+catalog_name+'.txt', /get_lun
  printf, outfile, 'RA,DEC,FLUX,EXTENDED,COMPONENTS,RA_1,DEC_1,FLUX_1,RA_2,DEC_2,FLUX_2,...'
  for i = 0, n_elements(catalog)-1 do begin
    data = number_formatter(catalog[i].ra)+','+number_formatter(catalog[i].dec)+','+number_formatter((catalog[i].flux).I)
    if catalog[i].extend ne ptr_new() then begin
      data = data + ',1'
      for j = 0, n_elements(*catalog[i].extend)-1 do begin
        data = data + ',' + number_formatter(n_elements(*catalog[i].extend)) +',' + number_formatter((*catalog[i].extend)[j].ra)+','+number_formatter((*catalog[i].extend)[j].dec)+','+number_formatter(((*catalog[i].extend)[j].flux).I)
      endfor
    endif else begin
      data = data + ',0'
    endelse
    printf, outfile, data
  endfor
  free_lun, outfile
  
end