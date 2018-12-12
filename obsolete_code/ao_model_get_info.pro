pro ao_model_get_info

  model_name = 'PicA-88comp'
  filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/anoko/mwa-reduce/models/model-'+model_name+'.txt'
  
  original_catalog_filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/GLEAMIDR4_181_consistent.sav'
  desination_catalog_filepath = '/nfs/eor-00/h1/rbyrne/catalogs/GLEAM+'+model_name+'.sav'
  restore, original_catalog_filepath, /relaxed
  
  ;;;;parameters I don't know so I'm setting arbitrarily;;;;
  id = 0
  x = 0.
  y = 0.
  ston = 0.
  alpha = -.8
  gain  = 1.
  flag = 0
  
  nlines = file_lines(filepath)
  file_contents = strarr(nlines)
  openr, datafile, filepath, /get_lun
  readf, datafile, file_contents
  free_lun, datafile
  
  source_start_lines = []
  component_start_lines = []
  for i = 0, nlines-1 do begin
    if strpos(file_contents[i], 'source {') ne -1 then begin
      source_start_lines = [source_start_lines, i]
    endif
    if strpos(file_contents[i], 'component {') ne -1 then begin
      component_start_lines = [component_start_lines, i]
    endif
  endfor
  source_start_lines = [source_start_lines, nlines-1]
  component_start_lines = [component_start_lines, nlines-1]
  
  catalog_addition = []
  for i = 0, n_elements(source_start_lines)-2 do begin
    component_array = []
    for j = 0, n_elements(component_start_lines)-2 do begin
      if (component_start_lines[j] gt source_start_lines[i]) and $
        (component_start_lines[j] lt source_start_lines[i+1]) then begin
        component_data = file_contents[component_start_lines[j]:min([component_start_lines[j+1],source_start_lines[i+1]])-1]
        
        ra_deg = 0
        dec_deg = 0
        freq_MHz = 0
        flux_I_Jy = 0
        flux_Q_Jy = 0
        flux_U_Jy = 0
        flux_V_Jy = 0
        
        for line = 0, n_elements(component_data)-1 do begin
        
          if strpos(component_data[line], 'position') ne -1 then begin
            pos_data = strsplit(component_data[line], /extract)
            ra = pos_data[1]
            ra_split = strsplit(ra,'h',/extract)
            ra_split = [ra_split[0], strsplit(ra_split[1],'m',/extract)]
            ra_split = [ra_split[0:1],strsplit(ra_split[2],'s',/extract)]
            ra_deg = (float(ra_split[0]) + float(ra_split[1])/60. + float(ra_split[2])/3600.)*360./24.
            dec = pos_data[2]
            dec_split = strsplit(dec,'d',/extract)
            dec_split = [dec_split[0], strsplit(dec_split[1],'m',/extract)]
            dec_split = [dec_split[0:1], strsplit(dec_split[2],'s',/extract)]
            dec_deg = abs(float(dec_split[0])) + float(dec_split[1])/60. + float(dec_split[2])/3600.
            if float(dec_split[0]) lt 0 then dec_deg = -dec_deg
          endif
          
          if strpos(component_data[line], 'frequency') ne -1 then begin
            freq_MHz = float((strsplit(component_data[line], /extract))[1])
          endif
          
          if strpos(component_data[line], 'fluxdensity') ne -1 then begin
            flux_data = strsplit(component_data[line], /extract)
            flux_I_Jy = float(flux_data[2])
            flux_Q_Jy = float(flux_data[3])
            flux_U_Jy = float(flux_data[4])
            flux_V_Jy = float(flux_data[5])
            
            flux_comp = catalog[0].flux
            flux_comp.XX = 0.
            flux_comp.YY = 0.
            flux_comp.XY = 0.
            flux_comp.YX = 0.
            flux_comp.I = flux_I_Jy
            flux_comp.Q = flux_Q_Jy
            flux_comp.U = flux_U_Jy
            flux_comp.V = flux_V_Jy
          ;flux_comp = {XX:0., YY:0., XY: 0., YX: 0., I:flux_I_Jy, Q:flux_Q_Jy, U:flux_U_Jy, V:flux_V_Jy}
          endif
          
        endfor
        
        component = catalog[0]
        component.id = id
        component.x = x
        component.y = y
        component.ra = ra_deg
        component.dec = dec_deg
        component.ston = ston
        component.freq = freq_MHz
        component.alpha = alpha
        component.gain = gain
        component.flag = flag
        component.extend = ptr_new()
        component.flux = flux_comp
        
        component_array = [component_array, component]
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
    
    print, "....."
    print, model_name + " RA: " + ra_deg_source + " Dec: " + dec_deg_source + " Flux: " + flux_sum
    print, "....."
  endfor

end