; OBSOLETE CODE
; Replaced with catalog_reformat_ao_to_fhd.pro 12/18

;Script that adds sources from the AO-compatible GLEAM extended source models to an FHD catalog.
;Limitations:
;  Only works when the model file only contains one measurement
;  ID, X, Y, STON, ALPHA, GAIN, FLAG are not calculated
;  Calculates flux in Stokes I only
;  Flux in XX, YY, XY, and YX is not calculated for the source or components.
;Written by R. Byrne 06/17

pro ao_model_adder

  model_name = 'PicA-88comp'
  filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/anoko/mwa-reduce/models/model-'+model_name+'.txt'

  target_freq_MHz = 180. ;If more than one measurement exists, use the one taken at the freq closest to this value
  ;original_catalog_filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/GLEAMIDR4_181_consistent.sav'
  original_catalog_filepath = '/nfs/eor-00/h1/rbyrne/MWA/IDL_code/FHD/catalog_data/GLEAM_plus_rlb2017.sav'
  desination_catalog_filepath = '/nfs/eor-00/h1/rbyrne/catalogs/GLEAM_plus.sav'
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
        freq_MHz = []
        flux_I_Jy = []

        for line = 0, n_elements(component_data)-1 do begin

          if strpos(component_data[line], 'position') ne -1 then begin
            pos_data = strsplit(component_data[line], /extract)
            ra = pos_data[1]
            ra_split = strsplit(ra,'h',/extract)
            ra_split = [ra_split[0], strsplit(ra_split[1],'m',/extract)]
            ra_split = [ra_split[0:1],strsplit(ra_split[2],'s',/extract)]
            ra_deg = (abs(float(ra_split[0])) + float(ra_split[1])/60. + float(ra_split[2])/3600.)*360./24.
            if strpos(ra_split[0],'-') ne -1 then begin
              ra_deg = -ra_deg
              ra_deg = ra_deg + 360. ;make all the RA values positive
            endif
            dec = pos_data[2]
            dec_split = strsplit(dec,'d',/extract)
            dec_split = [dec_split[0], strsplit(dec_split[1],'m',/extract)]
            dec_split = [dec_split[0:1], strsplit(dec_split[2],'s',/extract)]
            dec_deg = abs(float(dec_split[0])) + float(dec_split[1])/60. + float(dec_split[2])/3600.
            if strpos(dec_split[0],'-') ne -1 then dec_deg = -dec_deg
          endif

          if strpos(component_data[line], 'frequency') ne -1 then begin
            freq_MHz = [freq_MHz, float((strsplit(component_data[line], /extract))[1])]
          endif

          if strpos(component_data[line], 'fluxdensity Jy') ne -1 then begin
            flux_data = strsplit(component_data[line], /extract)
            flux_I_Jy = [flux_I_Jy, float(flux_data[2])]
          endif

        endfor

        if n_elements(freq_MHz) gt 1 then begin
          freq_dist = []
          for freq_index = 0, n_elements(freq_MHz)-1 do begin
            freq_dist = [freq_dist, abs(freq_MHz[freq_index]-target_freq_MHz)]
          endfor
          use_meas = value_locate(freq_dist, min(freq_dist))
        endif else use_meas = 0
        freq_MHz_use = freq_MHz[use_meas]
        flux_I_Jy_use = flux_I_Jy[use_meas]

        flux_comp = catalog[0].flux
        flux_comp.XX = 0.
        flux_comp.YY = 0.
        flux_comp.XY = 0.
        flux_comp.YX = 0.
        flux_comp.I = flux_I_Jy_use
        flux_comp.Q = 0.
        flux_comp.U = 0.
        flux_comp.V = 0.

        component = catalog[0]
        component.id = id
        component.x = x
        component.y = y
        component.ra = ra_deg
        component.dec = dec_deg
        component.ston = ston
        component.freq = freq_MHz_use
        component.alpha = alpha
        component.gain = gain
        component.flag = flag
        component.extend = ptr_new()
        component.flux = flux_comp

        component_array = [component_array, component]
      endif

    endfor

    if n_elements(component_array) gt 1 then extend = ptr_new(component_array) else extend = ptr_new()

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

    flux_source = catalog[0].flux
    flux_source.XX = 0.
    flux_source.YY = 0.
    flux_source.XY = 0.
    flux_source.YX = 0.
    flux_source.I = flux_sum
    flux_source.Q = 0.
    flux_source.U = 0.
    flux_source.V = 0.
    ;flux_source = {XX:0., YY:0., XY: 0., YX: 0., I:flux_sum, Q:0., U:0., V:0.}

    source = catalog[0]
    source.id = id
    source.x = x
    source.y = y
    source.ra = ra_deg_source
    source.dec = dec_deg_source
    source.ston = ston
    source.freq = freq_MHz_source
    source.alpha = alpha
    source.gain = gain
    source.flag = flag
    source.extend = extend
    source.flux = flux_source

    catalog_addition = [catalog_addition, source]

  endfor


  catalog = [catalog, catalog_addition]
  print, 'Saving new catalog to ' + desination_catalog_filepath
  save, filename = desination_catalog_filepath, catalog

end
