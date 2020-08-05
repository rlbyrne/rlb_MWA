; Script to convert .txt files output by pyradiosky to FHD-readable .sav files
; Written by Ruby Byrne, 11/19

pro pyradiosky_to_sav_catalog_convert, filename_txt, filename_sav

    flux_struct={flux,xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
    ;flux order is 0-3: xx, yy, xy, yx in apparent brightness; 4-7: I, Q, U, V in sky brightness
    ;flag type codes are 0: no flag, 1: low confidence 2: sidelobe contamination
    struct_base={id:'',x:0.,y:0.,ra:0.,dec:0.,ston:0.,freq:180.,alpha:0.,gain:1.,flag:0,extend:Ptr_new(),flux:flux_struct}

    nsources = file_lines(filename_txt)-1
    file_contents = strarr(nsources+1)
    openr, datafile, filename_txt, /get_lun
    readf, datafile, file_contents
    free_lun, datafile

    ; Assume header has the form:
    ; SOURCE_ID  RA_J2000 [deg]  Dec_J2000 [deg]  Flux [Jy]  Frequency [Hz]
    data = strarr(5, nsources)
    for line=0,nsources-1 do data[*, line]=strsplit(file_contents[line+1], /extract)
    undefine, file_contents
    source_ids = data[0, *]
    source_ras = float(data[1, *])*!RaDeg
    source_decs = float(data[2, *])
    source_fluxes_I = float(data[3, *])
    source_freqs = float(data[4, *])
    undefine, data

    catalog=Replicate(struct_base,nsources>1)
    for ind=0,nsources-1 do begin
        catalog[ind].id = source_ids[ind]
        catalog[ind].ra = source_ras[ind]
        catalog[ind].dec = source_decs[ind]
        catalog[ind].freq = source_freqs[ind]
        catalog[ind].flux.I = source_fluxes_I[ind]
    endfor
    
    save, catalog, filename=filename_sav

end
