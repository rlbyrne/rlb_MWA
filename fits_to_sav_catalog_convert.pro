; Script to convert fits images into source components and add to FHD-readable source catalogs
; Written by Ruby Byrne, 11/19

pro fits_to_sav_catalog_convert, filename_fits, output_filename

  filename_fits = '/Users/ruby/Downloads/75219_test_file.fits'
  filename_gleam = '/Users/ruby/EoR/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav'
  filename_sav = '/Users/ruby/EoR/kelcey_catalogs/GLEAM_plus_75219.sav'
  
  fits_read, filename_fits, data, header
  ;grab header information:
  crpix1 = sxpar(header, 'crpix1')
  crpix2 = sxpar(header, 'crpix2')
  crval1 = sxpar(header, 'crval1')
  crval2 = sxpar(header, 'crval2')
  cdelt1 = sxpar(header, 'cdelt1')
  cdelt2 = sxpar(header, 'cdelt2')
  naxis1 = sxpar(header, 'naxis1')
  naxis2 = sxpar(header, 'naxis2')
  gleamnum = sxpar(header, 'gleamnum')
  
  ;confirm expected formatting
  if strtrim(sxpar(header, 'ctype1')) ne 'RA' or strtrim(sxpar(header, 'ctype2')) ne 'DEC' $
    or strtrim(sxpar(header, 'cunit1')) ne 'DEG' or strtrim(sxpar(header, 'cunit2')) ne 'DEG' then $
    print, 'ERROR: Unexpected formatting'
    
  xarr = findgen(naxis1)*cdelt1
  xarr += crval1-xarr[crpix1-1]
  yarr = findgen(naxis2)*cdelt2
  yarr += crval2-yarr[crpix2-1]
  
  null = where(data ne 0., ncomps)

  flux_struct = {flux,xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
  ;flux order is 0-3: xx, yy, xy, yx in apparent brightness; 4-7: I, Q, U, V in sky brightness
  ;flag type codes are 0: no flag, 1: low confidence 2: sidelobe contamination
  struct_base = {id:'',x:0.,y:0.,ra:0.,dec:0.,ston:0.,freq:180.,alpha:0.,gain:1.,flag:0,extend:Ptr_new(),flux:flux_struct}

  comps = Replicate(struct_base, ncomps>1)
  compind = 0
  for xind=0,naxis1-1 do begin
    for yind=0,naxis2-1 do begin
      if data[xind, yind] ne 0. then begin
        comps[compind].ra = xarr[xind]
        comps[compind].dec = yarr[yind]
        comps[compind].flux.I = data[xind, yind]
        compind += 1
      endif
    endfor
  endfor

  source_flux = total(data)
  data_xscaled = data
  for xind=0,naxis1-1 do data_xscaled[xind, *] *= xarr[xind]
  source_ra = total(data_xscaled)/source_flux
  data_yscaled = data
  for yind=0,naxis2-1 do data_yscaled[*, yind] *= yarr[yind]
  source_dec = total(data_yscaled)/source_flux
  
  print, source_ra
  print, source_dec
  print, source_flux
  
  source = Replicate(struct_base, 1)
  source.flux.I = source_flux
  source.ra = source_ra
  source.dec = source_dec
  source.extend = ptr_new(comps)
  
  restore, filename_gleam
  print, catalog[gleamnum].ra 
  print, catalog[gleamnum].dec 
  print, catalog[gleamnum].flux.I 
  catalog[gleamnum] = source

  ;save, catalog, filename=filename_sav

end
