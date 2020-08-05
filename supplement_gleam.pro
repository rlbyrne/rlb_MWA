; Script to add additional sources to the GLEAM catalog
; Written by Ruby Byrne, 8/20

pro supplement_gleam
  
  outfile_name = '/Users/ruby/EoR/gleam_catalogs/GLEAM_v2_plus_rlb2019_with_labels.sav'
  path_to_gleam = '/Users/ruby/EoR/gleam_catalogs/gleam.sav'
  gleam_cat = getvar_savefile(path_to_gleam, 'catalog', /compatibility_mode)
  
  supplement_source_path = '/Users/ruby/EoR/source_models_for_supplementing_GLEAM'
  supplement_source_names = ['3C161', '3C409', 'CassiopeiaA',$
    'CentaurusA', 'FornaxA', 'HeraA', 'HydraA',$
    'PictorA', 'VirgoA']
  supplement_source_filenames = ['3C161_from_anoko_1comp.sav', '3C409_from_anoko_1comp.sav', 'CasA_from_anoko_1comp.sav',$
    'CenA_gridded_model_5874comp.sav', 'FornaxA_from_FHD.sav', 'HerA_from_anoko_27comp.sav', 'HydA_from_anoko_58comp.sav',$
    'PicA_from_anoko_88comp.sav', 'VirA_from_anoko_171comp.sav']
    
  flux_struct={flux,xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
  struct_base={id:'',x:0.,y:0.,ra:0.,dec:0.,ston:0.,freq:180.,alpha:0.,gain:1.,flag:0,extend:Ptr_new(),flux:flux_struct}
  cat_supplement=Replicate(struct_base,n_elements(supplement_source_filenames)>1)
  
  for source_ind=0,n_elements(supplement_source_filenames)-1 do begin
    source = getvar_savefile(supplement_source_path+'/'+supplement_source_filenames[source_ind], 'catalog', /compatibility_mode)
    cat_supplement[source_ind].id = supplement_source_names[source_ind]
    cat_supplement[source_ind].ra = source.ra
    cat_supplement[source_ind].dec = source.dec
    cat_supplement[source_ind].freq = source.freq
    cat_supplement[source_ind].flux.I = source.flux.I
    if ptr_valid(source.extend) then begin
      comps = Replicate(struct_base,n_elements(*source.extend)>1)
      for comp_ind=0,n_elements(*source.extend)-1 do begin
        comps[comp_ind].id = supplement_source_names[source_ind]+'-'+strtrim(string(comp_ind+1), 2)
        comps[comp_ind].ra = (*source.extend)[comp_ind].ra
        comps[comp_ind].dec = (*source.extend)[comp_ind].dec
        comps[comp_ind].freq = (*source.extend)[comp_ind].freq
        comps[comp_ind].flux.I = (*source.extend)[comp_ind].flux.I
      endfor
      cat_supplement[source_ind].extend = ptr_new(comps)
    endif
  endfor
  
  catalog = [gleam_cat, cat_supplement]
  save, catalog, filename=outfile_name
  
end