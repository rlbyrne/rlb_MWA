pro rlb_versions_gillyweed
  except=!except
  !except=0
  heap_gc
  
  obs_id = '1061316296'
  vis_file_list = '/Users/Shared/uvfits/4.1/'+obs_id+'.uvfits'
  output_directory = '/Volumes/Bilbo/rlb_fhd_outputs'
  version = 'rlb_get_uvf_cubes_Apr2020'

  case version of

    'rlb_get_uvf_cubes_Apr2020': begin
      recalculate_all = 1
      uvfits_version = 4
      uvfits_subversion = 1
      max_sources = 200000
      calibration_catalog_file_path = filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
      diffuse_calibrate = 0
      smooth_width = 32
      pad_uv_image = 1
      return_sidelobe_catalog = 1
      dft_threshold = 0
      ring_radius = 0
      debug_region_grow = 0
      n_pol = 2
      time_cut = -4 ;flag an extra 4 seconds from the end of each obs
      save_uvf = 1
    end

  endcase

  if ~keyword_set(vis_file_list) then begin
    if platform eq 'aws' then begin
      vis_file_list = '/uvfits/' + STRING(obs_id) + '.uvfits'
    endif else begin
      SPAWN, 'read_uvfits_loc.py -v ' + STRING(uvfits_version) + ' -s ' + $
        STRING(uvfits_subversion) + ' -o ' + STRING(obs_id), vis_file_list
    endelse
  endif

  undefine, uvfits_subversion, uvfits_version

  fhd_file_list=fhd_path_setup(vis_file_list,version=version,output_directory=output_directory)
  healpix_path=fhd_path_setup(output_dir=output_directory,subdir='Healpix',output_filename='Combined_obs',version=version)


  ; Set global defaults and bundle all the variables into a structure.
  ; Any keywords set on the command line or in the top-level wrapper will supercede these defaults
  eor_wrapper_defaults,extra
  fhd_depreciation_test, _Extra=extra

  print,""
  print,"Keywords set in wrapper:"
  print,structure_to_text(extra)
  print,""

  general_obs,_Extra=extra

end
