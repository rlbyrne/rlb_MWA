pro run_quickview_alone

  fhd_path = '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Oct2019/fhd_rlb_gaussian_diffuse_small_test_short_baselines_Jan2019'
  obs = getvar_savefile(fhd_path+'/metadata/1130773144_obs.sav', 'obs')
  status_str = getvar_savefile(fhd_path+'/metadata/1130773144_status.sav', 'status_str')
  psf = getvar_savefile(fhd_path+'/beams/1130773144_beams.sav', 'psf')
  cal = 0
  jones = getvar_savefile(fhd_path+'/beams/1130773144_jones.sav', 'jones')
  skymodel = getvar_savefile(fhd_path+'/output_data/1130773144_skymodel.sav', 'skymodel')
  fhd_params = getvar_savefile(fhd_path+'/metadata/1130773144_params.sav', 'params')
  image_uv_arr = ptrarr(2)
  im_uv_xx = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_XX.sav', 'grid_uv')
  im_uv_yy = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_YY.sav', 'grid_uv')
  image_uv_arr[0] = ptr_new(im_uv_xx)
  image_uv_arr[1] = ptr_new(im_uv_yy)
  weights_arr = ptrarr(2)
  weights_xx = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_weights_XX.sav', 'weights_uv')
  weights_yy = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_weights_YY.sav', 'weights_uv')
  weights_arr[0] = ptr_new(weights_xx)
  weights_arr[1] = ptr_new(weights_yy)
  model_uv_arr = ptrarr(2)
  mod_uv_xx = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_model_XX.sav', 'grid_uv_model')
  mod_uv_yy = getvar_savefile(fhd_path+'/grid_data/1130773144_uv_model_YY.sav', 'grid_uv_model')
  model_uv_arr[0] = ptr_new(mod_uv_xx)
  model_uv_arr[1] = ptr_new(mod_uv_yy)
  image_filter_fn = 'filter_uv_weighted'
  file_path_fhd = fhd_path+'/1130773144'
  
  fhd_quickview,obs,status_str,psf,cal,jones,skymodel,fhd_params,image_uv_arr=image_uv_arr,weights_arr=weights_arr,$
    model_uv_arr=model_uv_arr,image_filter_fn=image_filter_fn,file_path_fhd=file_path_fhd
    
  ;weighted_uv_xx = filter_uv_weighted(*model_uv_arr[0],obs=obs,psf=psf,params=params,weights=*weights_arr[0])
  ;weighted_uv_yy = filter_uv_weighted(*model_uv_arr[1],obs=obs,psf=psf,params=params,weights=*weights_arr[1])
  stop

end