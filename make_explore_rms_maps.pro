pro make_explore_rms_maps

  map_orig_path = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/average_map_1131454296_rm_undone_IQU.sav'
  outpath = '/Users/rubybyrne/diffuse_survey_plotting_Aug2020/explore_rms_maps'
  map_orig = getvar_savefile(map_orig_path, 'model_arr')
  hpx_inds = getvar_savefile(map_orig_path, 'hpx_inds')
  nside = getvar_savefile(map_orig_path, 'nside')
  
  angles = ['45', '90', '135', '180', '225', '270', '315']
  for ang_ind=0,n_elements(angles)-1 do begin
    model_arr = ptrarr(4)
    model_arr[0] = ptr_new(*map_orig[0])
    model_arr[3] = ptr_new(*map_orig[3])
    model_arr[1] = ptr_new(cos(float(angles[ang_ind])*!DtoR)*(*map_orig[1]) - sin(float(angles[ang_ind])*!DtoR)*(*map_orig[2]))
    model_arr[2] = ptr_new(sin(float(angles[ang_ind])*!DtoR)*(*map_orig[1]) + cos(float(angles[ang_ind])*!DtoR)*(*map_orig[2]))
    save_path = outpath+'/average_map_rot_angle_'+angles[ang_ind]+'.sav
    print, 'Saving map to '+save_path
    save, model_arr, hpx_inds, nside, filename=save_path
  endfor
  
end