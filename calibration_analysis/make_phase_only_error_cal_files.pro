pro make_phase_only_error_cal_files

  cal_file_save_path = '/Users/rubybyrne/array_simulation_331/cal_files'
  obsids = ['hex_array_sim_331', 'split_hex_array_sim_331', 'random1_array_sim_331']
  for obs =0,n_elements(obsids)-1 do begin
    obsid = obsids[obs]
    cal = getvar_savefile('/Users/rubybyrne/array_simulation_331/cal_files/'+obsid+'_cal_abs_errors_only.sav', 'cal')
    for pol = 0,1 do begin
      for i = 0, (size(*cal.gain[pol]))[1]-1 do begin
        for j = 0, (size(*cal.gain[pol]))[2]-1 do begin
          (*cal.gain[pol])[i,j] = (*cal.gain[pol])[i,j]/abs((*cal.gain[pol])[i,j])
        endfor
      endfor
    endfor
    
    print, 'Saving file to '+cal_file_save_path+'/'+obsid+'_cal_phase_errors_only.sav'
    save, cal, filename=cal_file_save_path+'/'+obsid+'_cal_phase_errors_only.sav'

  endfor

end