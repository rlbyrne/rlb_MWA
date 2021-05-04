pro fix_calibration_normalization

  base_path = '/Volumes/Bilbo/rlb_fhd_outputs/calibration_tests_Feb2021/'
  reference_run = 'fhd_rlb_perfreq_cal_Feb2021'
  correction_run = 'fhd_rlb_perfreq_all_baseline_cal_with_polarized_diffuse_Mar2021'
  
  ; Get obsids
  filenames = file_search(base_path+correction_run+'/calibration/*_cal.sav')
  obslist = [strmid(filenames, strlen(base_path+correction_run+'/calibration/'), 10)]
  
  for obsind=0,n_elements(obslist)-1 do begin
    obsid = obslist(obsind)
    cal_reference = getvar_savefile(base_path+reference_run+'/calibration/'+obsid+'_cal.sav', 'cal')
    cal_reference_amp = [mean(abs(*cal_reference.gain[0])), mean(abs(*cal_reference.gain[1]))]
    cal_to_correct = getvar_savefile(base_path+correction_run+'/calibration/'+obsid+'_cal.sav', 'cal')
    cal_to_correct_amp = [mean(abs(*cal_to_correct.gain[0])), mean(abs(*cal_to_correct.gain[1]))]
    correction_factor = cal_reference_amp/cal_to_correct_amp
    
    cal = cal_to_correct
    *cal.gain[0] = correction_factor[0]*(*cal_to_correct.gain[0])
    *cal.gain[1] = correction_factor[1]*(*cal_to_correct.gain[1])
    save, cal, filename=base_path+correction_run+'/calibration/calibration_gain_amp_adjusted/'+obsid+'_cal.sav'
  endfor

end