; This function performs absolute calibration based on a chi-squared minimization where only the abs cal parameters are varied
; It assumes the model visibilities are already relatively calibrated
; It does not support flagging or weighting
; Written by R. Byrne, 1/19

pro absolute_calibration

  obsid = 'hex_array_sim_331'
  model_run_path = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_Jun2018'
  reference_run_path = '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_reference_Jun2018'

  restore, model_run_path+'/vis_data/'+obsid+'_vis_model_XX.sav'
  model_vis = *vis_model_ptr
  undefine, vis_model_ptr
  obs_model = getvar_savefile(model_run_path+'/metadata/'+obsid+'_obs.sav', 'obs')
  restore, reference_run_path+'/vis_data/'+obsid+'_vis_model_XX.sav'
  data_vis = *vis_model_ptr
  undefine, vis_model_ptr
  obs_data = getvar_savefile(reference_run_path+'/metadata/'+obsid+'_obs.sav', 'obs')
  
  antenna_a = (*obs_model.baseline_info).tile_a
  antenna_b = (*obs_model.baseline_info).tile_b
  
  antenna_a_data = (*obs_data.baseline_info).tile_a
  antenna_b_data = (*obs_data.baseline_info).tile_b
  
  ;if max(antenna_a ne antenna_a_data) gt 0 or max(antenna_b ne antenna_b_data) gt 0 then print, 'ERROR: antenna arrays do not match!'
  undefine, antenna_a_data
  undefine, antenna_b_data
  
  ; Remove autos
  keep_inds = where(antenna_a ne antenna_b)
  model_vis = model_vis[*, keep_inds]
  data_vis = data_vis[*, keep_inds]
  antenna_a = antenna_a[keep_inds]
  antenna_b = antenna_b[keep_inds]

  ; Get antenna positions
  csv_data = read_csv('/Users/rubybyrne/array_simulation_331/'+obsid+'_antenna_locs.csv', header=csv_header)
  ant_pos = make_array(n_elements(csv_data.field1), 2, /float, value=0.)
  for ant_ind=0, n_elements(csv_data.field1)-1 do begin
    ant_pos[fix(csv_data.field1[ant_ind]), 0] = double(csv_data.field2[ant_ind])
    ant_pos[fix(csv_data.field1[ant_ind]), 1] = double(csv_data.field3[ant_ind])
  endfor
  
  baseline_lengths = make_array(n_elements(antenna_a), 2, /float, value=0.)
  for baseline=0,n_elements(antenna_a)-1 do begin
    baseline_lengths[baseline, *] = ant_pos[antenna_b[baseline]-1, *] - ant_pos[antenna_a[baseline]-1, *]
  endfor
  
  n_frequencies = (size(model_vis))[1]
  delta_x = make_array(n_frequencies, /float, value=0.)
  delta_y = make_array(n_frequencies, /float, value=0.)
  amp2 = make_array(n_frequencies, /float, value=0.)
  for freq=0,n_frequencies-1 do begin
    grad_matrix_element_1 = total(baseline_lengths[*, 0]^2.*real_part(reform(conj(data_vis[freq, *])*model_vis[freq, *])))
    grad_matrix_element_2 = total(baseline_lengths[*, 0]*baseline_lengths[*, 1]*real_part(reform(conj(data_vis[freq, *])*model_vis[freq, *])))
    grad_matrix_element_3 = total(baseline_lengths[*, 1]^2.*real_part(reform(conj(data_vis[freq, *])*model_vis[freq, *])))
    grad_matrix = [[grad_matrix_element_1, grad_matrix_element_2], [grad_matrix_element_2, grad_matrix_element_3]]
    grad_vector_element_1 = total(baseline_lengths[*, 0]*imaginary(conj(data_vis[freq, *])*model_vis[freq, *]))
    grad_vector_element_2 = total(baseline_lengths[*, 1]*imaginary(conj(data_vis[freq, *])*model_vis[freq, *]))
    grad_vector = [grad_vector_element_1, grad_vector_element_2]
    result = matrix_multiply(invert(grad_matrix), grad_vector)
    delta_x[freq] = result[0]
    delta_y[freq] = result[1]
    amp2[freq] = total(real_part(exp(complex(0,1)*reform(result[0]*baseline_lengths[*, 0]+result[1]*baseline_lengths[*, 1]))*reform(conj(data_vis[freq, *])*model_vis[freq, *]))) $
      / total(abs(reform(model_vis[freq, *]))^2.)
  endfor
  stop
  print, delta_x
  print, delta_y
  print, amp2 
  
end