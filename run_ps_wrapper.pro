pro run_ps_wrapper

  ;obsids = ['hex_array_sim_20m', 'random1_array_sim_20m', 'random2_array_sim_20m', 'random3_array_sim_20m']
  obsids = ['hex_array_sim_331', 'split_hex_array_sim_331', 'random1_array_sim_331']
  
  for obs = 0, n_elements(obsids)-1 do begin

    ;ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_amp_errors_only_Jul2018', obsids[obs], /png;, image_window_name = 'Blackman-Harris', image_window_frac_size=0.5;, /refresh_ps
    ;ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_Jul2018', obsids[obs], /png;, image_window_name = 'Blackman-Harris', image_window_frac_size=0.5;, /refresh_ps
    ;ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_phase_errors_only_Jul2018', obsids[obs], /png;, image_window_name = 'Blackman-Harris', image_window_frac_size=0.5;, /refresh_ps
    ;ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_Jul2018', obsids[obs], /png;, image_window_name = 'Blackman-Harris', image_window_frac_size=0.5;, /refresh_ps
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_large_window_Jul2018', 'random1_array_sim_331',/png, uvf_input=0, /refresh_info
    
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /refresh_info, /refresh_ps
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_traditional_cal_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /refresh_info, /refresh_ps
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /refresh_info, /refresh_ps
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /plot_2d_masked
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_traditional_cal_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /plot_2d_masked
    ps_wrapper, '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_large_window_Jul2018', obsids[obs], /uvf_input, /png, image_window_name='Blackman-Harris', /plot_2d_masked


    ;ps_diff_wrapper, ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_amp_errors_only_Jul2018','/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_Jul2018'], $
    ;  diff_save_path='/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_'+obsids[obs]+'__amp_errors_only_minus_abs_errors_only/', $
    ;  obsids[obs], /png, /refresh
    ;ps_diff_wrapper, ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_Jul2018','/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_Jul2018'], $
    ;  diff_save_path='/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_'+obsids[obs]+'__perfect_cal_minus_abs_errors_only/', $
    ;  obsids[obs], /png, /refresh
    ;ps_diff_wrapper, ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_Jul2018','/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_amp_errors_only_Jul2018'], $
    ;  diff_save_path='/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_'+obsids[obs]+'__perfect_cal_minus_amp_errors_only/', $
    ;  obsids[obs], /png, /refresh
    
    ;ps_diff_wrapper, ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_abs_errors_only_large_window_Jul2018','/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_perfect_cal_large_window_Jul2018'], diff_save_path='/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_hex_array_sim_331__abs_errors_minus_perfect/', 'hex_array_sim_331', /png, /refresh, image_window_name='Blackman-Harris', /invert_colorbar
  
  endfor
end