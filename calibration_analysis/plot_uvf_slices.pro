pro plot_uvf_slices

  uv_arr = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_traditional_cal_large_window_Jul2018/hex_array_sim_331_even_gridded_uvf.sav', 'weights_uv_arr')
  v_val = 200
  pol = 0
  slice = make_array((size(*uv_arr[pol,0]))[2], (size(uv_arr))[2], /float, value=0.)
  for freq =0, (size(uv_arr))[2]-1 do begin
    slice[*, freq] = reform(abs((*uv_arr[pol,freq])[*,v_val]))
  endfor
  stop
end