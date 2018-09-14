pro ps_compare_cut_plots

  yrange = [1.e0, 1.e10]
  datafiles = ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_hex_array_sim_331__abs_errors_minus_perfect/data/1d_binning/hex_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_split_hex_array_sim_331__abs_errors_minus_perfect/data/1d_binning/split_hex_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_random1_array_sim_331__abs_errors_minus_perfect/data/1d_binning/random1_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_hex_array_sim_331__traditional_minus_perfect/data/1d_binning/hex_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_split_hex_array_sim_331__traditional_minus_perfect/data/1d_binning/split_hex_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_random1_array_sim_331__traditional_minus_perfect/data/1d_binning/random1_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave']
  colors = ['blue', 'red', 'orange', 'blue', 'red', 'orange']
  linestyles = [0,0,0,2,2,2]
  
  eor_file = '/Users/rubybyrne/FHD/catalog_data/eor_power_1d.idlsave'
  eor_k_centers = getvar_savefile(eor_file, 'k_centers')
  eor_power = getvar_savefile(eor_file, 'power')
  eor_plot_x = []
  eor_plot_y = []
  for datapoint = 0,n_elements(eor_power)-1 do begin
    eor_plot_y = [eor_plot_y, eor_power[datapoint], eor_power[datapoint]]
    if (datapoint ne 0) and (datapoint ne n_elements(eor_power)-1) then begin
      eor_plot_x = [eor_plot_x, (eor_k_centers[datapoint]+eor_k_centers[datapoint-1])/2, (eor_k_centers[datapoint]+eor_k_centers[datapoint+1])/2]
    endif else begin
      if datapoint eq 0 then eor_plot_x = [eor_plot_x, eor_k_centers[datapoint], (eor_k_centers[datapoint]+eor_k_centers[datapoint+1])/2]
      if datapoint eq n_elements(eor_power)-1 then eor_plot_x = [eor_plot_x, (eor_k_centers[datapoint]+eor_k_centers[datapoint-1])/2, eor_k_centers[datapoint]]
    endelse
  endfor
      
  cgps_open, '/Users/rubybyrne/array_simulation_331/1d_ps_plots/testplot.png'
  cgplot, eor_plot_x, eor_plot_y, /xlog, /ylog, yrange=yrange, xrange=[7e-2, 7e-1], $
    linestyle=0, color='black', thick=4, $
    ytitle=textoidl('P_k (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k (!8h!X Mpc^{-1})')
  for file_ind = 0,n_elements(datafiles)-1 do begin
    k_edges = getvar_savefile(datafiles[file_ind], 'k_edges')
    power = getvar_savefile(datafiles[file_ind], 'power')
    power[where(power lt yrange[0]/100.)] = yrange[0]/100.  ; Remove negative values for log plot

    plot_x = []
    plot_y = []
    for datapoint = 0,n_elements(power)-1 do begin
      plot_y = [plot_y, power[datapoint], power[datapoint]]
      plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
    endfor
    cgplot, plot_x, plot_y, /xlog, /ylog, $
      linestyle=linestyles[file_ind], color=colors[file_ind], thick=4, /overplot
  endfor
  cglegend, title=['Hex. array, per ant. cal.', 'Offset hex. array, per ant. cal.', 'Random array, per ant. cal.',$
    'Hex. array, errors in abs. cal. params.', 'Offset hex. array, errors in abs. cal. params.', 'Random array, errors in abs. cal. params.', 'Predicted EoR'], $
    linestyle=[2,2,2,0,0,0,0], thick=4, $
    color=[colors, 'black'], length=0.03, /center_sym, location=[.55,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=600
    

end