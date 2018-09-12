pro ps_compare_cut_plots

  yrange = [1.e0, 1.e10]
  datafiles = ['/Volumes/Bilbo/rlb_fhd_outputs/array_simulation/fhd_rlb_array_sim_Barry_effect_hex_array_sim_331__abs_errors_minus_perfect/data/1d_binning/hex_array_sim_331_cubeXX__even_odd_joint_bh_bh_averemove_swbh_res_xx_dencorr_1dkpower.idlsave']
  colors = ['blue', 'blue', 'cyan', 'cyan', 'dark green', 'dark green', 'black']
  linestyles = [0,2,0,2,0,2,0]
  colors = colors[0]
  linestyles = linestyles[0]
  
  cgps_open, '/Users/rubybyrne/array_simulation_331/1d_ps_plots/testplot.png'
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
    
    cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=[1e-1, 1e0], $
      ytitle=textoidl('P_k (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k (!8h!X Mpc^{-1})'), $
      linestyle=linestyles[file_ind], color=colors[file_ind], thick=4
  endfor
  cglegend, title=['Hex. array, errors in abs. cal. params.'], linestyle=linestyles, thick=4, $
    color=colors, length=0.03, /center_sym, location=[.55,.8], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps;, density=plot_resolution
    

end