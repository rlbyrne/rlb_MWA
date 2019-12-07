pro kpar0_plots_for_salf

  yrange = [1.e13, 1.e16]
  datafiles = ['/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_dirty_xx_averemove_swbh_dencorr_k0power.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_model_xx_averemove_swbh_dencorr_k0power.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_dirty_yy_averemove_swbh_dencorr_k0power.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_model_yy_averemove_swbh_dencorr_k0power.idlsave']
  ;colors = ['blue', 'red', 'orange', 'blue', 'red', 'orange']
  colors = ['blue', 'blue', 'red', 'red']
  ;colors = ['purple']
  ;linestyles = [0,0,0,2]
  ;linestyles = [0,0,0,2,2,2]
  linestyles = [0,2,0,2]

  cgps_open, '/Users/rubybyrne/kpar0_plot_for_salf.png'
  ;cgplot, eor_plot_x, eor_plot_y, /xlog, /ylog, yrange=yrange, xrange=[7e-2, 7e-1], $
  ;  linestyle=0, color='black', thick=4, $
  ;  ytitle=textoidl('P_k (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), /nodata
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
    if file_ind eq 0 then overplot=0 else overplot=1
    cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=[8e-4, 2e-1], $
      linestyle=linestyles[file_ind], color=colors[file_ind], thick=4, overplot=overplot, title='k-parallel=0 Power Spectrum', $ 
      ytitle=textoidl('P_k (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})')
  endfor
  cglegend, title=['Data E-W pol', 'Model E-W pol', 'Data N-S pol', 'Model N-S pol'], $
    linestyle=[0,2,0,2], thick=4, $
    color=[colors], length=0.03, /center_sym, location=[.7,.9], charsize=.8, /box, background='white'
  cgps_close, /png, /delete_ps, density=600


end