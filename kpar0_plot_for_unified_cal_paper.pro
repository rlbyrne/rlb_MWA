pro kpar0_plot_for_unified_cal_paper

  yrange = [2.e13, 3.e15]
  datafiles = ['/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_dirty_xx_averemove_swbh_dencorr_k0power.idlsave', $
    '/Volumes/Bilbo/rlb_fhd_outputs/fhd_bug_testing_Nov2019/fhd_rlb_master_reference_Nov2019/ps/data/1d_binning/1061316296_cubeXX__even_odd_joint_noimgclip_model_xx_averemove_swbh_dencorr_k0power.idlsave']
  colors = ['blue', 'blue']
  linestyles = [0,2]

  cgps_open, '/Users/rubybyrne/kpar0_plot_for_unified_cal_paper.png'
  for file_ind = 0,n_elements(datafiles)-1 do begin
    k_edges = getvar_savefile(datafiles[file_ind], 'k_edges')
    xrange=[8e-4, 2e-1]
    power = getvar_savefile(datafiles[file_ind], 'power')
    power[where(power lt yrange[0]/100.)] = yrange[0]/100.  ; Remove negative values for log plot

    plot_x = []
    plot_y = []
    for datapoint = 0,n_elements(power)-1 do begin
      plot_y = [plot_y, power[datapoint], power[datapoint]]
      plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
    endfor
    if file_ind eq 0 then overplot=0 else overplot=1
    cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
      linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.25,$
      ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
      xstyle=8 ;draw only the main axis, don't draw the top axis
  endfor
  cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.25
  xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
  ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
  cgText, xlocation, ylocation, 'baseline length (wavelengths)', $
    /Normal, Alignment=0.5, Charsize=1.25
  cglegend, title=['Data', 'Model'], $
    linestyle=[0,2], thick=8, $
    color=[colors], length=0.03, /center_sym, location=[.7,.8], charsize=1, /box, background='white'
  cgps_close, /png, /delete_ps, density=800


end