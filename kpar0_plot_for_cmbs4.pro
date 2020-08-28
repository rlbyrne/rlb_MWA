pro kpar0_plot_for_cmbs4

  path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/'
  legend_title = ['Data', 'Point Source Model']
  

  pol = 'xx'

  yrange = [1.e11, 3.e15]
  datafiles = [path+'fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_dencorr_k0power.idlsave', $
    path+'fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
    ;path+'fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
  colors = ['black', 'blue']
  linestyles = [0,0]

  cgps_open, '/Users/rubybyrne/diffuse_survey_plotting_for_CMBS4/pt_source_model.png'
  cgDisplay, 900, 650
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
      linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.,$
      ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
      xstyle=8 ;draw only the main axis, don't draw the top axis
  endfor
  cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
  xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
  ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
  cgText, xlocation, ylocation+.03, 'baseline length (wavelengths)', $
    /Normal, Alignment=0.5, Charsize=1.
  cglegend, title=legend_title, $
    linestyle=linestyles, thick=8, $
    color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=2.
  cgControl, Resize=[800,800]
  cgps_close, /png, /delete_ps, density=800



end