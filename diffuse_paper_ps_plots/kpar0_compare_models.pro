pro kpar0_compare_models

  ;version_names = ['fhd_rlb_GLEAM_calibration_reference_Aug2020', 'fhd_rlb_subtract_StokesI_diffuse_and_GLEAM_Aug2020', $
  ;  'fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020', 'fhd_rlb_subtract_diffuse_only_Aug2020']
  version_names = ['fhd_rlb_GLEAM_calibration_reference_Aug2020', $
    'fhd_rlb_subtract_diffuse_only_Aug2020', 'fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020']
  path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/'
  pols = ['xx', 'yy']
  yrange = [1e11, 3.e15]
  xrange=[2e-3, 2e-1]
  split_xloc_wl = 10.
  bl_range = [4.25194, 50.]

  colors = ['black', 'blue', 'turquoise', 'black']
  linestyles = [1,0,0,2]
  linewidths = [6,6,6,6]
  legend_labels = ['Compact Model', 'Diffuse Model', 'Compact and Diffuse Models', 'Data']

  for pol_ind=0,n_elements(pols)-1 do begin
    pol = pols[pol_ind]

    datafiles = []
    for run_ind=0,n_elements(version_names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
    datafiles = [datafiles, path+'fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']

    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Sept2020/models_'+pol+'.png'
    cgDisplay, 900, 650
    plotting_order = [1, 2, 0, 3]
    for ind = 0,n_elements(datafiles)-1 do begin
      file_ind = plotting_order[ind]
      k_edges = getvar_savefile(datafiles[file_ind], 'k_edges')
      power = getvar_savefile(datafiles[file_ind], 'power')
      power[where(power lt yrange[0]/100.)] = yrange[0]/100.  ; Remove negative values for log plot

      plot_x = []
      plot_y = []
      for datapoint = 0,n_elements(power)-1 do begin
        plot_y = [plot_y, power[datapoint], power[datapoint]]
        plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
      endfor
      if ind eq 0 then begin
        cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=linewidths[file_ind], title='', Charsize=1.5,$
          ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
          xstyle=8, /nodata
          cgcolorfill, [xrange[0], bl_range[0]*1e-3, bl_range[0]*1e-3, xrange[0]], [yrange[0], yrange[0], yrange[1], yrange[1]], $
          color='BLK2'
        cgcolorfill, [xrange[1], bl_range[1]*1e-3, bl_range[1]*1e-3, xrange[1]], [yrange[0], yrange[0], yrange[1], yrange[1]], $
          color='BLK2'
      endif
      
      cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
        linestyle=linestyles[file_ind], color=colors[file_ind], thick=linewidths[file_ind], /overplot, title='', Charsize=1.5,$
        ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=8 ;draw only the main axis, don't draw the top axis
    endfor
    ; Draw and redraw axes
    cgAxis, xaxis=0, xrange=xrange, xstyle=1, xtitle=textoidl(''), Charsize=1.5
    cgaxis, yaxis=0, yrange=yrange, ystyle=1, ytitle=textoidl(''), Charsize=1.5
    cgaxis, yaxis=1, yrange=yrange, ystyle=1, ytitle=textoidl(''), Charsize=1.5, yTICKFORMAT="(A1)"
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.5
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.5
    legend_ordering = [3,0,1,2]
    cglegend, title=legend_labels[legend_ordering], $
      linestyle=linestyles[legend_ordering], thick=6, $
      color=colors[legend_ordering], length=0.03, /center_sym, location=[.5,.87], charsize=1.2, /box, background='white', vspace=1.5
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800
    
  endfor

end