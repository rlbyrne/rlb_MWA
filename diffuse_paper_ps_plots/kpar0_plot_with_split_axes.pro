pro kpar0_plot_with_split_axes

  names = ['GLEAM_only', 'GLEAM_and_diffuse_I', 'GLEAM_and_diffuse']
  version_names = ['fhd_rlb_GLEAM_calibration_reference_Aug2020', 'fhd_rlb_subtract_StokesI_diffuse_and_GLEAM_Aug2020', $
    'fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020']
  path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/'
  pols = ['xx', 'yy']
  yrange = [1.e10, 3.e15]
  xrange=[8e-4, 2e-1]
  split_xloc_wl = 10.

  colors = ['black', 'black', 'blue', 'red']
  linestyles = [2, 0,0,0]
  legend_labels = ['Data', 'Compact Model', 'Compact and Unpol. Diffuse Model', 'Compact and Diffuse Model']

  for pol_ind=0,n_elements(pols)-1 do begin
    pol = pols[pol_ind]

    ;;;Plot models;;;

    datafiles = [path+'fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
    for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']

    ;plot normally
    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Sept2020/models_'+pol+'.png'
    cgDisplay, 900, 650
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
      cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
        linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.,$
        ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=8 ;draw only the main axis, don't draw the top axis
    endfor
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.
    cglegend, title=legend_labels, $
      linestyle=linestyles, thick=8, $
      color=colors, length=0.03, /center_sym, location=[.5,.87], charsize=1., /box, background='white', vspace=1.
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800
    
    if 0 then begin; Plot split
      cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Sept2020/models_'+pol+'_split.png'
      cgDisplay, 900, 650
  
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
  
        cgplot, plot_x, plot_y, Position=[0.15, 0.15, 0.4, 0.9], xStyle=9, YStyle=9, $
          xrange=[xrange[0], split_xloc_wl/1.e3], /ylog, yrange=yrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1., $
          ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)')
      endfor
      cgAxis, XAxis=1.0, XStyle=1, xtitle=textoidl(''), xrange=[xrange[0]*1.e3, split_xloc_wl], Charsize=1.
  
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
  
        cgplot, plot_x, plot_y, Position=[0.4, 0.15, 0.9, 0.9], xStyle=9, YStyle=5, $
          xrange=[split_xloc_wl/1.e3, xrange[1]], /xlog, /ylog, yrange=yrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1., /noerase
      endfor 
      cgAxis, XAxis=1.0, XStyle=1, xtitle=textoidl(''), xrange=[split_xloc_wl, xrange[1]*1.e3], /xlog, Charsize=1.
      cgaxis, yaxis=1, ytickformat='(A1)', ystyle=1
  
      xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
      ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
      cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
        /Normal, Alignment=0.5, Charsize=1.
      cglegend, title=legend_labels, $
        linestyle=linestyles, thick=8, $
        color=colors, length=0.03, /center_sym, location=[.5,.87], charsize=1., /box, background='white', vspace=1.
      cgControl, Resize=[800,800]
      cgps_close, /png, /delete_ps, density=800
    endif

    ;;;Plot residuals;;;

    datafiles = [path+'fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
    for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']

    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Sept2020/residuals_'+pol+'.png'
    cgDisplay, 900, 650
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
      cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
        linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.,$
        ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=8 ;draw only the main axis, don't draw the top axis
    endfor
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.
    cglegend, title=legend_labels, $
      linestyle=linestyles, thick=8, $
      color=colors, length=0.03, /center_sym, location=[.5,.87], charsize=1., /box, background='white', vspace=1.
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800


    ;;;Plot ratios;;;

    dirty_data = []
    res_data = []
    for run_ind=0,n_elements(names)-1 do begin
      dirty_data = [dirty_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
      res_data = [res_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
    endfor

    yrange_res = [0, 1]

    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Sept2020/res_ratios_'+pol+'.png'
    cgDisplay, 900, 650
    for file_ind = 0,n_elements(dirty_data)-1 do begin
      k_edges = getvar_savefile(dirty_data[file_ind], 'k_edges')
      power_dirty = getvar_savefile(dirty_data[file_ind], 'power')
      power_res = getvar_savefile(res_data[file_ind], 'power')
      power = power_res/power_dirty
      power[where(power_dirty eq 0.)] = yrange_res[1]*2
      power[where(power lt yrange_res[0])] = yrange_res[0]/2

      plot_x = []
      plot_y = []
      for datapoint = 0,n_elements(power)-1 do begin
        plot_y = [plot_y, power[datapoint], power[datapoint]]
        plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
      endfor
      if file_ind eq 0 then overplot=0 else overplot=1
      cgplot, plot_x, plot_y, /xlog, yrange=yrange_res, xrange=xrange, $
        linestyle=linestyles[file_ind+1], color=colors[file_ind+1], thick=8, overplot=overplot, title='', Charsize=1.,$
        ytitle=textoidl('k-parallel=0 Residual Power/Dirty Power'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=8 ;draw only the main axis, don't draw the top axis
    endfor
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.
    cglegend, title=legend_labels[1:*], $
      linestyle=linestyles[1:*], thick=8, $
      color=colors[1:*], length=0.03, /center_sym, location=[.5,.87], charsize=1., /box, background='white', vspace=1.
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800
    

  endfor

end