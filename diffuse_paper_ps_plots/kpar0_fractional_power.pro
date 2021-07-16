pro kpar0_fractional_power

  version_names = ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/fhd_rlb_GLEAM_calibration_reference_Aug2020', $
    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_subtract_StokesI_diffuse_and_GLEAM_Jul2021', $
    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_May2021/fhd_rlb_subtract_diffuse_and_GLEAM_Jul2021']
  pols = ['xx', 'yy']
  yrange = [0,100]
  ;xrange=[2e-3, 2e-1]
  xrange=[2.5e-3, 2.5e-1]
  bl_range = [6.1, 50.]

  colors = ['Violet Red', 'Cornflower Blue', 'YGB7']
  linestyles = [0,0,0]
  linewidths = [6,6,6]
  legend_labels = ['Compact Sources', 'Compact and Stokes I Diffuse', 'Compact and Polarized Diffuse']

  for pol_ind=0,n_elements(pols)-1 do begin
    pol = pols[pol_ind]
    
    dirty_file = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/fhd_rlb_GLEAM_calibration_reference_Aug2020/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave'
    dirty_power = getvar_savefile(dirty_file, 'power')
 
    datafiles = []
    for run_ind=0,n_elements(version_names)-1 do datafiles = [datafiles, version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']

    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Jul2021/frac_power_recovered_'+pol+'.png'
    cgDisplay, 900, 650
    for file_ind = 0,n_elements(datafiles)-1 do begin
      k_edges = getvar_savefile(datafiles[file_ind], 'k_edges')
      power = getvar_savefile(datafiles[file_ind], 'power')
      frac_power = (1-sqrt(power/dirty_power))*100.

      plot_x = []
      plot_y = []
      for datapoint = 0,n_elements(frac_power)-1 do begin
        plot_y = [plot_y, frac_power[datapoint], frac_power[datapoint]]
        plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
      endfor
      if file_ind eq 0 then begin
        cgplot, plot_x, plot_y, /xlog, yrange=yrange, xrange=xrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=linewidths[file_ind], title='', Charsize=1.5,$
          ytitle=textoidl('Fraction of Signal Modeled (%)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
          xstyle=4, /nodata, Position=[0.1, 0.22, 0.97, 0.9]
          cgcolorfill, [xrange[0], bl_range[0]*1e-3, bl_range[0]*1e-3, xrange[0]], [yrange[0], yrange[0], yrange[1], yrange[1]], $
            color='BLK2'
          cgcolorfill, [xrange[1], bl_range[1]*1e-3, bl_range[1]*1e-3, xrange[1]], [yrange[0], yrange[0], yrange[1], yrange[1]], $
            color='BLK2'
      endif
      cgplot, plot_x, plot_y, /xlog, yrange=yrange, xrange=xrange, $
        linestyle=linestyles[file_ind], color=colors[file_ind], thick=linewidths[file_ind], /overplot, title='', Charsize=1.5,$
        ytitle=textoidl('Fraction of Signal Modeled (%)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=4 ;suppress horizontal axes
    endfor
    ; Draw lower x axis
    tick_angles_to_label = [0.4, 0.6, 0.8, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30, 40, 60]
    tick_angles = [.1*findgen(9, start=1), indgen(9, start=1), 10*indgen(9, start=1), 100, 180]
    tick_angles = reverse(tick_angles)
    tick_names = []
    for tick_ind=0,n_elements(tick_angles)-1 do begin
      null = where(tick_angles_to_label eq tick_angles[tick_ind], count)
      if count eq 0 then begin
        tick_names = [tick_names, ' ']
      endif else begin
        if tick_angles[tick_ind] ge 1 then use_string=strtrim(fix(tick_angles[tick_ind]), 1) $
        else use_string=STRING(tick_angles[tick_ind], FORMAT='(F3.1)')
        tick_names = [tick_names, use_string]
      endelse
    endfor
    tick_pos = 1/sin(tick_angles/180.*!Pi)
    cgAxis, 0.1, 0.1, /normal, xAxis=0, /Save, Color='black', Title='angular scale (degrees)', xRange=xrange*1.e3, xstyle=1, Charsize=1.5, xlog=1,$
      xtickv=tick_pos, xticks=n_elements(tick_angles), xtickname=tick_names

    cgAxis, xaxis=1, xrange=xrange, xstyle=1, xtitle=textoidl(''), Charsize=1.5 ;draw top axis
    cgaxis, yaxis=1, yrange=yrange, ystyle=1, ytitle=textoidl(''), Charsize=1.5, yTICKFORMAT="(A1)" ;draw right axis
    cgAxis, XAxis=0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl('baseline length (wavelengths)'), Charsize=1.5 ;draw bottom axis
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
      /Normal, Alignment=0.5, Charsize=1.5
    cglegend, title=legend_labels, $
      linestyle=linestyles, thick=6, $
      color=colors, length=0.03, /center_sym, location=[.5,.4], charsize=1.1, /box, background='white', vspace=1.5
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800

  endfor

end