pro plot_diffuse_cal_kpar0

  path = '/Volumes/Bilbo/rlb_fhd_outputs/calibration_tests_Feb2021'
  cubes = ['Aug23_all', 'Aug23_minustwo', 'Aug23_minusone', 'Aug23_zenith', 'Aug23_plusone', 'Aug23_plustwo']
  ;cubes = ['Aug23_all', 'Aug23_minustwo', 'Aug23_zenith', 'Aug23_plusone', 'Aug23_plustwo']
  cube_type = 'model'
  
  ;version_names = ['fhd_rlb_perfreq_cal_Feb2021', 'fhd_rlb_polarized_diffuse_with_amp_corrected_calibration_May2021']
  version_names = ['fhd_rlb_perfreq_cal_Feb2021', 'fhd_rlb_perfreq_all_baseline_cal_with_polarized_diffuse_Mar2021']

  pols = ['xx','yy']
  yrange = [1e14, 1e17]
  xrange=[2.5e-3, 1e-1]
  bl_range = [6.1, 50.]

  colors = ['Violet Red', 'Cornflower Blue', 'Violet Red', 'Cornflower Blue']
  linestyles = [0,0,2,2]
  linewidths = [6,6,6,6]
  legend_labels = ['GLEAM cal XX', 'Diffuse and GLEAM cal XX', 'GLEAM cal YY', 'Diffuse and GLEAM cal YY']

  for cube_ind=0,n_elements(cubes)-1 do begin
    cube = cubes[cube_ind]

    datafiles = []
    for pol_ind=0,n_elements(pols)-1 do begin
      pol = pols[pol_ind]
      for version_ind=0,n_elements(version_names)-1 do begin
        datafiles = [datafiles, path+'/'+version_names[version_ind]+'/ps/data/1d_binning/Combined_obs_'+cube+'_cubeXX__even_odd_joint_noimgclip_'+cube_type+'_'+pol+'_dft_averemove_swbh_dencorr_k0power.idlsave']
      endfor
    endfor
  
    cgps_open, '/Users/rubybyrne/diffuse_calibration_test_plotting_May2021/'+cube_type+'_'+cube+'_reference_run.png'
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
      if file_ind eq 0 then begin
        cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=linewidths[file_ind], title='', Charsize=1.5,$
          ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
          xstyle=8, /nodata, Position=[0.1, 0.22, 0.97, 0.9]
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
  
    cgAxis, xaxis=0, xrange=xrange, xstyle=1, xtitle=textoidl(''), Charsize=1.5
    cgaxis, yaxis=0, yrange=yrange, ystyle=1, ytitle=textoidl(''), Charsize=1.5
    cgaxis, yaxis=1, yrange=yrange, ystyle=1, ytitle=textoidl(''), Charsize=1.5, yTICKFORMAT="(A1)"
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.5
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.5
    cglegend, title=legend_labels, $
      linestyle=linestyles, thick=6, $
      color=colors, length=0.03, /center_sym, location=[.15,.8], charsize=1.1, /box, background='white', vspace=1.5
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800
  
  endfor

end