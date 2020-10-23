pro kpar0_rms_explore_plot


  paths = ['/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020', $
    '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/fhd_rlb_subtract_StokesI_diffuse_and_GLEAM_Aug2020']
  angles = ['45', '90', '135', '180', '225', '270', '315']
  for angle_ind=0,n_elements(angles)-1 do begin
    paths = [paths, '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/rotate_diffuse_Aug2020/fhd_rlb_subtract_rotated_diffuse_and_GLEAM_angle_'+angles[angle_ind]+'_Aug2020']
  endfor
  pols = ['xx', 'yy']
  yrange = [1.e11, 3.e15]

  for pol_ind=0,n_elements(pols)-1 do begin
    pol = pols[pol_ind]

    datafiles = [paths[0]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_dencorr_k0power.idlsave', $
      paths[1]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave',$
      paths[0]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
    legend_titles = ['Dirty', 'Unpolarized Res.', 'Rot. Ang. 0']
    for angle_ind=0,n_elements(angles)-1 do begin
      datafiles = [paths[angle_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave', datafiles]
      legend_titles = ['Rot. Ang. '+angles[angle_ind], legend_titles]
    endfor
    
    colors = ['orange', 'blue', 'cyan', 'purple', 'pur8', 'green','magenta', 'black', 'black', 'red']
    linestyles = [make_array(n_elements(angles), value=0), 0, 1, 0]

    cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/rm_exploration_'+pol+'.png'
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
        linestyle=linestyles[file_ind], $
        color=colors[file_ind], $
        thick=8, overplot=overplot, title='', Charsize=1.,$
        ytitle=textoidl('k-parallel=0 Power (mK^2 !8h!X^{-3} Mpc^3)'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
        xstyle=8 ;draw only the main axis, don't draw the top axis
    endfor
    cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
    xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
    ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
    cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
      /Normal, Alignment=0.5, Charsize=1.
    cglegend, title=legend_titles, $
      linestyle=linestyles, thick=8, $
      color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
    cgControl, Resize=[800,800]
    cgps_close, /png, /delete_ps, density=800


  endfor

end