pro kpar0_plot_for_diffuse_paper

  names = ['GLEAM_only', 'GLEAM_and_diffuse_I', 'GLEAM_and_diffuse', 'diffuse_only']
  version_names = ['fhd_rlb_GLEAM_calibration_reference_Aug2020', 'fhd_rlb_subtract_StokesI_diffuse_and_GLEAM_Aug2020', $
    'fhd_rlb_subtract_diffuse_and_GLEAM_Aug2020', 'fhd_rlb_subtract_diffuse_only_Aug2020']
  path = '/Volumes/Bilbo/rlb_fhd_outputs/diffuse_survey_Aug2020/'
  pols = ['xx', 'yy']
  yrange = [1.e11, 3.e15]
  
  plot_by_versions = 1
  plot_all_models = 1
  plot_all_residuals = 1
  plot_ratios = 0
  plot_ratios_liny = 0
  use_log_binned = 1
  
  for pol_ind=0,n_elements(pols)-1 do begin
    pol = pols[pol_ind]
       
    if keyword_set(plot_by_versions) then begin
      for run_ind =0,n_elements(names)-1 do begin
        
        if keyword_set(use_log_binned) then begin
          datafiles = [path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave', $
            path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave', $
            path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
        endif else begin
          datafiles = [path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_dencorr_k0power.idlsave', $
            path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_dencorr_k0power.idlsave', $
            path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
        endelse
        
        colors = ['blue', 'red', 'orange']
        linestyles = [0,0,0]
      
        cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/'+names[run_ind]+'_'+pol+'.png'
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
        cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
          /Normal, Alignment=0.5, Charsize=1.
        cglegend, title=['Data', 'Model', 'Res'], $
          linestyle=linestyles, thick=8, $
          color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
        cgControl, Resize=[800,800]
        cgps_close, /png, /delete_ps, density=800
        
      endfor
    endif
    
    if keyword_set(plot_all_models) then begin
      
      datafiles = []
      if keyword_set(use_log_binned) then begin
        for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
      endif else begin
        for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_model_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
      endelse
      
      colors = ['black', 'blue', 'red', 'orange']
      linestyles = [0,0,0,0]
      
      cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/models_'+pol+'.png'
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
      cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
        /Normal, Alignment=0.5, Charsize=1.
      cglegend, title=['GLEAM', 'GLEAM and diffuse I', 'GLEAM and diffuse', 'diffuse only'], $
        linestyle=linestyles, thick=8, $
        color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
      cgControl, Resize=[800,800]
      cgps_close, /png, /delete_ps, density=800
      
    endif
    
    if keyword_set(plot_all_residuals) then begin
  
      datafiles = []
      if keyword_set(use_log_binned) then begin
        for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_logkperp_dencorr_k0power.idlsave']
      endif else begin
        for run_ind=0,n_elements(names)-1 do datafiles = [datafiles, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
      endelse
  
      colors = ['black', 'blue', 'red', 'orange']
      linestyles = [0,0,0,0]
  
      cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/residuals_'+pol+'.png'
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
      cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
        /Normal, Alignment=0.5, Charsize=1.
      cglegend, title=['GLEAM', 'GLEAM and diffuse I', 'GLEAM and diffuse', 'diffuse only'], $
        linestyle=linestyles, thick=8, $
        color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
      cgControl, Resize=[800,800]
      cgps_close, /png, /delete_ps, density=800
  
    endif
    
    if keyword_set(plot_ratios) then begin
    
      dirty_data = []
      res_data = []
      for run_ind=0,n_elements(names)-1 do begin
        dirty_data = [dirty_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
        res_data = [res_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
      endfor
      
      colors = ['black', 'blue', 'red', 'orange']
      linestyles = [0,0,0,0]
      yrange = [.001, 4.]

      cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/res_ratios_'+pol+'.png'
      cgDisplay, 900, 650
      for file_ind = 0,n_elements(dirty_data)-1 do begin
        k_edges = getvar_savefile(dirty_data[file_ind], 'k_edges')
        xrange=[8e-4, 2e-1]
        power_dirty = getvar_savefile(dirty_data[file_ind], 'power')
        power_res = getvar_savefile(res_data[file_ind], 'power')
        power = power_res/power_dirty
        power[where(power_dirty eq 0.)] = yrange[1]*2
        power[where(power lt yrange[0])] = yrange[0]/2

        plot_x = []
        plot_y = []
        for datapoint = 0,n_elements(power)-1 do begin
          plot_y = [plot_y, power[datapoint], power[datapoint]]
          plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
        endfor
        if file_ind eq 0 then overplot=0 else overplot=1
        cgplot, plot_x, plot_y, /xlog, /ylog, yrange=yrange, xrange=xrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.,$
          ytitle=textoidl('k-parallel=0 Residual Power/Dirty Power'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
          xstyle=8 ;draw only the main axis, don't draw the top axis
      endfor
      cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
      xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
      ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
      cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
        /Normal, Alignment=0.5, Charsize=1.
      cglegend, title=['GLEAM', 'GLEAM and diffuse I', 'GLEAM and diffuse', 'diffuse only'], $
        linestyle=linestyles, thick=8, $
        color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
      cgControl, Resize=[800,800]
      cgps_close, /png, /delete_ps, density=800
    endif
    
    if keyword_set(plot_ratios_liny) then begin

      dirty_data = []
      res_data = []
      for run_ind=0,n_elements(names)-1 do begin
        dirty_data = [dirty_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_dirty_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
        res_data = [res_data, path+version_names[run_ind]+'/ps/data/1d_binning/1131454296_gridded_uvf__even_odd_joint_noimgclip_res_'+pol+'_averemove_swbh_dencorr_k0power.idlsave']
      endfor

      colors = ['black', 'blue', 'red', 'orange']
      linestyles = [0,0,0,0]
      yrange = [.001, 1.1]

      cgps_open, '/Users/rubybyrne/diffuse_survey_paper_plotting_Aug2020/res_ratios_liny_'+pol+'.png'
      cgDisplay, 900, 650
      for file_ind = 0,n_elements(dirty_data)-1 do begin
        k_edges = getvar_savefile(dirty_data[file_ind], 'k_edges')
        xrange=[8e-4, 2e-1]
        power_dirty = getvar_savefile(dirty_data[file_ind], 'power')
        power_res = getvar_savefile(res_data[file_ind], 'power')
        power = power_res/power_dirty
        power[where(power_dirty eq 0.)] = yrange[1]*2
        power[where(power lt yrange[0])] = yrange[0]/2

        plot_x = []
        plot_y = []
        for datapoint = 0,n_elements(power)-1 do begin
          plot_y = [plot_y, power[datapoint], power[datapoint]]
          plot_x = [plot_x, k_edges[datapoint], k_edges[datapoint+1]]
        endfor
        if file_ind eq 0 then overplot=0 else overplot=1
        cgplot, plot_x, plot_y, /xlog, yrange=yrange, xrange=xrange, $
          linestyle=linestyles[file_ind], color=colors[file_ind], thick=8, overplot=overplot, title='', Charsize=1.,$
          ytitle=textoidl('k-parallel=0 Residual Power/Dirty Power'), xtitle=textoidl('k-perpendicular (!8h!X Mpc^{-1})'), $
          xstyle=8 ;draw only the main axis, don't draw the top axis
      endfor
      cgAxis, XAxis=1.0, XRange=xrange*1.e3, XStyle=1, xtitle=textoidl(''), Charsize=1.
      xlocation = (!X.Window[1] - !X.Window[0]) / 2  + !X.Window[0]
      ylocation = !Y.Window[1] + 2.75 * (!D.Y_CH_Size / Float(!D.Y_Size))
      cgText, xlocation, ylocation+.01, 'baseline length (wavelengths)', $
        /Normal, Alignment=0.5, Charsize=1.
      cglegend, title=['GLEAM', 'GLEAM and diffuse I', 'GLEAM and diffuse', 'diffuse only'], $
        linestyle=linestyles, thick=8, $
        color=colors, length=0.03, /center_sym, location=[.7,.8], charsize=1., /box, background='white', vspace=1.
      cgControl, Resize=[800,800]
      cgps_close, /png, /delete_ps, density=800
    endif
    
  endfor

end