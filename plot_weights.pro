pro plot_weights

  weights = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/fhd_rlb_get_uvf_cubes_Apr2020/1061316296_even_gridded_uvf.sav', 'weights_uv_arr')
  obs = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/fhd_rlb_get_uvf_cubes_Apr2020/metadata/1061316296_obs.sav', 'obs')
  pol = 0 ;use XX pol
  freq_res = obs.freq_res*2
  freq_center = obs.freq_center
  freq_center_ind = (size(weights))[2]/2
  uv_pix_size = 0.5 ;assume half-wavelength resolution
  n_xvals = (size(*weights[pol, freq_center_ind]))[1]
  n_yvals = (size(*weights[pol, freq_center_ind]))[2]
  
  n_freqs_use = [1,5,9,101,191]
  for n_freq_ind = 0,n_elements(n_freqs_use)-1 do begin
    n_freqs = n_freqs_use[n_freq_ind]
    n_freqs = floor((n_freqs-1)/2)*2+1 ;ensure n_freqs is odd
    freq_interval = n_freqs*freq_res
    xvals = (indgen(n_xvals)-n_xvals/2)*uv_pix_size
    yvals = (indgen(n_yvals)-n_yvals/2)*uv_pix_size
    
    plot_cube = make_array(n_xvals, n_yvals, n_freqs, value=0.)
    
    for freq_ind = 0, n_freqs-1 do begin
      freq_channel = freq_center_ind-(n_freqs-1)/2+freq_ind
      plot_cube[*,*,freq_ind] = abs(*weights[pol, freq_channel])
    endfor
    
    if n_freqs gt 1 then plot_image = mean(plot_cube, dimension=3) else plot_image = plot_cube
    yes_no_image = make_array(n_xvals, n_yvals, value=1)
    yes_no_image[where(plot_image eq 0)] = 0
    
    xrange = [-50,50]
    yrange = [-50,50]
    plot_image = plot_image[xrange[0]/uv_pix_size+n_xvals/2:xrange[1]/uv_pix_size+n_xvals/2, yrange[0]/uv_pix_size+n_yvals/2:yrange[1]/uv_pix_size+n_yvals/2]
    yes_no_image = yes_no_image[xrange[0]/uv_pix_size+n_xvals/2:xrange[1]/uv_pix_size+n_xvals/2, yrange[0]/uv_pix_size+n_yvals/2:yrange[1]/uv_pix_size+n_yvals/2]
    xvals = xvals[xrange[0]/uv_pix_size+n_xvals/2:xrange[1]/uv_pix_size+n_xvals/2]
    yvals = yvals[yrange[0]/uv_pix_size+n_yvals/2:yrange[1]/uv_pix_size+n_yvals/2]
    
    quick_image, plot_image, xvals, yvals, data_range=[0,1], $
      xtitle='UV (wavelengths)', ytitle = 'UV (wavelengths)', title = 'UV Weights Over '+string(freq_interval/1e6)+' MHz', $
      savefile='/Users/rubybyrne/plot_weights_Apr2020/uv_weights_'+strsplit(string(n_freqs),/extract)+'freq', /png
      
    quick_image, yes_no_image, xvals, yvals, data_range=[0,1], $
      xtitle='UV (wavelengths)', ytitle = 'UV (wavelengths)', title = 'UV Weights Over '+string(freq_interval/1e6)+' MHz', $
      savefile='/Users/rubybyrne/plot_weights_Apr2020/uv_nonzero_'+strsplit(string(n_freqs),/extract)+'freq', /png
  endfor

end