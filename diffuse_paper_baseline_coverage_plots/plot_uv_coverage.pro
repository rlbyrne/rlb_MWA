pro plot_uv_coverage

  weights = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/fhd_rlb_get_uvf_cubes_Apr2020/1061316296_even_gridded_uvf.sav', 'weights_uv_arr')
  obs = getvar_savefile('/Volumes/Bilbo/rlb_fhd_outputs/fhd_rlb_get_uvf_cubes_Apr2020/metadata/1061316296_obs.sav', 'obs')
  pol = 0 ;use XX pol
  freq_res = obs.freq_res*2
  freq_center = obs.freq_center
  ;weights = weights[*,95:96] ;use only the center frequency
  n_freqs = (size(weights))[2]
  freq_center_ind = floor(n_freqs/2)
  uv_pix_size = 0.5 ;assume half-wavelength resolution
  n_xvals = (size(*weights[pol, freq_center_ind]))[1]
  n_yvals = (size(*weights[pol, freq_center_ind]))[2]

  xvals = (indgen(n_xvals)-n_xvals/2)*uv_pix_size
  yvals = (indgen(n_yvals)-n_yvals/2)*uv_pix_size

  plot_cube = make_array(n_xvals, n_yvals, n_freqs, value=0.)

  for freq_ind = 0, n_freqs-1 do begin
    plot_cube[*,*,freq_ind] = abs(*weights[pol, freq_ind])
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
  xvals_mesh = make_array(n_elements(xvals), n_elements(yvals), value=0.)
  for ind=0,n_elements(yvals)-1 do begin
    xvals_mesh[*, ind] = xvals
  endfor
  yvals_mesh = make_array(n_elements(xvals), n_elements(yvals), value=0.)
  for ind=0,n_elements(xvals)-1 do begin
    yvals_mesh[ind, *] = yvals
  endfor
  distances_from_origin = sqrt(xvals_mesh^2.+yvals_mesh^2.)
  
  quick_image, plot_image, xvals, yvals, data_range=[0,1], $
    xtitle='U (wavelengths)', ytitle = 'V (wavelengths)', cb_title='UV Weight',$
    savefile='/Users/rubybyrne/plot_uv_coverate_May2020/uv_weights_log', /png, /log

  quick_image, yes_no_image, xvals, yvals, data_range=[0,1], $
    xtitle='U (wavelengths)', ytitle = 'V (wavelengths)', cb_title='UV Weight',$
    savefile='/Users/rubybyrne/plot_uv_coverate_May2020/uv_nonzero', /png

end