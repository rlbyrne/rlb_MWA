pro plot_weights

  pathname = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_nb_2013longrun/ps/'
  pol = 'XX'
  term = 'even'
  obsid_1 = '1061316176'
  obsid_2 = '1061316296'
  
  filename1_tukey = pathname + '1061316176_'+term+'_cube'+pol+'_tk_weights_uvf.idlsave'
  filename2_tukey = pathname + 'Combined_obs_1061316296_'+term+'_cube'+pol+'_tk_weights_uvf.idlsave'
  filename1_tophat = pathname + 'Combined_obs_1061316176_'+term+'_cube'+pol+'_weights_uvf.idlsave'
  filename2_tophat = pathname + 'Combined_obs_1061316296_'+term+'_cube'+pol+'_weights_uvf.idlsave'
  
  tukey_weights_1 = getvar_savefile(filename1_tukey, "weights_cube")
  tukey_weights_2 = getvar_savefile(filename2_tukey, "weights_cube")
  tophat_weights_1 = getvar_savefile(filename1_tophat, "weights_cube")
  tophat_weights_2 = getvar_savefile(filename2_tophat, "weights_cube")
  
  tukey_weights_1_mean = MEAN(abs(tukey_weights_1), DIMENSION = 3, /NAN)
  tukey_weights_2_mean = MEAN(abs(tukey_weights_2), DIMENSION = 3, /NAN)
  tophat_weights_1_mean = MEAN(abs(tophat_weights_1), DIMENSION = 3, /NAN)
  tophat_weights_2_mean = MEAN(abs(tophat_weights_2), DIMENSION = 3, /NAN)
  
  note_tukey_1 = obsid_1+' , tukey, pol '+pol+', '+term
  note_tukey_2 = obsid_2+' , tukey, pol '+pol+', '+term
  note_tophat_1 = obsid_1+' , no tukey, pol '+pol+', '+term
  note_tophat_2 = obsid_2+' , no tukey, pol '+pol+', '+term
  
  output_tukey_1 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/weights_'+obsid_1+'_tukey_'+pol+'_'+term
  output_tukey_2 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/weights_'+obsid_2+'_tukey_'+pol+'_'+term
  output_tophat_1 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/weights_'+obsid_1+'_tophat_'+pol+'_'+term
  output_tophat_2 = '/nfs/eor-00/h1/rbyrne/error_analysis_plots/weights_'+obsid_2+'_tophat_'+pol+'_'+term
  
  QUICK_IMAGE, tukey_weights_1_mean, $
    DATA_ASPECT = .5,$
    TITLE = 'Average Weights', $
    XTITLE = 'kx', YTITLE = 'ky',$
    note = note_tukey_1, $
    data_range = [-1e6,1e6], $
    savefile = output_tukey_1
    
  QUICK_IMAGE, tukey_weights_2_mean, $
    DATA_ASPECT = .5,$
    TITLE = 'Average Weights', $
    XTITLE = 'kx', YTITLE = 'ky',$
    note = note_tukey_2, $
    data_range = [-1e6,1e6], $
    savefile = output_tukey_2
    
    
  QUICK_IMAGE, tophat_weights_1_mean, $
    DATA_ASPECT = .5,$
    TITLE = 'Average Weights', $
    XTITLE = 'kx', YTITLE = 'ky',$
    note = note_tophat_1, $
    data_range = [-1e6,1e6], $
    savefile = output_tophat_1
    
  QUICK_IMAGE, tophat_weights_2_mean, $
    DATA_ASPECT = .5,$
    TITLE = 'Average Weights', $
    XTITLE = 'kx', YTITLE = 'ky',$
    note = note_tophat_2, $
    data_range = [-1e6,1e6], $
    savefile = output_tophat_2
    
  stop
  
end