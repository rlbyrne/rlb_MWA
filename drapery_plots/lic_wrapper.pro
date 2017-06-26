pro lic_wrapper

  ;create the vector field
  vx = make_array(500,500,value=1.)
  vx[0:100,200:300]=4
  vy = make_array(500,500,value=1.)
  vx[200:300,200:300]=10
  vy[200:300,200:300]=10
  
  ;create the line integral convolution map:
  niter=10
  len=20
  map = lic(vx,vy,niter=niter,len=len)
  
  ;plot:
  cgps_open, '/nfs/eor-00/h1/rbyrne/drapery_plotting/drapery_test_plot.png'
  cgloadct,rgb_table=palette ;create the grayscale color scheme
  map_scaled = cgscalevector(map, 0, 255, minvalue=min(map), maxvalue=max(map))
  cgimage, map_scaled, palette=palette, color='black', /keep_aspect_ratio
  ;cgcolorbar, range=[min(map),max(map)], palette=palette, color='black'
  cgps_close, /png, /delete_ps
  
end
