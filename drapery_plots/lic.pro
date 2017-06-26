;+
; NAME:
;       line integral convolution LIC
;
; PURPOSE:
;
;       This program generates a two dimensional scalar field
;       out from a vector field, representing texture of streamlines.
;
; CATEGORY:
;
;       Graphics
;
; CALLING SEQUENCE:
;
;       map = lic(vec1,vec2,niter=niter,length=length)
;
; OPTIONAL INPUTS:
;
;       /norm - to normalize the scalar field, showing only direction of vectors.
;
; OUTPUTS:
;
;       map - the 2-dimensional texture map that represents the field lines.
;
; INPUT PARAMETERS:
;
;       niter - number of interations (minimum value is 1)
;
;       vec1, vec2 - 2-dimensional arrays representing each component of the vector field
;
;       length - total maximum length to be integrated per interation
;
;******************************************************************************************;
;  Copyright (c) 2009, by Grzegorz Kowal and Diego Falceta Goncalves                       ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use allowed given that the above copyright                           ;
;        notice is retained                                                                ;
;******************************************************************************************;



function lic, vx, vy, length = len, niter = ni, normalize = normalize, amplitude = amplitude, level = level, scalar = scalar

  if ( n_elements(len) eq 0 ) then len = 8
  if ( n_elements(ni)  eq 0 ) then ni  = 1

  sz = size(vx, /dim)

  nx = sz[0]
  ny = sz[1]

  uu = sqrt(vx^2 + vy^2)
  ii = where(uu eq 0.0)
  if (ii[0] ne -1) then uu[ii] = 1.0
  if keyword_set(normalize) then begin
    ux = vx / uu
    uy = vy / uu
  endif else begin
    ux = vx / max(uu)
    uy = vy / max(uu)
  endelse

  dp = [0.0, 0.0]
  dm = [0.0, 0.0]

  vl = randomu(1, sz)

  for it = 0, ni-1 do begin

    texture = vl

    vv = fltarr(sz)
    pi = rebin(reform(indgen(sz[0]),sz[0],1),sz[0],sz[1])
    pj = rebin(reform(indgen(sz[1]),1,sz[1]),sz[0],sz[1])
    mi = pi
    mj = pj

    ppi = 1.*pi
    ppj = 1.*pj
    mmi = 1.*mi
    mmj = 1.*mj

    for l = 0, len do begin
      dpi = interpolate(ux, ppi, ppj)
      dpj = interpolate(uy, ppi, ppj)
      dmi = interpolate(ux, mmi, mmj)
      dmj = interpolate(uy, mmi, mmj)

      ppi = ppi + 0.25*dpi
      ppj = ppj + 0.25*dpj
      mmi = mmi - 0.25*dmi
      mmj = mmj - 0.25*dmj

      pi = (round(ppi) + sz[0]) mod sz[0]
      pj = (round(ppj) + sz[1]) mod sz[1]
      mi = (round(mmi) + sz[0]) mod sz[0]
      mj = (round(mmj) + sz[1]) mod sz[1]

      ppi = pi + (ppi - round(ppi))
      ppj = pj + (ppj - round(ppj))
      mmi = mi + (mmi - round(mmi))
      mmj = mj + (mmj - round(mmj))

      vv = vv + interpolate(texture, ppi, ppj) + interpolate(texture, mmi, mmj)
    endfor

    vl = 0.25 * vv / len

  endfor

  if keyword_set(amplitude) then begin
    if (n_elements(scalar) gt 1) then uu = scalar
    if (n_elements(level) ne 1) then level = 0.1
    if (ii[0] ne -1) then uu[ii] = 0.0
    level = max([0.0, min([level, 1.0])])
    uu = (1.0 - level) * uu / max(uu) + level
    vl = vl*uu
  endif
  if (ii[0] ne -1) then vl[ii] = 0.0

  return, vl
end
