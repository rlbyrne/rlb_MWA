FUNCTION split_vis_flags_vissig_freq,obs,flag_arr,bi_use=bi_use,preserve_flags=preserve_flags,even_only=even_only,odd_only=odd_only
;function to create split even and odd time samples with identical flagging


n_pol=N_Elements(flag_arr) ;number of polarizations
IF Keyword_Set(preserve_flags) THEN flag_arr_use=pointer_copy(flag_arr) ELSE flag_arr_use=flag_arr
bin_start=(*obs.baseline_info).bin_offset
nt=obs.n_time

stop
IF nt LT 2 THEN RETURN,flag_arr_use
nb=(size(*flag_arr_use[0],/dimension))[1] ;number of bins
bin_end=fltarr(nt)
bin_end[0:nt-2]=bin_start[1:nt-1]-1
bin_end[nt-1]=nb-1
bin_i=lonarr(nb)-1
nt2=Floor(nt/2)
FOR t_i=0,2*nt2-1 DO bin_i[bin_start[t_i]:bin_end[t_i]]=t_i

stop

time_use=(*obs.baseline_info).time_use
time_start_i=Min(where(time_use))
nt3=Floor((nt-time_start_i)/2)*2
time_use_0=time_use[time_start_i:time_start_i+nt3-1:2]
time_use_1=time_use[time_start_i+1:time_start_i+nt3-1:2]
time_use_01=time_use_0*time_use_1
time_use*=0
time_use[time_start_i:time_start_i+nt3-1:2]=time_use_01
time_use[time_start_i+1:time_start_i+nt3-1:2]=time_use_01
time_cut_i=where(time_use LE 0,nt_cut)
stop
IF nt_cut GT 0 THEN BEGIN
    bin_i_cut=where(bin_i EQ time_cut_i,n_cut)
    IF n_cut GT 0 THEN bin_i[bin_i_cut]=-1 ; will be skipped by using where(bin_i mod 2 EQ 0,1) below (-1 mod 2 is still -1)
ENDIF
stop

bi_use=Ptrarr(2,/allocate)
*bi_use[0]=where(bin_i mod 2 EQ 0,n_even)
*bi_use[1]=where(bin_i mod 2 EQ 1,n_odd)
stop

IF n_even LT n_odd THEN *bi_use[1]=(*bi_use[1])[0:n_even-1]
IF n_odd LT n_even THEN *bi_use[0]=(*bi_use[0])[0:n_odd-1]

FOR pol_i=0,n_pol-1 DO BEGIN
    flag_use0=0>(*flag_arr_use[pol_i])[*,*bi_use[0]]<(*flag_arr_use[pol_i])[*,*bi_use[1]]<1
    *flag_arr_use[pol_i]*=0
    IF ~Keyword_Set(odd_only) THEN (*flag_arr_use[pol_i])[*,*bi_use[0]]=flag_use0
    IF ~Keyword_Set(even_only) THEN (*flag_arr_use[pol_i])[*,*bi_use[1]]=flag_use0 
    flag_use0=0 ;free memory
ENDFOR
RETURN,flag_arr_use
END