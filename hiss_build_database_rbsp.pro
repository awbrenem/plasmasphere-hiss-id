;-Verify that it is hiss
;-Plot dN/N and dB/B. Maybe the hiss needs this to actually grow!!!
;-Plot MLT
;-Remove data at L<2?
;-Use Scott's Lshell spectra to get averaged plasmapause position
;-plot cyclotron resonance energies




;-Hiss file (involves only EFW and basic EMFISIS data):
;--Using plasmasphere file, save file for each day with hiss data:
;---For each time DURING ENTIRE DAY (need these times for particle
;---data interpolation)
;--Make sure to have one for EACH day. Especially b/c I don't want to
;--load a gigantic file for BARREL comparision
;-----freq of peak hiss power
;-----f/flh
;-----f/100Hz
;-----f/fce
;-----f/fci
;-----amplitude
;-----E/B (at peak?)


pro hiss_build_database_rbsp



  ;--------------------------------------------------
  ;SCRIPT INPUT

  args = command_line_args()

  print,'-----------------'
  print,'DATE = ',args[0]
  print,'SC = ',args[1]
  print,'Enterps = ',args[2]
  print,'Exitps = ',args[3]
  print,'----------------'


  date = args[0]
  probe = args[1]
  enterps = args[2]
  exitps = args[3]


  tplot_options,'xmargin',[20.,16.]
  tplot_options,'ymargin',[3,9]
  tplot_options,'xticklen',0.08
  tplot_options,'yticklen',0.02
  tplot_options,'xthick',2
  tplot_options,'ythick',2
  tplot_options,'labflag',-1



  ;Must change enterps and exitps from single strings to arrays
  ;enterps = "['2013-01-01/01:16:03', '2013-01-01/06:08:26', '2013-01-01/14:44:27']"
  enterps = strmid(enterps,2)
  exitps = strmid(exitps,2)

  goo = strsplit(enterps,',',/extract,escape=']')
  for j=0,n_elements(goo)-1 do $
  if j eq 0 then goo[j] = strmid(goo[j],0,19) else goo[j] = strmid(goo[j],2,19)
  boo = strsplit(exitps,',',/extract,escape=']')
  for j=0,n_elements(boo)-1 do $
  if j eq 0 then boo[j] = strmid(boo[j],0,19) else boo[j] = strmid(boo[j],2,19)


  enterps = goo
  exitps = boo

  print,'ENTERPS vals = ' + enterps
  help,enterps
  print,'----****----'



  ;; ;--------------------------------------------------
  ;; ;TEST input


  ;;   fn = 'plasmasphere_rbspa_database.txt'
  ;;   path = '~/Desktop/code/Aaron/RBSP/survey_programs_hiss/'

  ;;   openr,lun,path+fn,/get_lun
  ;;   junk = ''
  ;;   tt0 = ''
  ;;   tt1 = ''

  ;;   while not eof(lun) do begin
  ;;      readf,lun,junk
  ;;      tt0 = [tt0,strmid(junk,0,19)]
  ;;      tt1 = [tt1,strmid(junk,20,19)]
  ;;   endwhile

  ;;   tt0 = tt0[1:n_elements(tt0)-1]
  ;;   tt1 = tt1[1:n_elements(tt1)-1]

  ;;   for i=0,n_elements(tt0)-1 do print,tt0[i] + ' to ' + tt1[i]

  ;;   date = '2013-01-03'           ;BARREL event
  ;;   t0 = date + '/00:00'
  ;;   t1 = date + '/24:00'
  ;;   probe = 'a'


  ;;   tt0tmp = strmid(tt0,0,10)
  ;;   tt1tmp = strmid(tt1,0,10)

  ;;   goo = where(tt0tmp eq date)
  ;;   if goo[0] ne -1 then Enterps = tt0[goo]
  ;;   goo = where(tt1tmp eq date)
  ;;   if goo[0] ne -1 then Exitps = tt1[goo]


  ;;   close,lun
  ;;   free_lun,lun

  ;; ;--------------------------------------------------



  timespan,date
  rbspx = 'rbsp'+probe

  rbsp_load_efw_spec,type='calibrated',/pT,probe=probe
  rbsp_load_emfisis,probe=probe,coord='gse',cadence='4sec',level='l3'
  rbsp_load_efw_waveform_l3,probe=probe


  get_data,rbspx+'_emfisis_l3_4sec_gse_Magnitude',data=magnit
  fce = 28.*magnit.y
  store_data,'fce',data={x:magnit.x,y:fce}
  store_data,'fce_2',data={x:magnit.x,y:fce/2.}
  store_data,'fci',data={x:magnit.x,y:fce/1836.}
  store_data,'flh',data={x:magnit.x,y:sqrt(fce*fce/1836.)}


  get_data,rbspx+'_efw_64_spec0',data=dat
  get_data,rbspx +'_efw_64_spec0',data=Evals,dlimits=dlim,limits=lim



  timespan,date

  ;THESE ARE THE REFERENCE TIMES
  times = Evals.x


  tinterpol_mxn,rbspx+'_efw_flags_all',times,newname=rbspx+'_efw_flags_all'
  tinterpol_mxn,['fce','fce_2','fci','flh'],times,$
    newname=['fce','fce_2','fci','flh']


  split_vec,rbspx+'_efw_flags_all'


  copy_data,rbspx+'_efw_flags_all_0',rbspx+'_efw_flag_global'
  copy_data,rbspx+'_efw_flags_all_1',rbspx+'_efw_flag_eclipse'
  copy_data,rbspx+'_efw_flags_all_15',rbspx+'_efw_flag_charging'



  ;make sure units are in nT
  get_data,rbspx +'_efw_64_spec4',data=Bvals,dlimits=dlimB
  Bvals.y = Bvals.y/1000.
  E2B = Evals.y/Bvals.y         ;Velocity=1000*E/B (km/s)
  boo = where(finite(E2B) eq 0)
  if boo[0] ne -1 then E2B[boo] = !values.f_nan

  store_data,'E2B',data={x:Evals.x,y:E2B,v:Evals.v},dlimits=dlim,limits=lim
  store_data,'E2B_C',data=['E2B','fce','fce_2','fci','flh']
  store_data,rbspx+'_efw_64_spec0_C',$
    data=[rbspx+'_efw_64_spec0','fce','fce_2','fci','flh']
  store_data,rbspx+'_efw_64_spec2_C',$
    data=[rbspx+'_efw_64_spec2','fce','fce_2','fci','flh']
  store_data,rbspx+'_efw_64_spec3_C',$
    data=[rbspx+'_efw_64_spec3','fce','fce_2','fci','flh']
  store_data,rbspx+'_efw_64_spec4_C',$
    data=[rbspx+'_efw_64_spec4','fce','fce_2','fci','flh']
  ylim,rbspx+'_efw_64_spec?_C',10,1000,1
  ylim,'E2B_C',3,10000,1
  zlim,'E2B_C',0.001,1d2,1
  options,'E2B_C','ztitle','(mV/m / nT)!C!CE/B'
  options,'E2B_C','ytitle','E12/Bw!C[Hz]'
  options,'E2B_C','ysubtitle',''




  ;Create a tplot variable of plasmasphere start and end times

  ps = replicate(0.,n_elements(times))

  for i=0,n_elements(enterps)-1 do begin
    goo = where((times ge time_double(enterps[i])) and $
                (times le time_double(exitps[i])))

    if goo[0] ne -1 then ps[goo] = 1
  endfor

  store_data,'plasmasphere',data={x:times,y:ps}
  store_data,'plasmasphere_C',data=['plasmasphere',rbspx+'_efw_Vavg']
  options,'plasmasphere_C','colors',[250,0]
  ylim,'plasmasphere_C',-10,5



  ;NOW REDUCE THE SPECTRAL DATA TO POSSIBLE HISS
  get_data,rbspx+'_efw_64_spec0',data=spec0,dlim=sdlime0,lim=slime0
  get_data,rbspx+'_efw_64_spec2',data=spec2,dlim=sdlimbu,lim=slimbu
  get_data,rbspx+'_efw_64_spec3',data=spec3,dlim=sdlimbv,lim=slimbv
  get_data,rbspx+'_efw_64_spec4',data=spec4,dlim=sdlimbw,lim=slimbw

  ;Remove lowest 3 freq bins (up to 20 Hz)
  spec0.y[*,0:2] = !values.f_nan
  spec2.y[*,0:2] = !values.f_nan
  spec3.y[*,0:2] = !values.f_nan
  spec4.y[*,0:2] = !values.f_nan

  ;Remove spec data outside of plasmasphere
  get_data,'plasmasphere',data=ps
  goo = where(ps.y eq 0.)
  if goo[0] ne -1 then spec0.y[goo,*] = !values.f_nan
  if goo[0] ne -1 then spec2.y[goo,*] = !values.f_nan
  if goo[0] ne -1 then spec3.y[goo,*] = !values.f_nan
  if goo[0] ne -1 then spec4.y[goo,*] = !values.f_nan


  ;Remove data during eclipse times
  get_data,rbspx+'_efw_flag_eclipse',data=flageclipse
  boo = where(flageclipse.y eq 1)
  if boo[0] ne -1 then spec0.y[boo,*] = !values.f_nan
  if boo[0] ne -1 then spec2.y[boo,*] = !values.f_nan
  if boo[0] ne -1 then spec3.y[boo,*] = !values.f_nan
  if boo[0] ne -1 then spec4.y[boo,*] = !values.f_nan



  ;; ;Remove flagged data
  ;; get_data,rbspx+'_efw_flag_global',data=flagglobal
  ;; boo = where(flagglobal.y eq 1)
  ;; if boo[0] ne -1 then spec0.y[boo,*] = !values.f_nan
  ;; if boo[0] ne -1 then spec2.y[boo,*] = !values.f_nan
  ;; if boo[0] ne -1 then spec3.y[boo,*] = !values.f_nan
  ;; if boo[0] ne -1 then spec4.y[boo,*] = !values.f_nan

  ;Remove power below threshold minimum amplitude
  ampthresh_E = 1e-06
  ampthresh_b = 1.
  boo = where(spec0.y lt ampthresh_E)
  if boo[0] ne -1 then spec0.y[boo] = !values.f_nan
  boo = where(spec2.y lt ampthresh_B)
  if boo[0] ne -1 then spec2.y[boo] = !values.f_nan
  boo = where(spec3.y lt ampthresh_B)
  if boo[0] ne -1 then spec3.y[boo] = !values.f_nan
  boo = where(spec4.y lt ampthresh_B)
  if boo[0] ne -1 then spec4.y[boo] = !values.f_nan



  ;plasmasphere versions only
  store_data,rbspx+'_efw_64_spec0_ps',data=spec0,dlim=sdlime0,lim=slime0
  store_data,rbspx+'_efw_64_spec2_ps',data=spec2,dlim=sdlimbu,lim=slimbu
  store_data,rbspx+'_efw_64_spec3_ps',data=spec3,dlim=sdlimbv,lim=slimbv
  store_data,rbspx+'_efw_64_spec4_ps',data=spec4,dlim=sdlimbw,lim=slimbw
  options,rbspx+'_efw_64_spec?_ps','spec',1


  ;Track the peak hiss frequency
  peak_e0 = fltarr(n_elements(spec0.x))
  peak_bu = peak_e0 & peak_bv = peak_e0 & peak_bw = peak_e0



  for i=0L,n_elements(spec0.x)-1 do begin
    goo = max(spec0.y[i,*],wh,/nan)
    if goo[0] ne -1 then peak_e0[i] = spec0.v[wh]
  endfor
  for i=0L,n_elements(spec2.x)-1 do begin
    goo = max(spec2.y[i,*],wh,/nan)
    if goo[0] ne -1 then peak_bu[i] = spec2.v[wh]
  endfor
  for i=0L,n_elements(spec3.x)-1 do begin
    goo = max(spec3.y[i,*],wh,/nan)
    if goo[0] ne -1 then peak_bv[i] = spec3.v[wh]
  endfor
  for i=0L,n_elements(spec4.x)-1 do begin
    goo = max(spec4.y[i,*],wh,/nan)
    if goo[0] ne -1 then peak_bw[i] = spec4.v[wh]
  endfor


  boo = where(peak_e0 eq 4.)
  if boo[0] ne -1 then peak_e0[boo] = !values.f_nan
  boo = where(peak_bu eq 4.)
  if boo[0] ne -1 then peak_bu[boo] = !values.f_nan
  boo = where(peak_bv eq 4.)
  if boo[0] ne -1 then peak_bv[boo] = !values.f_nan
  boo = where(peak_bw eq 4.)
  if boo[0] ne -1 then peak_bw[boo] = !values.f_nan

  store_data,'peak_e0',data={x:spec0.x,y:peak_e0}
  store_data,'peak_bu',data={x:spec2.x,y:peak_bu}
  store_data,'peak_bv',data={x:spec3.x,y:peak_bv}
  store_data,'peak_bw',data={x:spec4.x,y:peak_bw}
  rbsp_detrend,'peak_' + ['e0','bu','bv','bw'],60.*1.



  store_data,rbspx+'_efw_64_spec0_ps_C',$
  data=[rbspx+'_efw_64_spec0_ps','fce','fce_2','fci','flh','peak_e0_smoothed']
  store_data,rbspx+'_efw_64_spec2_ps_C',$
  data=[rbspx+'_efw_64_spec2_ps','fce','fce_2','fci','flh','peak_bu_smoothed']
  store_data,rbspx+'_efw_64_spec3_ps_C',$
  data=[rbspx+'_efw_64_spec3_ps','fce','fce_2','fci','flh','peak_bv_smoothed']
  store_data,rbspx+'_efw_64_spec4_ps_C',$
  data=[rbspx+'_efw_64_spec4_ps','fce','fce_2','fci','flh','peak_bw_smoothed']

  ylim,rbspx+'_efw_64_spec?_ps_C',1,10000,1





  get_data,'flh',data=flh
  get_data,'fci',data=fci
  get_data,'fce',data=fce

  spec0L = spec0
  spec2L = spec2
  spec3L = spec3
  spec4L = spec4
  spec0U = spec0
  spec2U = spec2
  spec3U = spec3
  spec4U = spec4


  ;choose low and high frequencies for spectral plots
  ;spectral plot 1
  flow1 = fci.y
  fhig1 = flh.y
  ;spectral plot 2
  flow2 = flh.y
  fhig2 = fce.y



  for i=0l,n_elements(spec0L.x)-1 do begin
    boo = where((spec0L.v ge fhig1[i]) or (spec0L.v le flow1[i]))
    if boo[0] ne -1 then spec0L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec2L.x)-1 do begin
    boo = where((spec2L.v ge fhig1[i]) or (spec2L.v le flow1[i]))
    if boo[0] ne -1 then spec2L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec3L.x)-1 do begin
    boo = where((spec3L.v ge fhig1[i]) or (spec3L.v le flow1[i]))
    if boo[0] ne -1 then spec3L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec4L.x)-1 do begin
    boo = where((spec4L.v ge fhig1[i]) or (spec4L.v le flow1[i]))
    if boo[0] ne -1 then spec4L.y[i,boo] = !values.f_nan
  endfor

  ;Find hiss at flow2<f<fhig2
  for i=0l,n_elements(spec0U.x)-1 do begin
    boo = where((spec0U.v le flow2[i]) or (spec0U.v ge fhig2[i]))
    if boo[0] ne -1 then spec0U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec2U.x)-1 do begin
    boo = where((spec2U.v le flow2[i]) or (spec2U.v ge fhig2[i]))
    if boo[0] ne -1 then spec2U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec3U.x)-1 do begin
    boo = where((spec3U.v le flow2[i]) or (spec3U.v ge fhig2[i]))
    if boo[0] ne -1 then spec3U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec4U.x)-1 do begin
    boo = where((spec4U.v le flow2[i]) or (spec4U.v ge fhig2[i]))
    if boo[0] ne -1 then spec4U.y[i,boo] = !values.f_nan
  endfor

  store_data,rbspx+'_efw_64_spec0L',data=spec0L
  store_data,rbspx+'_efw_64_spec2L',data=spec2L
  store_data,rbspx+'_efw_64_spec3L',data=spec3L
  store_data,rbspx+'_efw_64_spec4L',data=spec4L
  store_data,rbspx+'_efw_64_spec0U',data=spec0U
  store_data,rbspx+'_efw_64_spec2U',data=spec2U
  store_data,rbspx+'_efw_64_spec3U',data=spec3U
  store_data,rbspx+'_efw_64_spec4U',data=spec4U

  options,rbspx+'_efw_64_spec??','spec',1
  ylim,rbspx+'_efw_64_spec??',1,10000,1
  zlim,rbspx+'_efw_64_spec0?',1d-6,1d-2,1
  zlim,rbspx+'_efw_64_' + ['spec2?','spec3?','spec4?'],0.01,100,1



  store_data,rbspx+'_efw_64_spec0L_C',$
    data=[rbspx+'_efw_64_spec0L','fce','fce_2','fci','flh','peak_e0_smoothed']
  store_data,rbspx+'_efw_64_spec2L_C',$
    data=[rbspx+'_efw_64_spec2L','fce','fce_2','fci','flh','peak_bu_smoothed']
  store_data,rbspx+'_efw_64_spec3L_C',$
    data=[rbspx+'_efw_64_spec3L','fce','fce_2','fci','flh','peak_bv_smoothed']
  store_data,rbspx+'_efw_64_spec4L_C',$
    data=[rbspx+'_efw_64_spec4L','fce','fce_2','fci','flh','peak_bw_smoothed']
  store_data,rbspx+'_efw_64_spec0U_C',$
    data=[rbspx+'_efw_64_spec0U','fce','fce_2','fci','flh','peak_e0_smoothed']
  store_data,rbspx+'_efw_64_spec2U_C',$
    data=[rbspx+'_efw_64_spec2U','fce','fce_2','fci','flh','peak_bu_smoothed']
  store_data,rbspx+'_efw_64_spec3U_C',$
    data=[rbspx+'_efw_64_spec3U','fce','fce_2','fci','flh','peak_bv_smoothed']
  store_data,rbspx+'_efw_64_spec4U_C',$
    data=[rbspx+'_efw_64_spec4U','fce','fce_2','fci','flh','peak_bw_smoothed']


  ylim,[rbspx+'_efw_64_spec0_C',rbspx+'_efw_64_spec0?_C',$
  rbspx+'_efw_64_spec2_C',rbspx+'_efw_64_spec2?_C',$
  rbspx+'_efw_64_spec3_C',rbspx+'_efw_64_spec3?_C',$
  rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4?_C'],10,10000,1
  ylim,'*flag*',-2,2


  spec0L = spec0
  spec2L = spec2
  spec3L = spec3
  spec4L = spec4
  spec0U = spec0
  spec2U = spec2
  spec3U = spec3
  spec4U = spec4


  ;choose low and high frequencies for spectral plots
  ;spectral plot 1
  flow1 = replicate(20.,n_elements(fci.y)) ;Wen Li's definition of "unusually low freq hiss"
  fhig1 = replicate(100.,n_elements(fci.y))
  ;spectral plot 2
  flow2 = replicate(100.,n_elements(fci.y)) ;Wen Li's definition of "unusually low freq hiss"
  fhig2 = fce.y


  for i=0l,n_elements(spec0L.x)-1 do begin
    boo = where((spec0L.v ge fhig1[i]) or (spec0L.v le flow1[i]))
    if boo[0] ne -1 then spec0L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec2L.x)-1 do begin
    boo = where((spec2L.v ge fhig1[i]) or (spec2L.v le flow1[i]))
    if boo[0] ne -1 then spec2L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec3L.x)-1 do begin
    boo = where((spec3L.v ge fhig1[i]) or (spec3L.v le flow1[i]))
    if boo[0] ne -1 then spec3L.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec4L.x)-1 do begin
    boo = where((spec4L.v ge fhig1[i]) or (spec4L.v le flow1[i]))
    if boo[0] ne -1 then spec4L.y[i,boo] = !values.f_nan
  endfor

  ;Find hiss at flow2<f<fhig2
  for i=0l,n_elements(spec0U.x)-1 do begin
    boo = where((spec0U.v le flow2[i]) or (spec0U.v ge fhig2[i]))
    if boo[0] ne -1 then spec0U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec2U.x)-1 do begin
    boo = where((spec2U.v le flow2[i]) or (spec2U.v ge fhig2[i]))
    if boo[0] ne -1 then spec2U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec3U.x)-1 do begin
    boo = where((spec3U.v le flow2[i]) or (spec3U.v ge fhig2[i]))
    if boo[0] ne -1 then spec3U.y[i,boo] = !values.f_nan
  endfor
  for i=0l,n_elements(spec4U.x)-1 do begin
    boo = where((spec4U.v le flow2[i]) or (spec4U.v ge fhig2[i]))
    if boo[0] ne -1 then spec4U.y[i,boo] = !values.f_nan
  endfor

  store_data,rbspx+'_efw_64_spec0L2',data=spec0L
  store_data,rbspx+'_efw_64_spec2L2',data=spec2L
  store_data,rbspx+'_efw_64_spec3L2',data=spec3L
  store_data,rbspx+'_efw_64_spec4L2',data=spec4L
  store_data,rbspx+'_efw_64_spec0U2',data=spec0U
  store_data,rbspx+'_efw_64_spec2U2',data=spec2U
  store_data,rbspx+'_efw_64_spec3U2',data=spec3U
  store_data,rbspx+'_efw_64_spec4U2',data=spec4U

  options,rbspx+'_efw_64_spec??2','spec',1
  ylim,rbspx+'_efw_64_spec??2',1,10000,1
  zlim,rbspx+'_efw_64_spec0?2',1d-6,1d-2,1
  zlim,rbspx+'_efw_64_spec2?2',0.01,100,1
  zlim,rbspx+'_efw_64_spec3?2',0.01,100,1
  zlim,rbspx+'_efw_64_spec4?2',0.01,100,1


  store_data,rbspx+'_efw_64_spec0L2_C',$
    data=[rbspx+'_efw_64_spec0L2','fce','fce_2','fci','flh','peak_e0_smoothed']
  store_data,rbspx+'_efw_64_spec2L2_C',$
    data=[rbspx+'_efw_64_spec2L2','fce','fce_2','fci','flh','peak_bu_smoothed']
  store_data,rbspx+'_efw_64_spec3L2_C',$
    data=[rbspx+'_efw_64_spec3L2','fce','fce_2','fci','flh','peak_bv_smoothed']
  store_data,rbspx+'_efw_64_spec4L2_C',$
    data=[rbspx+'_efw_64_spec4L2','fce','fce_2','fci','flh','peak_bw_smoothed']
  store_data,rbspx+'_efw_64_spec0U2_C',$
    data=[rbspx+'_efw_64_spec0U2','fce','fce_2','fci','flh','peak_e0_smoothed']
  store_data,rbspx+'_efw_64_spec2U2_C',$
    data=[rbspx+'_efw_64_spec2U2','fce','fce_2','fci','flh','peak_bu_smoothed']
  store_data,rbspx+'_efw_64_spec3U2_C',$
    data=[rbspx+'_efw_64_spec3U2','fce','fce_2','fci','flh','peak_bv_smoothed']
  store_data,rbspx+'_efw_64_spec4U2_C',$
    data=[rbspx+'_efw_64_spec4U2','fce','fce_2','fci','flh','peak_bw_smoothed']

  ylim,[rbspx+'_efw_64_spec0_C',rbspx+'_efw_64_spec0?2_C',$
  rbspx+'_efw_64_spec2_C',rbspx+'_efw_64_spec2?2_C',$
  rbspx+'_efw_64_spec3_C',rbspx+'_efw_64_spec3?2_C',$
  rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4?2_C'],10,10000,1









  ;Create boolean hiss variable showing when peak hiss is gt or lt the
  ;cutoff

  ;; get_data,'peak_e0_smoothed',data=peak_e0
  ;; boolgoo = replicate(!values.f_nan,n_elements(times))
  ;; cutoff = flh.y
  ;; gooL = where(peak_e0.y le cutoff)
  ;; if gooL[0] ne -1 then boolgoo[gooL] = -1
  ;; gooU = where(peak_e0.y gt cutoff)
  ;; if gooU[0] ne -1 then boolgoo[gooU] = 1
  ;; store_data,'hiss_bool_flh',data={x:times,y:boolgoo}

  ;; boolgoo = replicate(!values.f_nan,n_elements(times))
  ;; cutoff = replicate(100.,n_elements(times)) ;Hz
  ;; gooL = where(peak_e0.y le cutoff)
  ;; if gooL[0] ne -1 then boolgoo[gooL] = -1
  ;; gooU = where(peak_e0.y gt cutoff)
  ;; if gooU[0] ne -1 then boolgoo[gooU] = 1
  ;; store_data,'hiss_bool_100Hz',data={x:times,y:boolgoo}

  ;; ylim,['hiss_bool_flh','hiss_bool_100Hz'],-2,2
  ;; options,'*bool*','panel_size',0.5





  ;; tplot,['plasmasphere_C',$
  ;;        rbspx+'_efw_64_spec4_C',$
  ;;        'hiss_bool_flh','hiss_bool_100Hz',$
  ;;        'dn_n',$
  ;;                               ;rbspx+'_efw_64_spec0_C',$
  ;;                               ;rbspx+'_efw_64_spec0U_C',$
  ;;                               ;rbspx+'_efw_64_spec0L_C',$
  ;;        rbspx+'_efw_64_spec4U_C',$
  ;;        rbspx+'_efw_64_spec4L_C']



  ;stop


  ;; tplot,['plasmasphere_C',$
  ;;        rbspx+'_efw_64_spec4_C',$
  ;;        rbspx+'_efw_64_spec4U_C',$
  ;;        rbspx+'_efw_64_spec4L_C',$
  ;;        rbspx+'_efw_64_spec4U2_C',$
  ;;        rbspx+'_efw_64_spec4L2_C',$
  ;;        'hiss_bool_flh','hiss_bool_100Hz']
  ;;  stop

  ;----------------------------------------------------------------


  ;Integrate total hiss amplitude

  fcals = rbsp_efw_get_gain_results()
  fbinsL = fcals.cal_spec.freq_spec64L
  fbinsH = fcals.cal_spec.freq_spec64H
  bandw = fbinsH - fbinsL

  get_data,'rbsp'+probe+'_efw_64_spec2L',data=bu2
  get_data,'rbsp'+probe+'_efw_64_spec3L',data=bv2
  get_data,'rbsp'+probe+'_efw_64_spec4L',data=bw2
  bu2.y[*,0:2] = 0. & bv2.y[*,0:2] = 0. & bw2.y[*,0:2] = 0.
  bu2.y[*,45:63] = 0. & bv2.y[*,45:63] = 0. & bw2.y[*,45:63] = 0.

  nelem = n_elements(bu2.x)
  bt = fltarr(nelem)
  ;Add up amplitude in all 64 bins and multiply each bin by a bandwidth
  ball = bu2.y + bv2.y + bw2.y

  for j=0L,nelem-1 do bt[j] = sqrt(total(ball[j,*]*bandw,/nan)) ;RMS method (Malaspina's method)
  store_data,'Bfield_hissintL',data={x:bu2.x,y:bt}


  get_data,'rbsp'+probe+'_efw_64_spec2U',data=bu2
  get_data,'rbsp'+probe+'_efw_64_spec3U',data=bv2
  get_data,'rbsp'+probe+'_efw_64_spec4U',data=bw2
  bu2.y[*,0:2] = 0. & bv2.y[*,0:2] = 0. & bw2.y[*,0:2] = 0.
  bu2.y[*,45:63] = 0. & bv2.y[*,45:63] = 0. & bw2.y[*,45:63] = 0.

  nelem = n_elements(bu2.x)
  bt = fltarr(nelem)
  ;Add up amplitude in all 64 bins and multiply each bin by a bandwidth
  ball = bu2.y + bv2.y + bw2.y

  for j=0L,nelem-1 do bt[j] = sqrt(total(ball[j,*]*bandw,/nan)) ;RMS method (Malaspina's method)
  store_data,'Bfield_hissintU',data={x:bu2.x,y:bt}




  get_data,'rbsp'+probe+'_efw_64_spec2L2',data=bu2
  get_data,'rbsp'+probe+'_efw_64_spec3L2',data=bv2
  get_data,'rbsp'+probe+'_efw_64_spec4L2',data=bw2
  bu2.y[*,0:2] = 0. & bv2.y[*,0:2] = 0. & bw2.y[*,0:2] = 0.
  bu2.y[*,45:63] = 0. & bv2.y[*,45:63] = 0. & bw2.y[*,45:63] = 0.

  nelem = n_elements(bu2.x)
  bt = fltarr(nelem)
  ;Add up amplitude in all 64 bins and multiply each bin by a bandwidth
  ball = bu2.y + bv2.y + bw2.y

  for j=0L,nelem-1 do bt[j] = sqrt(total(ball[j,*]*bandw,/nan)) ;RMS method (Malaspina's method)
  store_data,'Bfield_hissintL2',data={x:bu2.x,y:bt}


  get_data,'rbsp'+probe+'_efw_64_spec2U2',data=bu2
  get_data,'rbsp'+probe+'_efw_64_spec3U2',data=bv2
  get_data,'rbsp'+probe+'_efw_64_spec4U2',data=bw2
  bu2.y[*,0:2] = 0. & bv2.y[*,0:2] = 0. & bw2.y[*,0:2] = 0.
  bu2.y[*,45:63] = 0. & bv2.y[*,45:63] = 0. & bw2.y[*,45:63] = 0.

  nelem = n_elements(bu2.x)
  bt = fltarr(nelem)
  ;Add up amplitude in all 64 bins and multiply each bin by a bandwidth
  ball = bu2.y + bv2.y + bw2.y

  for j=0L,nelem-1 do bt[j] = sqrt(total(ball[j,*]*bandw,/nan)) ;RMS method (Malaspina's method)
  store_data,'Bfield_hissintU2',data={x:bu2.x,y:bt}



  tplot_options,'title','from hiss_build_database.pro'


  ylim,'plasmasphere',0,1.5

  options,['Bfield_hissintU','Bfield_hissintL'],'panel_size',0.6
  options,['Bfield_hissintU2','Bfield_hissintL2'],'panel_size',0.6

  options,'Bfield_hissintU','ytitle','Bw RMS!C(nT)!C!Cf > flh'
  options,'Bfield_hissintL','ytitle','Bw RMS!C(nT)!C!Cf < flh'
  options,'Bfield_hissintU2','ytitle','Bw RMS!C(nT)!C!Cf > 100Hz'
  options,'Bfield_hissintL2','ytitle','Bw RMS!C(nT)!C!Cf < 100Hz'
  options,rbspx+'_efw_64_spec4U_C','ytitle','Bw [Hz]!C!Cf > flhr'
  options,rbspx+'_efw_64_spec4L_C','ytitle','Bw [Hz]!C!Cf < flhr'
  options,rbspx+'_efw_64_spec4U2_C','ytitle','Bw [Hz]!C!Cf > 100Hz'
  options,rbspx+'_efw_64_spec4L2_C','ytitle','Bw [Hz]!C!Cf < 100Hz'
  options,'plasmasphere_C','ytitle','(V1+V2)/2!Cwith PS!Cindicator'


  d2 = strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)

  fname = 'hiss_plot_for_RBSP'+probe+'_'+d2+'.ps'
  popen,'~/Desktop/code/Aaron/RBSP/survey_programs_hiss/' + fname

  !p.charsize = 0.7
  tplot,[rbspx+'_efw_64_spec4_C',$
    'Bfield_hissintU',rbspx+'_efw_64_spec4U_C',$
    'Bfield_hissintL',rbspx+'_efw_64_spec4L_C',$
    'Bfield_hissintU2',rbspx+'_efw_64_spec4U2_C',$
    'Bfield_hissintL2',rbspx+'_efw_64_spec4L2_C',$
    'plasmasphere_C','E2B_C' $
    ]

  pclose




  get_data,'peak_bw_smoothed',data=fbpeak
  get_data,'peak_e0_smoothed',data=fepeak
  get_data,'flh',data=flh
  get_data,'fce',data=fce
  get_data,'fci',data=fci

  get_data,'Bfield_hissintL',data=bamp_fcitoflh  ;b/t fci and flh
  get_data,'Bfield_hissintU',data=bamp_flhtofce  ;b/t flh and fce
  get_data,'Bfield_hissintL2',data=bamp_20to100  ;b/t 20 Hz and 100 Hz
  get_data,'Bfield_hissintU2',data=bamp_100tofce ;b/t 100 Hz and fce


  fn = 'hiss_record_for_RBSP'+probe+'_'+d2+'.txt'
  openw,lun,'~/Desktop/Research/RBSP_low_freq_hiss/idl/'+fn,/get_lun


  printf,lun,'from hiss_build_database.pro'
  printf,lun,'f = freq of peak hiss power'
  printf,lun,'B1 = hiss RMS amp (nT) in range fci<f<flh'
  printf,lun,'B2 = hiss RMS amp (nT) in range flh<f<fce'
  printf,lun,'B3 = hiss RMS amp (nT) in range f20<f<f100 Hz'
  printf,lun,'B4 = hiss RMS amp (nT) in range f100<f<fce'
  printf,lun,'time                  f/fce      f/flh      f/fci      f/100Hz  B1       B2       B3        B4'

  for i=0,n_elements(times)-1 do printf,lun,times[i],$
  fbpeak.y[i]/fce.y[i],$
  fbpeak.y[i]/flh.y[i],$
  fbpeak.y[i]/fci.y[i],$
  fbpeak.y[i]/100.,$
  bamp_fcitoflh.y[i],$
  bamp_flhtofce.y[i],$
  bamp_20to100.y[i],$
  bamp_100tofce.y[i],format='(e15.9,f11.6,f11.4,f11.4,f11.4,4f9.2)'



  close,lun
  free_lun,lun






  ;; get_data,'Bfield_hissintL',data=hissL
  ;; get_data,'Bfield_hissintU',data=hissU

  ;; ;Create boolian hiss variable
  ;; hissbool = replicate(0B,n_elements(times))
  ;; hiss_thresh_nt = 10.
  ;; goo = where(hissL.y ge hiss_thresh_nt)
  ;; if goo[0] ne -1 then hissbool[goo] = 1B
  ;; store_data,'hiss_boolL',data={x:times,y:hissbool}
  ;; ylim,'hiss_boolL',0,2

  ;; hissbool = replicate(0B,n_elements(times))
  ;; goo = where(hissU.y ge hiss_thresh_nt)
  ;; if goo[0] ne -1 then hissbool[goo] = 1B
  ;; store_data,'hiss_boolU',data={x:times,y:hissbool}
  ;; ylim,'hiss_boolU',0,2





  ;; tplot,['hiss_bool_flh','hiss_bool_100Hz','Bfield_hissintL','hiss_boolL',rbspx+'_efw_64_spec0U_C',$
  ;;        rbspx+'_efw_64_spec0L_C',rbspx+'_efw_64_spec4U_C',rbspx+'_efw_64_spec4L_C',$
  ;;        'plasmasphere_C',rbspx+'_efw_flag_global']
  ;; stop





  ;; ;**************************************************


  ;; ;tplot,['hiss_bool_100Hz',rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4U_C',rbspx+'_efw_64_spec4L_C']
  ;;   tplot,['hiss_bool_flh',rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4U_C',rbspx+'_efw_64_spec4L_C']
  ;;   stop

  ;; ;; get_data,rbspx+'_efw_64_spec4L',data=BvalsL
  ;; ;; get_data,rbspx+'_efw_64_spec4U',data=BvalsU
  ;;   get_data,rbspx+'_efw_64_spec4',data=Bvals
  ;; ;get_data,'hiss_bool_100Hz',data=hbool
  ;;   get_data,'hiss_bool_flh',data=hbool


  ;; ;Plot freq spec
  ;;   freqs = bvals.v
  ;;   title = t0z + ' - ' + t1z
  ;;   t0z = time_double(t0z)
  ;;   t1z = time_double(t1z)


  ;;   gooL = where(hbool.y eq -1)
  ;;   gooU = where(hbool.y eq 1)

  ;;   specL = Bvals
  ;;   specU = Bvals
  ;;   specL.y[*,*] = 0.
  ;;   specU.y[*,*] = 0.



  ;;   if gooL[0] ne -1 then specL.y[gooL,*] = Bvals.y[gooL,*]
  ;;   if gooU[0] ne -1 then specU.y[gooU,*] = Bvals.y[gooU,*]

  ;;   store_data,'specL_tmp',data=specL
  ;;   store_data,'specU_tmp',data=specU
  ;;   options,['specL_tmp','specU_tmp'],'spec',1
  ;;   ylim,['specL_tmp','specU_tmp'],10,10000,1
  ;;   zlim,['specL_tmp','specU_tmp'],0.1,100,1

  ;; ;; tplot,['specU_tmp','specL_tmp','hiss_bool_100Hz',rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4U_C',$
  ;; ;;        rbspx+'_efw_64_spec4L_C']
  ;;   tplot,['specU_tmp','specL_tmp','hiss_bool_flh',rbspx+'_efw_64_spec4_C',rbspx+'_efw_64_spec4U_C',$
  ;;          rbspx+'_efw_64_spec4L_C']
  ;;   stop





  ;;   yvalsL = tsample('specL_tmp',time_double([t0z,t1z]),times=tms)
  ;;   yvalsU = tsample('specU_tmp',time_double([t0z,t1z]),times=tms)

  ;; ;put in nT^2/Hz
  ;;   yvalsL = yvalsL/1000./1000.
  ;;   yvalsU = yvalsU/1000./1000.


  ;;   !p.multi = [0,0,2]
  ;;   plot,freqs,yvalsL[0,*],/nodata,xtitle='freq (Hz)',ytitle='Hiss amp (nT^2/Hz)',xrange=[20,2000],/xlog,/ylog,yrange=[1d-10,1d-2],title=title,xstyle=1
  ;;   for i=0,n_elements(tms)-1 do oplot,freqs,yvalsL[i,*]
  ;;   for i=0,n_elements(tms)-1 do oplot,freqs,yvalsU[i,*],color=250
  ;;   plot,freqs,yvalsL[0,*],/nodata,xtitle='freq (Hz)',ytitle='Hiss amp (nT^2/Hz)',xrange=[20,2000],/ylog,yrange=[1d-10,1d-2],title=title,xstyle=1
  ;;   for i=0,n_elements(tms)-1 do oplot,freqs,yvalsL[i,*]
  ;;   for i=0,n_elements(tms)-1 do oplot,freqs,yvalsU[i,*],color=250





end
