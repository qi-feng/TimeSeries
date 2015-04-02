;This program plots a rebinned psd from input LC
;provide lc file name and psd index in command line, 
;LC.dat has 2 columns: time, flux;;;, flux error
;change the value of binwidth below correspondingly 

COMPILE_OPT idl2, HIDDEN

;arg 0: power law index of PSD
;arg 1: how much longer

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

;fname = strtrim(STRMID(args[0],0,strlen(args)-4),1)
beta = double(args[0])

m=300.
rms=60.
N=1024
Nlong=N*uint(args[1])

; have to specify bin width in the same unit as LC time
;binwidth=(600.D)/(86400.D)
binwidth=(50.D)/(86400.D)


funcSimLCbool=simLC_func(Nlong,m,rms,beta,binwidth,x)
Vfluxlong=real_part(x)
VLtimelong=findgen(Nlong)*binwidth

Vflux=Vfluxlong[0:N-1]
VLtime=VLtimelong[0:N-1]

set_plot,'ps'
device,filename='simLCandPSD_beta'+string(beta,format='(f0.2)')+'_m'+strtrim(string(m,format='(i)'),1)+'_u'+strtrim(string(rms,format='(i)'),1)+'_RNL'+strtrim(string(args[1],format='(i)'),1)+'.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

device, helvetica=1
device, isolatin1=1
;set_plot,'x'

device, decomposed=1
index=100;

Red=GETCOLOR('red',/true)
Blue=GETCOLOR('blue',/true)
Green=GETCOLOR('green',/true)
Cyan=GETCOLOR('cyan',/true)
Magenta=GETCOLOR('magenta',/true)
Yellow=GETCOLOR('yellow',/true)
Pink=GETCOLOR('pink',/true)
Orange=GETCOLOR('orange',/true)
Gray=GETCOLOR('gray',/true)
Navy=GETCOLOR('navy',/true)

psize=0.6
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;
;
;
!p.multi=[0,1,2]
!p.position=[0.15,0.57,0.95,0.97] 

plot,VLtimelong,Vfluxlong,ystyle=3,/xstyle,psym=1,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;

oplot,VLtimelong,Vfluxlong,color=Orange,psym=1
oplot,VLtime,Vflux,color=Blue,psym=4

xyouts,[max(VLtimelong)*0.75,max(VLtimelong)*0.75],[(max(Vflux)-min(Vflux))*0.9+min(Vflux),min(Vflux)+(max(Vflux)-min(Vflux))*0.9],'short ',font=0,color=Blue
xyouts,[max(VLtimelong)*0.75,max(VLtimelong)*0.75],[(max(Vflux)-min(Vflux))*0.85+min(Vflux),min(Vflux)+(max(Vflux)-min(Vflux))*0.85],'long ',font=0,color=Orange

!p.position=[0.15,0.08,0.95,0.48]

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freq=dblarr(N/2)
psd1=dblarr(N/2)
freqlong=dblarr(Nlong/2)
psdlong=dblarr(Nlong/2)

psdbool=PSD_fft_func(VLtime,Vflux,sqrt(Vflux),binwidth,freq,psd1)
psdlongbool=PSD_fft_func(VLtimelong,Vfluxlong,sqrt(Vfluxlong),binwidth,freqlong,psdlong)

  binfactor = 10.
  m = long(1.0D* (N/2L) / binfactor)
  mlong = long(1.0D* (Nlong/2L) / binfactor)

  nf = m*binfactor
  nflong = mlong*binfactor
  p = psd1[0:nf-1]
  f = freq[0:nf-1]
  pA = psdlong[0:nflong-1]
  fA = freqlong[0:nflong-1]

  psd = rebin(p,m)
  ;frq = rebin(f,m)
  psdA = rebin(pA,mlong)
  ;frqA = rebin(fA,mAlias)
  frq = 10^rebin(alog10(f),m)
  frqA = 10^rebin(alog10(fA),mlong)

  err = sqrt(1.0/binfactor)*psd
  errA = sqrt(1.0/binfactor)*psdA
  errf = sqrt(1.0/binfactor)*frq
  errfA = sqrt(1.0/binfactor)*frqA

  plot,frqA,psdA,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(frqA)*0.5,max(frq)*2],yrange=[min(psdA)*5e-2,max(psd)*4.],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'
  oploterror,frq,psd,errf,err,color=Blue,errcolor=Blue,/nohat,psym=8
  oploterror,frqA,psdA,errfA,errA,color=Orange,errcolor=Orange,/nohat,psym=8
  psdlongfit=LINFIT(alog10(freqlong),alog10(psdlong),CHISQ=psdlongChisq)
  psdfit=LINFIT(alog10(freq),alog10(psd1),CHISQ=psdutChisq)
  oplot,frqA,10^((psdlongfit[0])+(psdlongfit[1])*alog10(frqA)),thick=5,linestyle=3,color=Orange
  oplot,frq,10^((psdfit[0])+(psdfit[1])*alog10(frq)),thick=5,linestyle=3,color=Blue

  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psdA)*1.,max(psdA)*1.],'short:'+string(psdlongfit[1], format='(f0.3)'),font=0,color=Blue
  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psdA)*0.4,max(psdA)*0.4],'long:'+string(psdfit[1], format='(f0.3)'),font=0,color=Orange


;1. Read a real LC, careful with unit in LC (assuming days)
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
;readcol,args[0],F='D,D',VLtime,Vflux;
;readcol,'simLC_beta1_m300u60_6.dat',F='D,D',VLtime,Vflux;
;readcol,'pn_all_0.5_1keV_50s.txt',F='D,D,D',VLtime,Vflux,dVflux;
;readcol,'out_xmm0501_0.5_10keV.dat',F='D,D,D',VLtime,Vflux,dVflux;

;readcol,'NightlyVERLC.txt',F='D,D,D',VLtime,Vflux,dVflux;

;VLtime=VLtime/86400.D
;binwidth=1.
; duration of real LC
;tspan=0L
;binratio=0; so that sim LC has smaller binwidth/2^3 
;NV=n_elements(VLtime)
;tspan=(VLtime[NV-1]-VLtime[0])/binwidth

;RESOLVE_ROUTINE, 'suf_plotpsd_func', /IS_FUNCTION
;RESOLVE_ROUTINE, 'powspec', /IS_FUNCTION

;sufoutstruct=suf_plotpsd_func(VLtime, Vflux, binwidth, binratio,beta,3 )
;powspecoutstruct=powspec(VLtime*86400.D, Vflux)


;Nout=n_elements(sufoutstruct.realFreqrebin)
;outPSD=dblarr(3,Nout);
;outPSD[0,*]=sufoutstruct.realFreqrebin
;outPSD[1,*]=sufoutstruct.realPSDrebin
;outPSD[2,*]=sufoutstruct.drealPSDrebin
;
;openw,lun,arg_arr[0]+'_PSD.dat',/get_lun
;printf,lun,outPSD
;Free_lun,lun


device, /close_file

end

