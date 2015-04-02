;This program plots a rebinned psd from input LC
;provide lc file name and psd index in command line, 
;LC.dat has 2 columns: time, flux;;;, flux error
;change the value of binwidth below correspondingly 

COMPILE_OPT idl2, HIDDEN

;arg 0: power law index of PSD
;arg 1: how long is the flare in unit of time bins

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

;fname = strtrim(STRMID(args[0],0,strlen(args)-4),1)
beta = double(args[0])

m=300.
rms=60.
gAmp=rms*5
N=4096
; have to specify bin width in the same unit as LC time
;binwidth=(600.D)/(86400.D)
binwidth=(50.D)/(86400.D)

gSig=uint(args[1])

funcSimLCbool=simLC_func(N,m,rms,beta,binwidth,x)
Vflux=real_part(x)
VLtime=findgen(N)*binwidth

gauss_peak=gAmp*exp(-(findgen(N)-N/2)^2/(2*gsig^2))

VfluxFlare=Vflux+gauss_peak

set_plot,'ps'
device,filename='simLCandPSD_beta'+string(beta,format='(f0.2)')+'_m'+strtrim(string(m,format='(i)'),1)+'_u'+strtrim(string(rms,format='(i)'),1)+'_5rmsPeakWidth'+strtrim(string(gSig,format='(i)'),1)+'bins.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

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

plot,VLtime,VfluxFlare,ystyle=3,/xstyle,psym=1,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;

oplot,VLtime,Vflux,color=Orange,psym=1
oplot,VLtime,VfluxFlare,color=Blue,psym=4

xyouts,max(VLtime)*0.75,(max(VfluxFlare)-min(VfluxFlare))*0.9+min(Vflux),'with flare ',font=0,color=Blue
xyouts,max(VLtime)*0.75,(max(VfluxFlare)-min(VfluxFlare))*0.85+min(Vflux),'sim ',font=0,color=Orange

!p.position=[0.15,0.08,0.95,0.48]

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freq=dblarr(N/2)
psd1=dblarr(N/2)
freqFlare=dblarr(N/2)
psdFlare=dblarr(N/2)

psdbool=PSD_fft_func(VLtime,Vflux,sqrt(Vflux),binwidth,freq,psd1)
psdFlarebool=PSD_fft_func(VLtime,VfluxFlare,sqrt(VfluxFlare),binwidth,freqFlare,psdFlare)

  binfactor = 10.
  m = long(1.0D* (N/2L) / binfactor)
  mFlare = long(1.0D* (N/2L) / binfactor)

  nf = m*binfactor
  nfFlare = mFlare*binfactor
  p = psd1[0:nf-1]
  f = freq[0:nf-1]
  pA = psdFlare[0:nfFlare-1]
  fA = freqFlare[0:nfFlare-1]

  psd = rebin(p,m)
  ;frq = rebin(f,m)
  psdA = rebin(pA,mFlare)
  ;frqA = rebin(fA,mFlare)
  frq = 10^rebin(alog10(f),m)
  frqA = 10^rebin(alog10(fA),mFlare)

  err = sqrt(1.0/binfactor)*psd
  errA = sqrt(1.0/binfactor)*psdA
  errf = sqrt(1.0/binfactor)*frq
  errfA = sqrt(1.0/binfactor)*frqA

  plot,frqA,psdA,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(frq)*0.5,max(frq)*2],yrange=[min(psdA)*5e-2,max(psdA)*2.],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'
  oploterror,frq,psd,errf,err,color=Orange,errcolor=Orange,/nohat,psym=8
  oploterror,frqA,psdA,errfA,errA,color=Blue,errcolor=Blue,/nohat,psym=8
  psdFlarefit=LINFIT(alog10(freqFlare),alog10(psdFlare),CHISQ=psdFlareChisq)
  psdfit=LINFIT(alog10(freq),alog10(psd1),CHISQ=psdutChisq)
  oplot,frqA,10^((psdFlarefit[0])+(psdFlarefit[1])*alog10(frqA)),thick=5,linestyle=3,color=Blue
  oplot,frq,10^((psdfit[0])+(psdfit[1])*alog10(frq)),thick=5,linestyle=3,color=Orange

  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psd)*1.,max(psd)*1.],'with flare:'+string(psdFlarefit[1], format='(f0.3)'),font=0,color=Blue
  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psd)*0.6,max(psd)*0.6],'sim:'+string(psdfit[1], format='(f0.3)'),font=0,color=Orange


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

