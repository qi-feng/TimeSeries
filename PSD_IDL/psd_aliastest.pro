;This program plots a rebinned psd from input LC
;provide lc file name and psd index in command line, 
;LC.dat has 2 columns: time, flux;;;, flux error
;change the value of binwidth below correspondingly 

;e.g.:
;idl -e ".run psd_aliastest.pro" -args 1.0 8

COMPILE_OPT idl2, HIDDEN

;arg 0: power law index of PSD
;arg 1: alias ratio

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

;fname = strtrim(STRMID(args[0],0,strlen(args)-4),1)
beta = double(args[0])

m=300.
rms=60.
N=4096
; have to specify bin width in the same unit as LC time
;binwidth=(600.D)/(86400.D)
binwidth=(50.D)/(86400.D)

aliasRatio=uint(args[1])

funcSimLCbool=simLC_func(N,m,rms,beta,binwidth,x)
Vflux=real_part(x)
VLtime=findgen(N)*binwidth
VLtimeAlias=findgen(N/aliasRatio)*binwidth*aliasRatio
VfluxAlias=dblarr(N/aliasRatio)

for j=0,N/aliasRatio-1,1 do begin
  VfluxAlias[j]=Vflux[j*aliasRatio]
endfor

;data=dblarr(2,N)
;data[0,*]=VLtime
;data[1,*]=Vflux
;openw,lun,'simLC_N'+strtrim(string(N,format='(i)'),1)+'_aliasRatio'+strtrim(string(aliasRatio,format='(i)'),1)+'_beta'+strtrim(string(beta,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'.dat',/get_lun
;printf,lun,data
;Free_lun,lun

data=dblarr(2,N/aliasRatio)
data[0,*]=VLtimeAlias
data[1,*]=VfluxAlias
openw,lun,'simLC_N'+strtrim(string(N,format='(i)'),1)+'_with_alias_'+strtrim(string(aliasRatio,format='(i)'),1)+'_beta'+strtrim(string(beta,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'.dat',/get_lun
printf,lun,data
Free_lun,lun

;print,fname,beta
;fname='simLC_N1024_beta1.00m300u60'
;fname='simLC_N1024_beta1.10m300u30'

;readcol,args[0],F='D,D',VLtime,Vflux;
;beta=2.

set_plot,'ps'
device,filename='simLCandPSD_beta'+string(beta,format='(f0.2)')+'_m'+strtrim(string(m,format='(i)'),1)+'_u'+strtrim(string(rms,format='(i)'),1)+'_alias'+strtrim(string(aliasRatio,format='(i)'),1)+'.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

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

plot,VLtime,Vflux,ystyle=3,/xstyle,psym=1,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;

oplot,VLtime,Vflux,color=Orange,psym=1
oplot,VLtimeAlias,VfluxAlias,color=Blue,psym=4

xyouts,[max(VLtime)*0.75,max(VLtime)*0.75],[(max(Vflux)-min(Vflux))*0.9+min(Vflux),min(Vflux)+(max(Vflux)-min(Vflux))*0.9],'with aliasing ',font=0,color=Blue
xyouts,[max(VLtime)*0.75,max(VLtime)*0.75],[(max(Vflux)-min(Vflux))*0.85+min(Vflux),min(Vflux)+(max(Vflux)-min(Vflux))*0.85],'sim ',font=0,color=Orange

!p.position=[0.15,0.08,0.95,0.48]

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freq=dblarr(N/2)
psd1=dblarr(N/2)
freqAlias=dblarr(N/aliasRatio/2)
psdAlias=dblarr(N/aliasRatio/2)

psdbool=PSD_fft_func(VLtime,Vflux,sqrt(Vflux),binwidth,freq,psd1)
psdAliasbool=PSD_fft_func(VLtimeAlias,VfluxAlias,sqrt(VfluxAlias),aliasRatio*binwidth,freqAlias,psdAlias)

  binfactor = 10.
  m = long(1.0D* (N/2L) / binfactor)
  mAlias = long(1.0D* (N/aliasRatio/2L) / binfactor)

  nf = m*binfactor
  nfAlias = mAlias*binfactor
  p = psd1[0:nf-1]
  f = freq[0:nf-1]
  pA = psdAlias[0:nfAlias-1]
  fA = freqAlias[0:nfAlias-1]

  psd = rebin(p,m)
  ;frq = rebin(f,m)
  psdA = rebin(pA,mAlias)
  ;frqA = rebin(fA,mAlias)
  frq = 10^rebin(alog10(f),m)
  frqA = 10^rebin(alog10(fA),mAlias)

  err = sqrt(1.0/binfactor)*psd
  errA = sqrt(1.0/binfactor)*psdA
  errf = sqrt(1.0/binfactor)*frq
  errfA = sqrt(1.0/binfactor)*frqA

  plot,frqA,psdA,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(frq)*0.5,max(frq)*2],yrange=[min(psdA)*5e-2,max(psd)*2.],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'
  oploterror,frq,psd,errf,err,color=Orange,errcolor=Orange,/nohat,psym=8
  oploterror,frqA,psdA,errfA,errA,color=Blue,errcolor=Blue,/nohat,psym=8
  psdAliasfit=LINFIT(alog10(freqAlias),alog10(psdAlias),CHISQ=psdAliasChisq)
  psdfit=LINFIT(alog10(freq),alog10(psd1),CHISQ=psdutChisq)
  oplot,frqA,10^((psdAliasfit[0])+(psdAliasfit[1])*alog10(frqA)),thick=5,linestyle=3,color=Blue
  oplot,frq,10^((psdfit[0])+(psdfit[1])*alog10(frq)),thick=5,linestyle=3,color=Orange

  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psdA)*1.,max(psdA)*1.],'with aliasing:'+string(psdAliasfit[1], format='(f0.3)'),font=0,color=Blue
  xyouts,[max(frqA)*0.1,max(frqA)*0.1],[max(psdA)*0.6,max(psdA)*0.6],'sim:'+string(psdfit[1], format='(f0.3)'),font=0,color=Orange


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

