;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs, put it in a file
 
FUNCTION LWA_plot_func, realtime, realflux, binwidth, binratio,beta
COMPILE_OPT idl2, HIDDEN


tspan=0L
NV=long(n_elements(realtime))
tspan=(realtime[NV-1]-realtime[0])/binwidth

;3. Simulate a longer LC with finer sampling rate
;and with underlying power-law psd of index -beta
RESOLVE_ROUTINE, 'simLongLC_func', /IS_FUNCTION
simLong=simLongLC_func(realtime,realflux, beta, binwidth,binratio)

N=simLong.simN ;; N: total number of evenly-spaced points for sim time series, 
               ;; it is 2^int, longer than real x(t) with finer sampling
scatter=1.0
finebin=1.0D/(2L^binratio) ;;  time interval (bin width)
x0=simLong.simx
t=simLong.simt
meanf=mean(realflux)
;uf=sigma(realflux)
uf=stddev(realflux)

;print,"simLong.simN: ",simLong.simN," n_elements(simLong.simx): ",n_elements(simLong.simx);

;duplicate this long fine sim LC for contamination... sad...
;x0=x0*stddev(realflux-mean(realflux))/stddev(x0)
;x0=(x0-mean(x0))*stddev(realflux)/stddev(x0)+meanf
simx=x0

;3. Cut a piece of the simulated LC of the same duration as real LC,
;   and then bin the cut sim LC to same bin width as observations
;   cutting from the tail of vector for best randomness
datapadN=2L^(uint(alog(tspan)/alog(2))+1)
Ncut=datapadN*(2L^binratio)
Ncutrebin=datapadN

print,'N: ',N,'Ncut: ',Ncut,'Ncutrebin: ',Ncutrebin

tcut=dblarr(Ncut)
xcut=dblarr(Ncut)
tcut=t[0:Ncut-1]
xcut=x0[N-Ncut:N-1]
;xcut=xcut-mean(xcut)
;xcut=xcut*stddev(realflux-mean(realflux))/stddev(xcut)

if binratio ne 0 then begin
  RESOLVE_ROUTINE, 'binner', /IS_FUNCTION
  simxrebin=binner(tcut,xcut,binsize=binwidth)
  if n_elements(simxrebin.xmean) ne Ncutrebin then begin
    print,'something strange with rebinning sim LC!'
    print,'simx rebin from binner: ',n_elements(simxrebin.xmean)
    print,'Ncutrebin from Ncut*binwidth: ',Ncutrebin
    ;return
  endif
  xrebin=simxrebin.ymean
  dxrebin=simxrebin.yerr
  trebin=simxrebin.xmean-binwidth/2.
endif else begin
  xrebin=xcut
  trebin=tcut
  dxrebin=sqrt(xrebin)
endelse 

;need to readjust stddev and mean of rebinned LC:
;xrebin=meanf+(xrebin-mean(xrebin))*(uf/stddev(xrebin))
;dxrebin=dxrebin*(uf/stddev(xrebin))

print,'std dev of unbinned and binned sims are: ',stddev(x0),stddev(xrebin)

;4. Pad real LC with zeros
RESOLVE_ROUTINE, 'Padding_func', /IS_FUNCTION
Vpadded=padding_func(realtime,realflux,binwidth)
VPadpsd=dblarr(Vpadded.padN/2)
VPadfreq=dblarr(Vpadded.padN/2)

;5. Apply window function to sim LC
RESOLVE_ROUTINE, 'Window_func', /IS_FUNCTION
Windowfuncbool=Window_func(realtime,binwidth,trebin,xrebin)

device,filename='simLCandPSD_beta'+string(beta,format='(f0.2)')+'_u'+strtrim(string(uf,format='(i)'),1)+'.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

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

realtime=realtime*86.400
tcut=tcut*86.400
trebin=trebin*86.400

;plot,realtime-realtime[0],(realflux-mean(realflux)),psym=1,color=Red
;plot,realtime-realtime[0],(realflux),ystyle=3,/xstyle,psym=1,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;,color=Red
plot,tcut,xcut,ystyle=3,/xstyle,psym=1,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;

;oplot,t,simx,color=Orange,psym=3
oplot,tcut,xcut,color=Orange,psym=1
oplot,trebin,xrebin,color=Blue,psym=4
;oplot,realtime-realtime[0],(realflux),color=Red,psym=1

;xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[max(realflux)*0.96,max(realflux)*0.96],'sim LC ',font=0,color=Blue
;xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[max(realflux)*0.91,max(realflux)*0.91],'input LC ',font=0,color=Red
xyouts,[max(tcut)*0.75,max(tcut)*0.75],[(max(xcut)-min(xcut))*0.9+min(xcut),min(xcut)+(max(xcut)-min(xcut))*0.9],'rebinned ',font=0,color=Blue
xyouts,[max(tcut)*0.75,max(tcut)*0.75],[(max(xcut)-min(xcut))*0.85+min(xcut),min(xcut)+(max(xcut)-min(xcut))*0.85],'fast sampled ',font=0,color=Orange

realtime=realtime/86.400D
tcut=tcut/86.400D
trebin=trebin/86.400D


;device, /close_file
;;end
;
;device,filename='simPSD.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8
;
;device, helvetica=1
;device, isolatin1=1
;;set_plot,'x'
;
;device, decomposed=1
;index=100;

!p.position=[0.15,0.08,0.95,0.48]

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freqFunc=dblarr(Ncutrebin/2)
psdFunc=dblarr(Ncutrebin/2)
freq2=dblarr(Ncutrebin/2)
psd2=dblarr(Ncutrebin/2)

simfreqFunc=dblarr(N/2)
simpsdFunc=dblarr(N/2)
fftfuncPSDbool=PSD_fft_func(trebin,xrebin,dxrebin,binwidth,freqFunc,psdFunc)
simfftfuncPSDbool=PSD_fft_func(t,simx,sqrt(simx),finebin*binwidth,simfreqFunc,simpsdFunc)

;print,simfreqFunc,simpsdFunc

fftPSDPadded=PSD_fft_func(Vpadded.padt,Vpadded.padx,sqrt(Vpadded.padx),binwidth,VPadfreq,VPadpsd)

cutfreq=dblarr(Ncut/2)
cutpsd=dblarr(Ncut/2)

cutPSDbool=PSD_fft_func(tcut,xcut,sqrt(xcut),finebin*binwidth,cutfreq,cutpsd)

; set binning factor
  binfactor = 10.
  mcutrebin = long(1.0* (Ncutrebin/2) / binfactor)  
  mcut = long(1.0D* (Ncut/2L) / binfactor)
  m = long(1.0D* (N/2L) / binfactor)
  mV = long(1.0* (Vpadded.padN/2) / binfactor)

  print,m

  nf = m*binfactor
  nfcut = mcut*binfactor
  nfcutrebin = mcutrebin*binfactor
  nfV = mV*binfactor

; remove underfilled bin from end
  pV = VPadpsd[0:nfV-1]
  fV = VPadfreq[0:nfV-1]
  ps = simpsdFunc[0:nf-1]
  fs = simfreqFunc[0:nf-1]
  pc = cutpsd[0:nfcut-1]
  fc = cutfreq[0:nfcut-1]
  fcr = freqFunc[0:nfcutrebin-1]
  pcr = psdFunc[0:nfcutrebin-1]

; rebin data to m bins
  psdV = rebin(pV,mV)
  frqV = rebin(fV,mV)
  psds = rebin(ps,m)
  frqs = rebin(fs,m)
  psdc = rebin(pc,mcut)
  frqc = rebin(fc,mcut)
  psdcr = rebin(pcr,mcutrebin)
  frqcr = rebin(fcr,mcutrebin)

; theoretical error (from chi-square distribution)
  errV = sqrt(1.0/binfactor)*psdV
  errs = sqrt(1.0/binfactor)*psds
  errc = sqrt(1.0/binfactor)*psdc
  errcr = sqrt(1.0/binfactor)*psdcr
  errfV = sqrt(1.0/binfactor)*frqV
  errfs = sqrt(1.0/binfactor)*frqs
  errfc = sqrt(1.0/binfactor)*frqc       
  errfcr = sqrt(1.0/binfactor)*frqcr


;plot,VPadfreq,VPadpsd,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(VPadfreq)*1e-1,max(VPadfreq)*1e2],yrange=[min(VPadpsd)*1e-2,max(VPadpsd)*1e2]
;;,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];
;oplot,simfreqFunc,simpsdFunc,color=Orange
;oplot,cutfreq,cutpsd,color=Green
;oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0],linestyle=2
;oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta,linestyle=1
;oplot,freqFunc,psdFunc,color=Blue
;oplot,VPadfreq,VPadpsd,color=Red

;plot,frqV,psdV,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(VPadfreq)*0.9,max(VPadfreq)*5],yrange=[min(VPadpsd)*4e-1,max(VPadpsd)*1e1],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'
plot,frqcr,psdcr,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(frqc)*0.1,max(frqc)*2],yrange=[min(psdcr)*5e-2,max(psdc)*2.],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'

;,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];
;oplot,frqs,psds,color=Orange,psym=8
;oplot,frqc,psdc,color=Green,psym=8
;oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0],linestyle=2
;oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta,linestyle=1
;oplot,frqcr,psdcr,color=Blue
;oplot,frqV,psdV,color=Red

;oploterror,frqs,psds,errs,color=Orange,errcolor=Orange,/nohat,psym=8
oploterror,frqc,psdc,errfc,errc,color=Orange,errcolor=Orange,/nohat,psym=8
oploterror,frqcr,psdcr,errfcr,errcr,color=Blue,errcolor=Blue,/nohat,psym=8
;oploterror,frqV,psdV,errfV,errV,color=Red,errcolor=Red,/nohat,psym=8
;oploterror,frqcr,psdcr,errcr,color=Blue,errcolor=Blue,/nohat,psym=4

;Simpsdrebin=binner(alog10(freqFunc),alog10(psdFunc),nperbin=binfactor)
;Vpsdrebin=binner(alog10(VPadfreq),alog10(VPadpsd),  nperbin=binfactor)

Simpsdrebin=binner((freqFunc),(psdFunc),nperbin=binfactor)
Vpsdrebin=binner((VPadfreq),(VPadpsd),nperbin=binfactor)

psdcrfit=LINFIT(alog10(freqFunc),alog10(psdFunc),CHISQ=psdcrChisq)
;psdVfit=LINFIT(alog10(frqV),alog10(psdV),CHISQ=psdVChisq)
psdVfit=LINFIT(alog10(VPadfreq),alog10(VPadpsd),CHISQ=psdVChisq)
psdcutfit=LINFIT(alog10(cutfreq),alog10(cutpsd),CHISQ=psdcutChisq)

;oploterror,10^Simpsdrebin.xmean,alog10(Simpsdrebin.ymean),alog10(Simpsdrebin.yerr),color=Cyan,errcolor=Cyan,/nohat,psym=8
;oploterror,10^Vpsdrebin.xmean,alog10(Vpsdrebin.ymean),alog10(Vpsdrebin.yerr),color=Orange,errcolor=Orange,/nohat,psym=8

;oploterror,Simpsdrebin.xmean,(Simpsdrebin.ymean),(Simpsdrebin.yerr),color=Cyan,errcolor=Cyan,/nohat,psym=8
;oploterror,Vpsdrebin.xmean,(Vpsdrebin.ymean),(Vpsdrebin.yerr),color=Orange,errcolor=Orange,/nohat,psym=8

oplot,freqFunc,10^((psdcrfit[0])+(psdcrfit[1])*alog10(freqFunc)),thick=5,linestyle=3,color=Blue
;oplot,VPadfreq,10^((psdVfit[0])+(psdVfit[1])*alog10(VPadfreq)),thick=5,linestyle=3,color=Red
oplot,cutfreq,10^((psdcutfit[0])+(psdcutfit[1])*alog10(cutfreq)),thick=5,linestyle=3,color=Orange

;xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdcr)*11.,max(psdcr)*11.],'rebinned:'+string(psdcrfit[1], format='(f0.3)'),font=0,color=Blue
;xyouts,[max(frqcr)*0.15,max(frqcr)*0.15],[max(psdcr)*3.,max(psdcr)*3.],'input LC index:'+string(psdVfit[1], format='(f0.3)'),font=0,color=Red
;xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdcr)*4.,max(psdcr)*4.],'fast sampled:'+string(psdcutfit[1], format='(f0.3)'),font=0,color=Orange

xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdcr)*11.,max(psdcr)*11.],'gapped:'+string(psdcrfit[1], format='(f0.3)'),font=0,color=Blue
;xyouts,[max(frqcr)*0.15,max(frqcr)*0.15],[max(psdcr)*3.,max(psdcr)*3.],'input LC index:'+string(psdVfit[1], format='(f0.3)'),font=0,color=Red
xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdcr)*4.,max(psdcr)*4.],'sim:'+string(psdcutfit[1], format='(f0.3)'),font=0,color=Orange

;xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdc)*1.,max(psdc)*1.],'rebinned:'+string(psdcrfit[1], format='(f0.3)'),font=0,color=Blue
;xyouts,[max(frqcr)*0.1,max(frqcr)*0.1],[max(psdc)*0.8,max(psdc)*0.8],'fast sampled:'+string(psdcutfit[1], format='(f0.3)'),font=0,color=Orange


;device, /close_file

;data=dblarr(4,N);
;data=dblarr(2,N)
;data[0,*]=t
;data[1,*]=x
;data[2,*]=w
;data[3,*]=psd
;openw,lun,'simLC_VWindow.dat',/get_lun
;printf,lun,data
;Free_lun,lun

;psdFunc=psdFunc/(mean(xrebin^2))
;VPadpsd=VPadpsd/(mean(Vpadded.padx^2))
outstruct = {                  $
       srN: Ncutrebin,            $
       srt: trebin,               $
       srx: xrebin,               $
       srFreq: freqFunc,          $
       srPSD: psdFunc,            $
       realFreq: VPadfreq,         $
       realPSD: VPadpsd }
 
return, outstruct


end

