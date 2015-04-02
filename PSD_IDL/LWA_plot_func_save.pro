;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs, put it in a file
 
FUNCTION LWA_plot_func, realtime, realflux, binwidth, binratio,beta

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

tcut=dblarr(Ncut)
xcut=dblarr(Ncut)
tcut=t[0:Ncut-1]
xcut=x0[N-Ncut:N-1]
;xcut=xcut-mean(xcut)
;xcut=xcut*stddev(realflux-mean(realflux))/stddev(xcut)

RESOLVE_ROUTINE, 'binner', /IS_FUNCTION
simxrebin=binner(tcut,xcut,binsize=binwidth)
if n_elements(simxrebin.xmean) ne Ncutrebin then begin
  print,'something strange with rebinning sim LC!'
  print,'simx rebin from binner: ',n_elements(simxrebin)
  print,'Ncutrebin from Ncut*binwidth: ',Ncutrebin
  ;return
endif
xrebin=simxrebin.ymean
dxrebin=simxrebin.yerr
trebin=simxrebin.xmean-binwidth/2.
;need to readjust stddev and mean of rebinned LC:
xrebin=meanf+(xrebin-mean(xrebin))*(uf/stddev(xrebin))
dxrebin=dxrebin*(uf/stddev(xrebin))

print,'std dev of unbinned and binned sims are: ',stddev(x0),stddev(xrebin)

;4. Pad real LC with zeros
RESOLVE_ROUTINE, 'Padding_func', /IS_FUNCTION
Vpadded=padding_func(realtime,realflux,binwidth)
VPadpsd=dblarr(Vpadded.padN/2)
VPadfreq=dblarr(Vpadded.padN/2)

;5. Apply window function to sim LC
RESOLVE_ROUTINE, 'Window_func', /IS_FUNCTION
Windowfuncbool=Window_func(realtime,binwidth,trebin,xrebin)

device,filename='simLCandPSD.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

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

;plot,realtime-realtime[0],(realflux-mean(realflux)),psym=1,color=Red
plot,realtime-realtime[0],(realflux),ystyle=3,psym=1,xtitle='!17Time',ytitle='!17Input flux unit ',font=0;,color=Red
;oplot,t,simx,color=Orange,psym=3
;oplot,tcut,xcut,color=Green,psym=1
oplot,trebin,xrebin,color=Blue,psym=4
oplot,realtime-realtime[0],(realflux),color=Red,psym=1

;xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[max(realflux)*0.96,max(realflux)*0.96],'sim LC ',font=0,color=Blue
;xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[max(realflux)*0.91,max(realflux)*0.91],'input LC ',font=0,color=Red
xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[(max(realflux)-min(realflux))*0.95+min(realflux),min(realflux)+(max(realflux)-min(realflux))*0.95],'sim LC ',font=0,color=Blue
xyouts,[max(realtime-realtime[0])*0.8,max(realtime-realtime[0])*0.8],[(max(realflux)-min(realflux))*0.9+min(realflux),min(realflux)+(max(realflux)-min(realflux))*0.9],'input LC ',font=0,color=Red

;
;print,'sigma of sim lc: ',sigma(simx)
;print,'sigma of real lc: ',sigma(realflux)
;print,Ncutrebin, Vpadded.padN, 'cross check Nrebin and Padded LC N'
;print,N,Ncut

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
  binfactor = 20.
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

;plot,VPadfreq,VPadpsd,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(VPadfreq)*1e-1,max(VPadfreq)*1e2],yrange=[min(VPadpsd)*1e-2,max(VPadpsd)*1e2]
;;,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];
;oplot,simfreqFunc,simpsdFunc,color=Orange
;oplot,cutfreq,cutpsd,color=Green
;oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0],linestyle=2
;oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta,linestyle=1
;oplot,freqFunc,psdFunc,color=Blue
;oplot,VPadfreq,VPadpsd,color=Red

plot,frqV,psdV,psym=8,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(VPadfreq)*0.9,max(VPadfreq)*5],yrange=[min(VPadpsd)*2e-1,max(VPadpsd)*1e1],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0
;,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];
;oplot,frqs,psds,color=Orange,psym=8
;oplot,frqc,psdc,color=Green,psym=8
;oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0],linestyle=2
;oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta,linestyle=1
;oplot,frqcr,psdcr,color=Blue
;oplot,frqV,psdV,color=Red

;oploterror,frqs,psds,errs,color=Orange,errcolor=Orange,/nohat,psym=8
;oploterror,frqc,psdc,errc,color=Green,errcolor=Green,/nohat,psym=8
oploterror,frqcr,psdcr,errcr,color=Blue,errcolor=Blue,/nohat,psym=8
oploterror,frqV,psdV,errV,color=Red,errcolor=Red,/nohat,psym=8

psdcrfit=LINFIT(alog10(frqcr),alog10(psdcr),CHISQ=psdcrChisq)
psdVfit=LINFIT(alog10(frqV),alog10(psdV),CHISQ=psdVChisq)

oplot,frqcr,10^((psdcrfit[0])+(psdcrfit[1])*alog10(frqcr)),thick=5,linestyle=3,color=Blue
oplot,frqV,10^((psdVfit[0])+(psdVfit[1])*alog10(frqV)),thick=5,linestyle=3,color=Red

xyouts,[max(frqcr)*0.15,max(frqcr)*0.15],[max(psdcr)*9.,max(psdcr)*9.],'sim LC index:'+string(psdcrfit[1], format='(f0.3)'),font=0,color=Blue
xyouts,[max(frqcr)*0.15,max(frqcr)*0.15],[max(psdcr)*3.,max(psdcr)*3.],'input LC index:'+string(psdVfit[1], format='(f0.3)'),font=0,color=Red

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

