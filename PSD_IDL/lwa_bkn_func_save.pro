;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs, put it in a file
 
FUNCTION lwa_bkn_func, realtime, realflux, binwidth, binratio,beta1,beta2,Fbkn

tspan=0L
NV=n_elements(realtime)
tspan=(realtime[NV-1]-realtime[0])/binwidth

;3. Simulate a longer LC with finer sampling rate
;and with underlying power-law psd of index -beta
RESOLVE_ROUTINE, 'simbknlc_func', /IS_FUNCTION
simLong=simbknlc_func(realtime,realflux, beta1,beta2,Fbkn, binwidth,binratio)

N=simLong.simN ;; N: total number of evenly-spaced points for sim time series, 
               ;; it is 2^int, longer than real x(t) with finer sampling
scatter=1.0
finebin=1.0D/(2L^binratio) ;;  time interval (bin width)
x0=simLong.simx
t=simLong.simt

;duplicate this long fine sim LC for contamination... sad...
x0=x0*stddev(realflux-mean(realflux))/stddev(x0)
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
xcut=xcut-mean(xcut)
xcut=xcut*stddev(realflux-mean(realflux))/stddev(xcut)

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
tr2=trebin
xr2=xrebin
dxr2=dxrebin

;4. Pad real LC with zeros
RESOLVE_ROUTINE, 'Padding_func', /IS_FUNCTION
Vpadded=padding_func(realtime,realflux,binwidth)
VPadpsd=dblarr(Vpadded.padN/2)
VPadfreq=dblarr(Vpadded.padN/2)

;5. Apply window function to sim LC
RESOLVE_ROUTINE, 'Window_func', /IS_FUNCTION
Windowfuncbool=Window_func(realtime,binwidth,trebin,xrebin)

;plot,realtime-realtime[0],(realflux-mean(realflux)),psym=1;,color=Red
;oplot,t,simx,color=Orange,psym=3
;oplot,tcut,xcut,color=Green,psym=1
;oplot,tr2,xr2,color=Cyan,psym=7,symsize=0.4
;oplot,trebin,xrebin,color=Blue,psym=4
;oplot,realtime-realtime[0],(realflux-mean(realflux)),color=Red,psym=1

;print,Ncutrebin, Vpadded.padN, 'cross check Nrebin and Padded LC N'
;print,N,Ncut

;device, /close_file
;end

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freqFunc=dblarr(Ncutrebin/2)
psdFunc=dblarr(Ncutrebin/2)
freq2=dblarr(Ncutrebin/2)
psd2=dblarr(Ncutrebin/2)

simfreqFunc=dblarr(N/2)
simpsdFunc=dblarr(N/2)
fftfuncPSDbool=PSD_fft_func(trebin,xrebin,dxrebin,binwidth,freqFunc,psdFunc)
;simfftfuncPSDbool=PSD_fft_func(t,simx,sqrt(simx),finebin*binwidth,simfreqFunc,simpsdFunc)

fftPSDPadded=PSD_fft_func(Vpadded.padt,Vpadded.padx,sqrt(Vpadded.padx),binwidth,VPadfreq,VPadpsd)

cutfreq=dblarr(Ncut/2)
cutpsd=dblarr(Ncut/2)

;cutPSDbool=PSD_fft_func(tcut,xcut,sqrt(xcut),finebin*binwidth,cutfreq,cutpsd)
;PSD2bool=PSD_fft_func(tr2,xr2,dxr2,binwidth,freq2,psd2)

;plot,VPadfreq,VPadpsd/(mean(Vpadded.padx^2)),/xstyle,/xlog,/ystyle,/ylog,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];

;oplot,simfreqFunc,simpsdFunc/(mean(simx^2)),color=Orange
;oplot,cutfreq,cutpsd/(mean(xcut^2)),color=Green
;oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0]/(mean(xrebin^2)),linestyle=2
;oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta/(mean(xrebin^2)),linestyle=1
;oplot,freq2,psd2/(mean(xr2^2)),color=Cyan
;oplot,freqFunc,psdFunc/(mean(xrebin^2)),color=Blue
;oplot,VPadfreq,VPadpsd/(mean(Vpadded.padx^2)),color=Red

;data=dblarr(4,N);
;data=dblarr(2,N)
;data[0,*]=t
;data[1,*]=x
;data[2,*]=w
;data[3,*]=psd
;openw,lun,'simLC_VWindow.dat',/get_lun
;printf,lun,data
;Free_lun,lun

psdFunc=psdFunc/(mean(xrebin^2))
VPadpsd=VPadpsd/(mean(Vpadded.padx^2))
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

