;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs, put it in a file
 
FUNCTION PSDsim_histo_plot_func, Nsim, binwidth, binratio,beta

tspan=0L
;pseudo input of Nsim data points
pseudotime=findgen(Nsim)*binwidth
pseudoflux=RANDOMU(seed1,Nsim,/normal)

tspan=(pseudotime[Nsim-1]-pseudotime[0])/binwidth

;3. Simulate a longer LC with finer sampling rate
;and with underlying power-law psd of index -beta
RESOLVE_ROUTINE, 'simLongLC_func', /IS_FUNCTION
simLong=simLongLC_func(pseudotime,pseudoflux, beta, binwidth,binratio)

N=simLong.simN ;; N: total number of evenly-spaced points for sim time series, 
               ;; it is 2^int, longer than real x(t) with finer sampling
scatter=1.0
finebin=1.0D/(2L^binratio) ;;  time interval (bin width)
x0=simLong.simx
t=simLong.simt
meanf=mean(pseudoflux)
;uf=sigma(realflux)
uf=stddev(pseudoflux)

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

print,'N ',N,'Ncutrebin ',Ncutrebin
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
Vpadded=padding_func(pseudotime,pseudoflux,binwidth)
VPadpsd=dblarr(Vpadded.padN/2)
VPadfreq=dblarr(Vpadded.padN/2)

;5. Apply window function to sim LC
RESOLVE_ROUTINE, 'Window_func', /IS_FUNCTION
Windowfuncbool=Window_func(pseudotime,binwidth,trebin,xrebin)

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


psdcrfit=LINFIT(alog10(frqcr),alog10(psdcr),CHISQ=psdcrChisq)
psdVfit=LINFIT(alog10(frqV),alog10(psdV),CHISQ=psdVChisq)


outstruct = {                     $
       srN: Ncutrebin,            $
       srt: trebin,               $
       srx: xrebin,               $
       srFreq: freqFunc,          $
       srPSD: psdFunc,            $
       srIndex: psdcrfit[1],      $
       realFreq: VPadfreq,        $
       realPSD: VPadpsd }
 
return, outstruct


end

