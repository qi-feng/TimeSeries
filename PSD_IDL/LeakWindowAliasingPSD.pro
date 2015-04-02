;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real

set_plot, 'ps'
;1. Choose power-law index beta for PSD
beta=1.
device,filename='SimLWA2_PSDbeta1.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8
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

;2. Read a real LC, careful with unit in LC (assuming days)
readcol,'NightlyVERLC.txt',F='D,D,D',VLtime,Vflux,dVflux;

; have to specify bin width in the same unit as LC time
binwidth=1.
; duration of real LC
tspan=0L
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

;3. Simulate a longer LC with finer sampling rate
;and with underlying power-law psd of index -beta
RESOLVE_ROUTINE, 'simLongLC_func', /IS_FUNCTION
binratio=3; so that sim LC has smaller binwidth/2^3 
simLong=simLongLC_func(VLtime,Vflux, beta, binwidth,binratio)
;print,'N=',simLong.simN
;plot,simLong.simt,simLong.simx

;N=2L^14  
N=simLong.simN ;; N: total number of evenly-spaced points for sim time series, 
               ;; it is 2^int, longer than real x(t) with finer sampling
scatter=1.0
finebin=1.0D/(2L^binratio) ;;  time interval (bin width)
x0=simLong.simx
t=simLong.simt

;duplicate this long fine sim LC for contamination... sad...
x0=x0*stddev(Vflux-mean(Vflux))/stddev(x0)
simx=x0

;3. Cut a piece of the simulated LC of the same duration as real LC,
;   and then bin the cut sim LC to same bin width as observations
;   cutting from the tail of vector for best randomness
;Ncut=ulong(tspan/binwidth)
datapadN=2L^(uint(alog(tspan)/alog(2))+1)
Ncut=datapadN*(2L^binratio)
;Ncutrebin=2L^(uint(alog(tempNcutrebin)/alog(2))+1)
Ncutrebin=datapadN

tcut=dblarr(Ncut)
xcut=dblarr(Ncut)
tcut=t[0:Ncut-1]
;xcut=x0[0:Ncut-1]
;tcut=t[N-Ncut:N-1]
xcut=x0[N-Ncut:N-1]
xcut=xcut-mean(xcut)
xcut=xcut*stddev(Vflux-mean(Vflux))/stddev(xcut)
;;print,'Ncutrebin: ', Ncutrebin

RESOLVE_ROUTINE, 'binner', /IS_FUNCTION
;simxrebin=binner(t,x0,binsize=binwidth)
simxrebin=binner(tcut,xcut,binsize=binwidth)
;print,'rebin time lower edge: ',simxrebin.xmean-binwidth/2., ' rebin x: ',simxrebin.ymean,' rebin x error: ',simxrebin.yerr
;oploterror, 10^(simxrebin.xmean), 10^(simxrebin.ymean), 10^simxrebin.ymean*((simxrebin.yerr/simxrebin.ymean)), psym=8, color=Orange,/nohat,ERRCOLOR=Orange
if n_elements(simxrebin.xmean) ne Ncutrebin then begin
  print,'something strange with rebinning sim LC!'
  print,'simx rebin from binner: ',n_elements(simxrebin)
  print,'Ncutrebin from Ncut*binwidth: ',Ncutrebin
  ;return
endif
;xrebin=dblarr(rebinN)
;trebin=dblarr(rebinN)
xrebin=simxrebin.ymean
dxrebin=simxrebin.yerr
trebin=simxrebin.xmean-binwidth/2.
;print,'rebin time: ',trebin
;plot,trebin,xrebin
;oplot,VLtime-VLtime[0],(Vflux-mean(Vflux))*1e9,color=Red
tr2=trebin
xr2=xrebin
dxr2=dxrebin

;4. Pad real LC with zeros
RESOLVE_ROUTINE, 'Padding_func', /IS_FUNCTION
Vpadded=padding_func(VLtime,Vflux,binwidth)
VPadpsd=dblarr(Vpadded.padN/2)
VPadfreq=dblarr(Vpadded.padN/2)

;5. Apply window function to sim LC
RESOLVE_ROUTINE, 'Window_func', /IS_FUNCTION
Windowfuncbool=Window_func(VLtime,binwidth,trebin,xrebin)
;;;plot,t,x;,psym=8;,yrange=[8.0,8.6]

;plot,VLtime-VLtime[0],(Vflux-mean(Vflux)),psym=1;,color=Red
;oplot,t,simx,color=Orange,psym=3
;oplot,tcut,xcut,color=Green,psym=1
;;;;oplot,Vpadded.padt,Vpadded.padx,color=Red,psym=5
;oplot,tr2,xr2,color=Cyan,psym=7,symsize=0.4
;oplot,trebin,xrebin,color=Blue,psym=4
;oplot,VLtime-VLtime[0],(Vflux-mean(Vflux)),color=Red,psym=1

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
simfftfuncPSDbool=PSD_fft_func(t,simx,sqrt(simx),finebin*binwidth,simfreqFunc,simpsdFunc)

fftPSDPadded=PSD_fft_func(Vpadded.padt,Vpadded.padx,sqrt(Vpadded.padx),binwidth,VPadfreq,VPadpsd)

cutfreq=dblarr(Ncut/2)
cutpsd=dblarr(Ncut/2)

cutPSDbool=PSD_fft_func(tcut,xcut,sqrt(xcut),finebin*binwidth,cutfreq,cutpsd)
PSD2bool=PSD_fft_func(tr2,xr2,dxr2,binwidth,freq2,psd2)

;plot,freqFunc,psdFunc,/xstyle,/xlog,/ystyle,/ylog
;plot,simfreqFunc,simpsdFunc/(mean(simx^2)),/xstyle,/xlog,/ystyle,/ylog,title='Green-Cut sim; Cyan-Rebin sime with window, Blue-padded data',yrange=[1e0,1e9],xrange=[1e-9,1e-4];
;plot,simfreqFunc,simpsdFunc/(mean(simx^2)),/xstyle,/xlog,/ystyle,/ylog,yrange=[1e0,1e9],xrange=[1e-9,1e-4];
plot,VPadfreq,VPadpsd/(mean(Vpadded.padx^2)),/xstyle,/xlog,/ystyle,/ylog,yrange=[1e-4,1e9],xrange=[1e-9,1e-4];

oplot,simfreqFunc,simpsdFunc/(mean(simx^2)),color=Orange
oplot,cutfreq,cutpsd/(mean(xcut^2)),color=Green
oplot,freqFunc,1./freqFunc*psdFunc[0]*freqFunc[0]/(mean(xrebin^2)),linestyle=2
oplot,freqFunc,(1./freqFunc)^beta*psdFunc[0]*(freqFunc[0])^beta/(mean(xrebin^2)),linestyle=1
oplot,freq2,psd2/(mean(xr2^2)),color=Cyan
oplot,freqFunc,psdFunc/(mean(xrebin^2)),color=Blue
;oplot,VPadfreq,VPadpsd*Vpadded.padN/2,color=Pink
;print,freq2
;print,'psd2: ',psd2
oplot,VPadfreq,VPadpsd/(mean(Vpadded.padx^2)),color=Red
;print,simpsdFunc
;print,simfreqFunc
;data=dblarr(4,N);
;data=dblarr(2,N)
;data[0,*]=t
;data[1,*]=x
;data[2,*]=w
;data[3,*]=psd
;openw,lun,'simLC_VWindow.dat',/get_lun
;printf,lun,data
;Free_lun,lun

device, /close_file

end

