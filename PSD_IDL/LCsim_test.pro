;args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;assume a .xxx file appendix
;arg_arr = STRMID(args,0,strlen(args)-4)

;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

;fname='simLC_N1024_beta0.50m300u60'
;fname='simLC_N1024_beta0.00m300u60'
;fname='simLC_N1024_beta-0.50m300u60'
;fname='simLC_N1024_beta1.50m300u60'
;fname='simLC_N1024_beta2.00m300u60'
;fname='simLC_N1024_beta1.00m300u60'
;fname='simLC_N512_beta1.00000'
;fname='simLC_N512_beta1_10m300u60'

COMPILE_OPT idl2, HIDDEN 

beta=1.5
fname='PSD_histo_of_simLC_beta'+string(beta,format='(f0.2)')

set_plot,'ps'
device,filename=fname+'.eps', /encapsulated, xsize=5 ,ysize=7,/inches,/color,bits_per_pixel=8

device, helvetica=1
device, isolatin1=1

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

; have to specify bin width in the same unit as LC time
binwidth=(50.D)/(86400.D)
;binwidth=1.
; duration of real LC
binratio=0; so that sim LC has smaller binwidth/2^3 

;pseudo input of Nlc data points
Nlc=512
rms=30.
fmean=300.
pseudotime=findgen(Nlc)*binwidth
pseudoflux=RANDOMU(seed1,Nlc,/normal)*rms+fmean

;now sim Nsim LCs
Nsim=500

PSDindex=dblarr(Nsim)
PSDnorm=dblarr(Nsim)
psdsimindex=dblarr(Nsim)
psdsimnorm=dblarr(Nsim)

RESOLVE_ROUTINE, 'PSDsim_histo_plot_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'simLongLC_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION

for j=0,Nsim-1,1 do begin
    LWAstruct=PSDsim_histo_plot_func(Nlc, binwidth, binratio,beta )
    ;PSDindex[j]=LWAstruct.srIndex
    psdsimfuncfit=LINFIT(alog10(LWAstruct.srFreq),alog10(LWAstruct.srPSD))
    PSDindex[j]=psdsimfuncfit[1]
    PSDnorm[j]=psdsimfuncfit[0]
    simLong=simLongLC_func(pseudotime,pseudoflux, beta, binwidth,binratio)
    xsim=simLong.simx
    tsim=simLong.simt
    tcut=tsim[0:(Nlc*2^binratio-1)]
    xcut=xsim[0:(Nlc*2^binratio-1)]
    print,tcut,xcut
    if binratio ne 0 then begin
      simxrebin=binner(tcut,xcut,binsize=binwidth)
      freqsimrebin=dblarr(long(n_elements(simxrebin.xmean)/2))
      psdsimrebin=dblarr(long(n_elements(simxrebin.xmean)/2))
      fftPSDsim=PSD_fft_func(simxrebin.xmean,simxrebin.ymean,simxrebin.yerr,binwidth,freqsimrebin,psdsimrebin)
    endif else begin
      freqsimrebin=dblarr(long(n_elements(xcut)/2))
      psdsimrebin=dblarr(long(n_elements(xcut)/2))
      fftPSDsim=PSD_fft_func(tcut,xcut,sqrt(xcut),binwidth,freqsimrebin,psdsimrebin)
    endelse

    freqsim=dblarr(long(n_elements(xcut)/2))
    psdsim=dblarr(long(n_elements(xcut)/2))
    ;freqsim=dblarr(long((simLong.simN)/2))
    ;psdsim=dblarr(long((simLong.simN)/2))
    ;fftPSDsim=PSD_fft_func(tsim,xsim,sqrt(xsim),binwidth/(2^binratio),freqsim,psdsim)
    ;fftPSDsim=PSD_fft_func(tcut,xcut,sqrt(xcut),binwidth/(2^binratio),freqsim,psdsim)

    ;print,freqsim,psdsim,xsim
    ;psdsimfit=LINFIT(alog10(freqsim),alog10(psdsim))
    psdsimfit=LINFIT(alog10(freqsimrebin),alog10(psdsimrebin))
    psdsimindex[j]=psdsimfit[1]
    psdsimnorm[j]=psdsimfit[0]
endfor

!p.multi=[0,1,2]
!p.position=[0.15,0.57,0.95,0.97] 

plothist,PSDindex,bin=0.01,font=0,charsize=1.2,xtitle='!17 Simulated PSD index',ytitle='!17 Histogram',color=Red;,yrange=[0,60];xrange=[-2,0],yrange=[0,35]

plothist,psdsimindex,bin=0.01,/overplot,color=Green

!p.position=[0.15,0.08,0.95,0.48]

plothist,PSDnorm,bin=0.1,font=0,charsize=1.2,xtitle='!17 Simulated PSD normalization',ytitle='!17 Histogram',color=Red
plothist,psdsimnorm,bin=0.1,/overplot,color=Green

device, /close_file

end

