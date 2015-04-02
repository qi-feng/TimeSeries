;Simulate light curves with power law noise + a gauss flare, assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;w=2*pi*f
;also assume mean and rms

;1, Choose power-law index beta for PSD
beta1=1.0
beta2=1.5
Fbkn=0.007

set_plot,'x'

psize=0.6
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;

;2, For each frequency w_i, generate two Gaussian distributed random numbers re and im, multiply then by (1./w_i)^(beta/2.) and get the real and imaginary part of the Fourier transform of the data

N=2^12  ;; N: total number of evenly-spaced points for time series, 
       ;; must be an even number
rms=30.
m=300.0
;binwidth=1.0 ;;  time interval (bin width)
;binwidth=10.D/(24.*60.)
binwidth=50.D/(86400.D)

x=dblarr(N)
t=findgen(N)*binwidth

RESOLVE_ROUTINE, 'simshortbknlc_func', /IS_FUNCTION
funcSimLCbool=simshortbknlc_func(N,m,rms,beta1,beta2,Fbkn,binwidth,x)
;print,'X(t) should be real: ',real_part(x)
;plot,t,x;,psym=8;,yrange=[8.0,8.6]

;gauss_peak=dblarr(N)

gsig=10*binwidth
;gAmp=stddev(real_part(x))*10
gAmp=0
gauss_peak=gAmp*exp(-(findgen(N)-N/2)^2/(2*gsig^2))

x=real_part(x)+gauss_peak

;plot,x

;plot,gauss_peak

data=dblarr(2,N)
data[0,*]=t
data[1,*]=real_part(x)
;data[2,*]=w
;data[3,*]=psd
;openw,lun,'sim/simLC_beta'+strtrim(string(beta),1)+'+GaussPeak'+strtrim(string(gsig),1)+'.dat',/get_lun
;openw,lun,'sim/simLC_N512_beta'+strtrim(string(beta),1)+'.dat',/get_lun
openw,lun,'sim/simLCbkn_N'+strtrim(string(N,format='(i)'),1)+'_beta1_'+strtrim(string(beta1,format='(f0.2)'),1)+'_beta2_'+strtrim(string(beta2,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'_Fbkn'+strtrim(string(Fbkn,format='(E10.1)'),1)+'.dat',/get_lun

;openw,lun,'sim/simLCbkn_N'+strtrim(string(N,format='(i)'),1)+'_beta'+strtrim(string(beta,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'+GaussPeak_Width'+strtrim(string(gsig/binwidth,format='(i)'),1)+'.dat',/get_lun

;printf,lun,transpose(real_part(x))
printf,lun,data
Free_lun,lun

;device, /close_file

end

