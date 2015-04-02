;Simulate light curves with power law noise + a gauss flare, assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;w=2*pi*f
;also assume mean and rms

set_plot,'x'
;set_plot, 'ps'
;device,filename='SimLCtest.eps', /encapsulated, xsize=8 ,ysize=5,/inches

psize=0.6
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;

;1, Choose power-law index beta for PSD
beta=1.0

;2, For each frequency w_i, generate two Gaussian distributed random numbers re and im, multiply then by (1./w_i)^(beta/2.) and get the real and imaginary part of the Fourier transform of the data

N=2^10  ;; N: total number of evenly-spaced points for time series, 
       ;; must be an even number
rms=0.5
m=1.0
binwidth=1.0 ;;  time interval (bin width)
x=dblarr(N)
t=findgen(N)

RESOLVE_ROUTINE, 'simLC_func', /IS_FUNCTION
funcSimLCbool=simLC_func(N,m,rms,beta,binwidth,x)
;print,'X(t) should be real: ',real_part(x)
;plot,t,x;,psym=8;,yrange=[8.0,8.6]

;gauss_peak=dblarr(N)

gsig=1
;gAmp=stddev(real_part(x))*10
gAmp=0
gauss_peak=gAmp*exp(-(findgen(N)-N/2)^2/(2*gsig^2))

x=real_part(x)+gauss_peak

plot,x

;plot,gauss_peak

data=dblarr(2,N)
data[0,*]=t
data[1,*]=real_part(x)
;data[2,*]=w
;data[3,*]=psd
openw,lun,'sim/simLC_beta'+strtrim(string(beta),1)+'+GaussPeak'+strtrim(string(gsig),1)+'.dat',/get_lun
;printf,lun,transpose(real_part(x))
printf,lun,data
Free_lun,lun

;device, /close_file

end

