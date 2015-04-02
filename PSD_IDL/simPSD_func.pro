;Simulate light curves assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;w=2*pi*f 
;input parameters: Number of points to simulate (must be even), 
;index for power law noise (beta), and binwidth
; freq and psd is a dblarr[N/2] to store the PSD results
FUNCTION simPSD_func, N , beta, binwidth, freq, psd
;set_plot,'x'

; For each frequency w_i, generate two Gaussian distributed random numbers re and im, multiply then by (1./w_i)^(beta/2.) and get the real and imaginary part of the Fourier transform of the data

scatter=1.0
t=findgen(N)*binwidth ; 
Tdur=N*binwidth
freq_min=1./Tdur
freq_Nyq=N/(2.*Tdur)  ;; In frequency domain, only N/2 points
f=dblarr(N)  ;; Define f(-w)=f*(w) <here * means conjugate>, so f has length N again

;; frequency are arranged so that f[0] contains the zero frequency component
;; from f[1] to f[N-1]: 
;; freq_min, 2*freq_min, ... , (N/2-1)*freq_min, freq_Nyq <i.e. N/2*freq_min>, -(N/2-1)*freq_min, ..., -2*freq_min, -freq_min
;; 
;; below assign positive freq from freq_min to freq_Nyq, assume f[0]=0, valid?
;; store them  from f[1] to f[N/2], values are
;; freq_min, 2*freq_min, ... , (N/2-1)*freq_min, freq_Nyq <i.e. N/2*freq_min>
for j=0,N/2,1 do begin
    f[j]=freq_min*(j)
endfor
;; below assign negative freq from -freq_Nyq+freq_min <i.e. -(N/2-1)*freq_min> to -freq_min
;; store them  from f[N/2+1] to f[N-1]
;; values are: -(N/2-1)*freq_min, ..., -2*freq_min, -freq_min
for j=N/2+1,N-1,1 do begin
    f[j]=freq_min*(j-N)
endfor
w=2*!pi*f

;; wPos is the physical positive frequencies
wPos=dblarr(N/2)
for j=0,N/2-1,1 do begin
      wPos[j]=freq_min*(j+1)
endfor

;; x will store the time series, which is the inverse FFT of "comp" below
;; x is guaranteed to be real since comp(-f)=comp*(f)
;x=dblarr(N)
re=dblarr(N/2+2)
im=dblarr(N/2+2)

;; below generate the real and imaginary part of 
;; the Fourier transform of the LC, 
;; which is a white noise on top of a power-law spectrum
re=(RANDOMU(seed1,N/2+2,/normal))*((1./w)^(beta/2.))*scatter
im=(RANDOMU(seed2,N/2+2,/normal))*((1./w)^(beta/2.))*scatter

;; below is the complex Fourier transform of the LC
;; from freq_min to freq_Nyq, assume f[0]=0, valid?;;
comp=dcomplexarr(N)
comp[0]=complex(0.,0.)
for j=1,N/2,1 do begin
    ;f[j]=freq_min*(j)
    comp[j]=complex(re[j],im[j])
endfor
;; below calculate poer for negative freq 
;; from -freq_Nyq+freq_min <i.e. -(N/2-1)*freq_min> to -freq_min
;; so that comp(-f)=conj(comp(f))
;; store them  from comp[N/2+1] to comp[N-1]
for j=N/2+1,N-1,1 do begin
    ;f[j]=freq_min*(j-N)
    comp[j]=complex(re[N-j],-im[N-j])
endfor

;print,'complex Fourier transform',comp

;3, Now convert comp back to time domain,
;   and get time series x(t), 
;   through inverse FFT

comp_xt=dcomplexarr(N)
comp_xt=fft(comp,/inverse)
;x=real_part(fft(comp,/inverse))

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION

;freqFunc=dblarr(N/2)
;psdFunc=dblarr(N/2)
fftfuncPSDbool=PSD_fft_func(t,comp_xt,sqrt(comp_xt),1,freq,psd)

;plot,freqFunc,psdFunc,/xstyle,/xlog,/ystyle,/ylog
;oplot,w,1./w*psdFunc[0]*freqFunc[0]
;plot,t,x;,psym=8;,yrange=[8.0,8.6]

;oplot,freqFunc,psdFunc

;plot,w,psd,psym=8,/xstyle,/xlog,/ystyle,/ylog,yrange=[1e-22,1e-14],xrange=[1e-9,1e-5],xtitle='!17Frequency (Hz)',ytitle='!17Power',charsize=1.2,font=0

;data=dblarr(2,N)
;data[0,*]=t
;data[1,*]=x
;;data[2,*]=w
;;data[3,*]=psd
;openw,lun,'simLCtest.dat',/get_lun
;printf,lun,data
;Free_lun,lun

;device, /close_file

end

