;Simulate light curves assuming
;broken power law spectrum density: S(f)~(1/f)^beta1; when f<Fbkn
;                           S(f)~(1/f)^beta2; when f>Fbkn (Timmer & Konig 1995)
;w=2*pi*f 
;input parameters: Number of points to simulate (must be even), 
;index for power law noise (beta), and binwidth
; x is a dblarr[N] to store the results
; t will be t=(findgen(N))*binwidth
; m and u is the desired mean and rms of the sim LC

FUNCTION simshortbknlc_func, N , m, u, beta1,beta2,Fbkn, binwidth, x
COMPILE_OPT idl2, HIDDEN

; For each frequency w_i, generate two Gaussian distributed random numbers re and im, multiply then by (1./w_i)^(beta/2.) and get the real and imaginary part of the Fourier transform of the data

t=(findgen(N))*binwidth ; 
Tdur=N*binwidth
freq_min=1./Tdur
freq_Nyq=N/(2.*Tdur)  ;; In frequency domain, only N/2 points
f=dblarr(N)  ;; Define f(-w)=f*(w) <* is conjugate>, so f has length N again

;; frequency are arranged so that f[0] contains the zero frequency component
;; from f[1] to f[N-1]: 
;; freq_min, 2*freq_min, ... , (N/2-1)*freq_min, freq_Nyq <i.e. N/2*freq_min>, -(N/2-1)*freq_min, ..., -2*freq_min, -freq_min
;; 
;; below assign positive freq from freq_min to freq_Nyq, assume f[0]=0, valid?
;; store them  from f[1] to f[N/2], values are
;; freq_min, 2*freq_min, ... , (N/2-1)*freq_min, freq_Nyq <i.e. N/2*freq_min>
for j=0L,N/2,1 do begin
    f[j]=freq_min*(j)
endfor
;; below assign negative freq from -freq_Nyq+freq_min <i.e. -(N/2-1)*freq_min> to -freq_min
;; store them  from f[N/2+1] to f[N-1]
;; values are: -(N/2-1)*freq_min, ..., -2*freq_min, -freq_min
for j=N/2L+1,N-1,1 do begin
    f[j]=freq_min*(j-N)
endfor
w=2*!pi*f

;; wPos is the physical positive frequencies
wPos=dblarr(N/2)
for j=0L,N/2-1,1 do begin
      wPos[j]=2*!pi*freq_min*(j+1)
endfor

;; x will store the time series, which is the inverse FFT of "comp" below
;; x is guaranteed to be real since comp(-f)=comp*(f)
;x=dblarr(N)
re=dblarr(N/2+2)
im=dblarr(N/2+2)
re=(RANDOMU(seed1,N/2+2,/normal))
im=(RANDOMU(seed2,N/2+2,/normal))

;; below generate the real and imaginary part of 
;; the Fourier transform of the LC, 
;; which is a white noise on top of a power-law spectrum
for j=0L,N/2-1,1 do begin
  if wPos[j] lt 2*!pi*Fbkn then begin
    re[j]=re[j]*((1./wPos[j])^(beta1/2.))
    im[j]=im[j]*((1./wPos[j])^(beta1/2.))
  endif else begin
    re[j]=re[j]*((1./wPos[j])^(beta2/2.))
    im[j]=im[j]*((1./wPos[j])^(beta2/2.))
  endelse
endfor


;; below is the complex Fourier transform of the LC
;; from freq_min to freq_Nyq, assume f[0]=0, valid?;;
comp=dcomplexarr(N)
comp[0]=complex(0.,0.)
for j=1L,N/2,1 do begin
    ;f[j]=freq_min*(j)
    comp[j]=complex(re[j-1],im[j-1])
endfor
;; below calculate poer for negative freq 
;; from -freq_Nyq+freq_min <i.e. -(N/2-1)*freq_min> to -freq_min
;; so that comp(-f)=conj(comp(f))
;; store them  from comp[N/2+1] to comp[N-1]
for j=N/2L+1,N-1,1 do begin
    ;f[j]=freq_min*(j-N)
    comp[j]=complex(re[N-j],-im[N-j])
endfor

;3, Now convert comp back to time domain,
;   and get time series x(t), 
;   through inverse FFT

;x=fft(comp,/inverse)
;print,'searchkey',mean(x)

x=real_part(fft(comp,/inverse))   
;x=m+x*u/sigma(x)     
;the following change x to mean m and stddev u
x=m+(x-mean(x))*(u/stddev(x))   

;x=x/mean(x)
;x=2*Tdur/(N*N)*real_part(fft(comp*2*Tdur/(N*N),/inverse))
;x=real_part(fft(comp*2*Tdur/(N*N),/inverse))
;x=1/N*real_part(fft(comp,/inverse))

;x=real_part(fft(comp*(N*N)/(2*Tdur),/inverse))

end

