;Simulate light curves assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;w=2*pi*f 
;input parameters: Number of points to simulate (must be even), 
;index for power law noise (beta), and binwidth
; x is a dblarr[N] to store the results
; t will be t=(findgen(N))*binwidth
; m and u is the desired mean and rms of the sim LC, m is used for Leahy norm
FUNCTION simLC_func, N , m, u, beta, binwidth, x

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
re=(RANDOMU(seed1,N/2+2,/normal))*((1./w)^(beta/2.))
im=(RANDOMU(seed2,N/2+2,/normal))*((1./w)^(beta/2.))

;; below is the complex Fourier transform of the LC
;; from freq_min to freq_Nyq, assume f[0]=0, valid?;;
comp=dcomplexarr(N)
comp[0]=complex(0.,0.)
for j=1L,N/2,1 do begin
    ;f[j]=freq_min*(j)
    comp[j]=complex(re[j],im[j])
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
;inverse rms normalization P=2T/(m^2*N^2)|F|^2??
;so F=N*m*sqrt(P)/sqrt(2T)??
;x=m+real_part(fft(2*comp*N*abs(m)/sqrt(2*Tdur),/inverse))
;x=m+real_part(fft(2*comp*sqrt(2*Tdur)/(m*N),/inverse))

;x=m+u/(2.D)*real_part(fft(comp,/inverse))
x=real_part(fft(comp,/inverse))
;x=m+x*u/sigma(x)
;the following change x to mean m and stddev u
x=m+(x-mean(x))*(u/stddev(x))

;x=m+real_part(fft(comp,/inverse))

;print,'searchkey',mean(x)

;x=x/mean(x)
;x=real_part(fft(comp,/inverse))
;x=2*Tdur/(N*N)*real_part(fft(comp*2*Tdur/(N*N),/inverse))
;x=real_part(fft(comp*2*Tdur/(N*N),/inverse))
;x=1/N*real_part(fft(comp,/inverse))

;x=real_part(fft(comp*(N*N)/(2*Tdur),/inverse))

end

