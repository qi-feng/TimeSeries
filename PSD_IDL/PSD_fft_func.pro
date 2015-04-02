;This is an idl function to calculate power spectral density 
;from a discretely unevenly sampled light curve with bin=1, t=findgen(N)+1;
;following Uttley 2002
;QF 2014-01-01

;t,fx,dx are padded input light curve data points, (2^N points)
;t should be in the unit of days
;freq, psd are output frequencies and normalized power spectrum, 
;of a DIFFERENT dimension as light curve with N/2 points
;binwidth in the unit of days

FUNCTION PSD_fft_func, t, x, dx,binwidth, freq, psd
COMPILE_OPT idl2, HIDDEN

binunit=86400.D
m=mean(x)
u=sigma(x)
N=n_elements(x)
;print,'xoxoxox', N, m
;Here subtract the mean, the res has a mean of zero
res=x-m

;Tdur in the unit of seconds
Tdur=(t[N-1]-t[0])*binunit
;freq_min in the unit of Hz
freq_min=1D/Tdur
;freq_Nyq=N/(2D*Tdur)
;Nyquist freq in the unit of Hz
freq_Nyq=1D/(2D*binwidth*binunit)

;print,'freq_min', freq_min
;print,'freq_Nyq', freq_Nyq

;make sure Nfreq is 2^int instead of 2^int-1
Nfreq=long(freq_Nyq/freq_min+0.9D)
;Nfreq=N/2
if Nfreq ne N/2L then begin
   print,'PSD_fft found number of points in LC strange '
   print,'number of points in LC: ',N
   print,'number of freq if assuming N_freq=N/2 ', N/2L
   print,'Number of freq points from freq_Nyq/freq_min: ',Nfreq
endif

;freq=dblarr(Nfreq)
ft=complexarr(N)
;psd=dblarr(Nfreq)

;; frequency are arranged so that ft[0] contains the zero frequency component
;; from ft[1] to ft[N-1]: 
;; freq_min, 2*freq_min, ... , (N/2-1)*freq_min, freq_Nyq <i.e. N/2*freq_min>, -(N/2-1)*freq_min, ..., -2*freq_min, -freq_min

;fft of zero mean time series
ft=fft(res,1, dimension=1) 
;print,'fft of padded LC: ',ft

;get sum of fft modulus squared for positive frequencies
;to get periodogram as estimate of PSD

psd = abs(ft[1:Nfreq])^2
;print,'raw psd',psd[1]
;normalize PSD
;psd=(2*binwidth*binunit)/(N)*psd
;psd=2*(Tdur/N)/N*psd
;psd=2.0D*Tdur/(N*N*m*m)*psd
psd=2.0D*binwidth*binunit/(N*m*m)*psd
;psd=2.0*Tdur/(N*m*m)*psd

;print,'norm psd',psd[1]
;print,'N,Tdur,norm:',N,Tdur,2*(Tdur/N)/N*psd
for j=0L,Nfreq-1,1 do begin
  freq[j]=freq_min*(j+1)
endfor

END

