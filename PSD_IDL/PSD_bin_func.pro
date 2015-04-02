;This is an idl function to calculate power spectral density 
;from a discretely unevenly sampled light curve;
;following Uttley 2002
;QF 2014-01-01

;mjd,flux,dflux are input light curve data points
;freq, psd are output frequencies and normalized power spectrum, 
;of a DIFFERENT dimension as light curve

FUNCTION PSD_bin_func, mjd, flux, dflux,Tbin, freq, psd

onemjd=86400.

m=mean(flux)
N=n_elements(flux)
;Here subtract the mean, the res has a mean of zero
res=flux-m

Tdur=(mjd[N-1]-mjd[0])*onemjd
freq_min=1/Tdur
;freq_Nyq=N/(2*Tdur)
freq_Nyq=1/(2*Tbin)

Nfreq=long(Tdur/(2*Tbin))

print,'N: ',N
print,'N_from_bin: ',Nfreq

pi=3.14159
;freq=dblarr(N)
freq=dblarr(Nfreq)
;ft_real=dblarr(N)
;ft_imag=dblarr(N)
;psd=dblarr(N)
ft_real=dblarr(Nfreq)
ft_imag=dblarr(Nfreq)
psd=dblarr(Nfreq)

;for j=0,N-1,1 do begin
for j=0L,Nfreq-1,1 do begin
  freq[j]=freq_min*(j+1)
  for jj=0,N-1,1 do begin
    ft_real[j]+=cos(2*pi*freq[j]*mjd[jj]*onemjd)*res[jj]
    ft_imag[j]+=sin(2*pi*freq[j]*mjd[jj]*onemjd)*res[jj]
  endfor
  ;below is PSD: the modulus squared of LC's DFT
  psd[j]=ft_real[j]*ft_real[j]+ft_imag[j]*ft_imag[j]
  ;normalize PSD
  ;psd[j]=2*Tdur/(m*m*N*N)*psd[j]
  psd[j]=2*Tdur/(N*N)*psd[j]

endfor

;print,'freq: ',freq
;print,'psd, ',psd

print,'size of freq in procedure: ',size(freq)
END

