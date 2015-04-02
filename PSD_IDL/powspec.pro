;pro powspec
FUNCTION powspec,t,x

; -------------------------------------------------
; IDL procedure to compute the normalised power
; spectral density (PSD) of a time series using
; the evenly binned periodogram (|DFT|^2).
; -------------------------------------------------

; -------------------- LOAD DATA --------------

; data = read_table()
;  t = data[0,*]
;  x = data[1,*]
  n = n_elements(x)
  s=size(x)
  if (s[0] gt 1) then x=reform(x,n,/overwrite)
  s=size(t)
  if (s[0] gt 1) then t=reform(t,n,/overwrite)
 
; ----------- COMPUTE PERIODOGRAM ------------

; subtract mean value

  mean_x = mean(x)
  x = x - mean_x

; time sampling rate

  dT = t[1]-t[0]    
  df = 1.0/(dT*n)

; convert time resolution to exact 1/2^J

;  j = round(alog(1.0/dT)/alog(2.0))
;  dT = double(1.0/2^j)
;  df = double(2^j)/n

  print,'-- Time resolution: ',dT
  print,'-- Frequency resolution: ',df

; calculate DFT

  nf = n/2                      ; no. +ve frequencies
  f = (findgen(nf)+1) * df      ; the +ve Fourier frequencies
  dft = fft(x, 1, dimension=1)  ; (NB: 1/N factor when for forward DFT)

; extract only positive (non-zero) frequency components
; and square to get power

  per = abs(dft[1:nf])^2

; apply rms/mean normalisation (NB: fft applied 1/N factor already)

  per = per * (2.0 * dT / mean_x^2 / float(n))

; ------------ BIN DATA -----------

; set binning factor

  binfactor = 20
  m = long(1.0* nf / binfactor)  
  nf = m*binfactor

; remove underfilled bin from end

  per = per[0:nf-1]
  f = f[0:nf-1]

; rebin data to m bins

  psd = rebin(per,m)
  frq = rebin(f,m)

; theoretical error (from chi-square distribution)

  err = sqrt(1.0/binfactor)*psd

; subtract Poisson noise

;  psd = psd - 2.0/mean_x

; -------------- PLOT DATA ---------------

;  window,0,retain=2

;  plot, frq, psd, /xlog, /ylog, psym=10, min_value=0, $
;    xrange=[min(frq),max(frq)],yrange=[5e-5,0.2],/xstyle, $
;    xthick=3, ythick=3, charthick=3, $
;    position=[0.15,0.15,0.9,0.9], charsize=1.4, font=-6, $
;    xtitle="!6Frequency",ytitle="Power"

;  err_plot,frq,psd-err,psd+err,width=0, min_value=0

outstruct = {                  $
       freq:  frq,         $
       psd:   psd,         $
       dpsd:  err }
 
return, outstruct

end
