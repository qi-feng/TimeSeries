; This idl function deconvolves a window function from an padded observed LC
; to get a even sampled light curve;
; QF 2014-01-23
; t is a array containing each sampling time from a padded real light curve
; it is the output of Padding_func and the s
; x is the padded light curve points at t
; bin is the bin width of real LC in the *SAME unit* as t
; deconx is the deconvolved LC, evenly sampled at decont, 
; for a longer period than real LC, 

FUNCTION deconvolve_window_func, t, x, bin, decont, deconx

Nt=n_elements(t)
decont=dblarr(Nt)
deconx=dblarr(Nt)
decont=t

tspan=0L
tspan=(t[Nt-1]-t[0])/bin
;print,'observed time span in the unit of number of bins: ',tspan
  
;time since the first sample t[0]
t_0=dblarr(Nt)
t_0=t-t[0]

;print,'observed sampling times: ',t_0
window=dblarr(Nt)

for j=0L,Nt-1,1 do begin
  if ( x[j] ne 0. ) then begin 
       window[j]=1
  endif else begin
       window[j]=0
  endelse
endfor

fwindow=fft(window,-1)
fx=fft(x,-1)
deconx=fft((fx/fwindow),/inverse)

  outstruct = {                   $
        twindow: window,            $
        fwindow: fwindow,           $
        decont: decont,           $
        deconx: deconx  }
        
    return, outstruct


plot,decont,window,psym=8
oplot,decont,deconx,psym=8

;print,'freq: ',freq
;print,'psd, ',psd

END

