; This idl function convolves the FFT of a window function to a simulated LC
; to get a discrete uneven sampling fashion similar to real light curve;
; QF 2014-01-23
; t is a array containing each sampling time from real light curve
; bin is the bin width of real LC in the *SAME unit* as t
; simx is the simulated LC, evenly sampled at simt, 
; for a longer period than real LC, 
; it will be changed to after multiplication of the window function

FUNCTION Window_func, t, bin, simt, simx

Nt=n_elements(t)
Nsim=n_elements(simt)
m=mean(simx)

;print,' real LC number of points: ',Nt
;print,' sim LC number of points: ',Nsim

tspan=0L
tspan=(t[Nt-1]-t[0])/bin
;print,'observed time span in the unit of number of bins: ',tspan
if (Nsim lt tspan) then begin
  print,' Nsim is smaller than LC span, need more sims! '
  return,-1
endif
  
;time since the first sample t[0]
t_0=dblarr(Nt)
t_0=t-t[0]

;print,'observed sampling times: ',t_0
window=dblarr(Nsim)

jt=0

for j=0L,Nsim-1,1 do begin
  if ( jt lt Nt-1) then begin
     ;;if ( t_0[jt] ge j*bin ) and ( t_0[jt] lt (j+1)*bin ) then begin 
     if ( t_0[jt] ge (j-0.5D)*bin ) and ( t_0[jt] lt (j+0.5D)*bin ) then begin 

     ;;if ( (t_0[jt] - j*bin) ge -(bin*0.001D) ) and ( (t_0[jt] - (j+1)*bin) le -(0.0001D*bin) ) then begin

       window[j]=1
       ;print,'for window ',j,' at time ',t_0[jt], ' is sampled'
       jt++;
     endif else begin 
       window[j]=0
       simx[j]=m
       ;print,'for window ',j,' at time ',t_0[jt], ' is *NOT* sampled!'
     endelse 
  endif else begin
     window[j]=0
     simx[j]=m
     ;print,'for window ',j,' after time ',t_0[jt], 'NOT sampled!'
  endelse
endfor

fwindow=fft(window,-1)
    outstruct = {                   $
        twindow: window,               $
        fwindow: fwindow}

    return, outstruct


;plot,simt,window,psym=8
;plot,simt,simx,psym=8

;print,'freq: ',freq
;print,'psd, ',psd

END

