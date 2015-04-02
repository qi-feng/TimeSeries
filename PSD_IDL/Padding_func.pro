; This idl function subtract the mean from a LC and 
; add bins with the value of zero to the gaps in a LC series
; to get 2^N bins ready for fft;
; QF 2014-01-23
; t is a array containing each sampling time from real light curve
; x is the LC measurements at each t
; bin is the bin width of real LC in the *SAME unit* as t
; return values: padN, padt, padx: 
; padt is an array of evenly spaced sequential time bins, starting from 0 
; padt has 2^N elements, N is the smallest int satisfying 2^N > n_elements(t)
; padx is the padded LC, evenly sampled at padt, 
; padx=x when there is a measurement, otherwise padx=mean(x)

FUNCTION Padding_func, t, x, bin
COMPILE_OPT idl2, HIDDEN

Nt=n_elements(t)
;print,' original LC number of points: ',Nt
m=mean(x)

tspan=0L
tspan=(t[Nt-1]-t[0])/bin
;print,'observed time span in the unit of number of bins: ',tspan
padN=2L^(uint(alog(tspan)/alog(2))+1)
;print,' padded LC number of points: ',padN

;time since the first sample t[0]
t_0=dblarr(Nt)
t_0=t-t[0]

;print,'observed sampling times: ',t_0
;padt=arr(padN)
padt=findgen(padN)*bin
padx=dblarr(padN)

jt=0

for j=0L,padN-1,1 do begin
  if ( jt le Nt-1) then begin
     if ( t_0[jt] ge (j-0.5D)*bin ) and ( t_0[jt] lt (j+0.5D)*bin ) then begin 
       padx[j]=x[jt] ; original LC value
       ;print,'for padding bin ',j,' at real time ',t_0[jt], ' , padding time ',padt[j], ' is sampled'
       jt++;
     endif else begin 
       padx[j]=m ; here pad zero (mean)
       ;print,'for padding bin ',j,' at padding time ',padt[j], ' real time ',t_0[jt], ' is *NOT* sampled! Pad zero'
     endelse 
  endif else begin
     padx[j]=m
     ;print,'for padding bin ',j,' at padding time ',padt[j], ' after time ',t_0[jt-1], ' is NOT sampled! pad zero'
  endelse
endfor

    outstruct = {                 $
        padN: padN,               $
        padt: padt,               $
        padx: padx}

    return, outstruct


;plot,simt,window,psym=8
;plot,simt,simx,psym=8

;print,'freq: ',freq
;print,'psd, ',psd

END

