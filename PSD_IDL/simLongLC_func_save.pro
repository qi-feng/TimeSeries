;Simulate light curves assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;w=2*pi*f 
;input parameters: t,x is a real LC, 
;index for power law noise (beta), and binwidth
;simulated length is 2^6 to 2^7 (for boost=7) longer than
;the input LC, with smaller bin width: binwidth/2^binratio
; output will be a struct containing:
; return values: simN, simt, simx: 
; simt is an array of evenly spaced sequential time bins, starting from 0 
; simt has simN(=2^integer) elements
; simx is the sim LC, evenly sampled at simt, with zero mean

FUNCTION simLongLC_func,t,x, beta, binwidth,binratio

; For each frequency w_i, generate two Gaussian distributed random numbers re and im, multiply then by (1./w_i)^(beta/2.) and get the real and imaginary part of the Fourier transform of the data

scatter=1.0
Nt=long(n_elements(t))
tspan=0L
tspan=(t[Nt-1]-t[0])/binwidth
m=mean(x)
;u=sigma(x)
u=stddev(x)

boost=6; how much longer the simulated LC is, abount 2^boost times longer
;finebin=1.0/(2.^4) ;;  how much smaller the bins are
finebin=1.0D/(2.^binratio) ;;  how much smaller the bins are

simN=2L^(uint(alog(tspan)/alog(2))+boost)/finebin
print,'simN ',simN
simt=(findgen(simN))*binwidth*finebin
simx=dblarr(simN)

RESOLVE_ROUTINE, 'simLC_func', /IS_FUNCTION
;funcSimLCbool=simLC_func(simN,beta,binwidth*finebin,simx)
funcSimLCbool=simLC_func(simN,m,u,beta,binwidth*finebin,simx)

outstruct = {                 $
    simN: simN,               $
    simt: simt,               $
    simx: simx}

return, outstruct


end

