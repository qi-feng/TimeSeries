;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real

FUNCTION suf_brokenpsd_func, realtime, realflux, binwidth, binratio, beta1,beta2,Fbkn, Nsim

tspan=0L
NV=n_elements(realtime)
tspan=(realtime[NV-1]-realtime[0])/binwidth

;specify a filename base
LCfname='LWA_LC_bkn_'+strcompress(beta1,/REMOVE_ALL)+strcompress(beta2,/REMOVE_ALL)+strcompress(Fbkn,/REMOVE_ALL)+'sim_'
PSDfname='LWA_PSD_bkn_'+strcompress(beta1,/REMOVE_ALL)+strcompress(beta2,/REMOVE_ALL)+strcompress(Fbkn,/REMOVE_ALL)+'sim_'

RESOLVE_ROUTINE, 'lwa_bkn_func', /IS_FUNCTION

simChisq=dblarr(Nsim)
obsChisq=0.D
obsPSD=dblarr(Nsim)
 LWAstruct=lwa_bkn_func( realtime, realflux, binwidth, binratio,beta1,beta2,Fbkn )
;;;; some strange naming convention: PSDdata is simulation, realPSD is from observation
 LCdata=dblarr(2,Nsim,LWAstruct.srN);
 PSDdata=dblarr(2,Nsim,LWAstruct.srN/2)
 realPSD=dblarr(2,LWAstruct.srN/2)
 realPSD[0,*]=LWAstruct.realFreq      
 realPSD[1,*]=LWAstruct.realPSD
 simPSDavg=dblarr(LWAstruct.srN/2)
 simPSDstdev=dblarr(LWAstruct.srN/2)
for j=0,Nsim-1,1 do begin
    LWAstruct=lwa_bkn_func( realtime, realflux, binwidth, binratio,beta1,beta2,Fbkn )
    LCdata[0,j,*]=LWAstruct.srt
    LCdata[1,j,*]=LWAstruct.srx
    PSDdata[0,j,*]=LWAstruct.srFreq
    ;PSDdata[1,j,*]=LWAstruct.srPSD
    PSDdata[1,j,*]=LWAstruct.srPSD*mean(realPSD[1,*])/mean(LWAstruct.srPSD)
endfor

;now iterate through each frequency i for some calculation
for i=0,LWAstruct.srN/2-1,1 do begin
  simPSDavg[i]=mean(PSDdata[1,*,i])
  simPSDstdev[i]=stddev(PSDdata[1,*,i])
  obsChisq+=(realPSD[1,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
  ;print,'adding up obsChisq for freq ',i,' : ',obsChisq
  for ii=0,Nsim-1,1 do begin
    simChisq[ii]+=(PSDdata[1,ii,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
   ; print,'adding up simChisq for sim ',ii,' and freq ',i, ' :  ',simChisq[ii]
  endfor
endfor

m=0

for i=0,Nsim-1,1 do begin
  if obsChisq lt simChisq[i] then m++;
  ;print,obsChisq, simChisq[i], m
endfor

print,'success fraction is: ',double(m)/double(Nsim)
return, double(m)/double(Nsim)

end

