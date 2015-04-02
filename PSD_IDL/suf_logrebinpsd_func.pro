;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real

FUNCTION suf_logrebinpsd_func, realtime, realflux, binwidth, binratio, beta, Nsim

tspan=0L
NV=n_elements(realtime)
tspan=(realtime[NV-1]-realtime[0])/binwidth

;specify a filename base
LCfname='LWA_LC_beta_'+strcompress(beta,/REMOVE_ALL)+'sim_'
PSDfname='LWA_PSD_beta_'+strcompress(beta,/REMOVE_ALL)+'sim_'

RESOLVE_ROUTINE, 'LWA_logrebin_func', /IS_FUNCTION

simChisq=dblarr(Nsim)
obsChisq=0.D
obsPSD=dblarr(Nsim)
 LWAstruct=LWA_logrebin_func( realtime, realflux, binwidth, binratio,beta )
;;;;  naming convention: PSDsim is simulation, realPSD is from observation
 LCsim=dblarr(2,Nsim,LWAstruct.srN);
 PSDsim=dblarr(2,Nsim,n_elements(LWAstruct.srFreq))
 ;print,'1st LWA func: ',n_elements(LWAstruct.srFreq)
 ;realPSD=dblarr(2,LWAstruct.srN/2)
 RealPSD=dblarr(2,n_elements(LWAstruct.realFreq))
 RealPSD[0,*]=LWAstruct.realFreq      
 RealPSD[1,*]=LWAstruct.realPSD
 simPSDavg=dblarr(n_elements(LWAstruct.srFreq))
 simPSDstdev=dblarr(n_elements(LWAstruct.srFreq))
for j=0,Nsim-1,1 do begin
    LWAstruct=LWA_logrebin_func( realtime, realflux, binwidth, binratio,beta )
    ;print,j, 'th LWA func: ',n_elements(LWAstruct.srFreq)
    LCsim[0,j,*]=LWAstruct.srt
    LCsim[1,j,*]=LWAstruct.srx
    PSDsim[0,j,*]=LWAstruct.srFreq
    PSDsim[1,j,*]=LWAstruct.srPSD
    ;PSDsim[1,j,*]=LWAstruct.srPSD*mean(RealPSD[1,*])/mean(LWAstruct.srPSD)
endfor

;RESOLVE_ROUTINE, 'binner', /IS_FUNCTION
;realPSDrebin=dblarr(2,LWAstruct.srN/2)
;realPSD_log=alog10(LWAstruct.realFreq)
realPSD_log=(LWAstruct.realFreq)
;realPSDrebin[1,*]=LWAstruct.realPSD
realPSDrebin=LWAstruct.realPSD
;realPSDrebin=binner(realPSD_log,LWAstruct.realPSD,binsize=0.1)
;oploterror, 10^(Vrebin.xmean), 10^(Vrebin.ymean), 10^Vrebin.ymean*((Vrebin.yerr/Vrebin.ymean)), psym=8, color=Orange,/nohat,ERRCOLOR=Orange


;now iterate through each frequency i for some calculation
;for i=0,LWAstruct.srN/2-1,1 do begin
for i=0,n_elements(LWAstruct.srFreq)-1,1 do begin
  simPSDavg[i]=mean(PSDsim[1,*,i])
  simPSDstdev[i]=stddev(PSDsim[1,*,i])
  obsChisq+=(RealPSD[1,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
  ;print,'adding up obsChisq for freq ',i,' : ',obsChisq
  for ii=0,Nsim-1,1 do begin
    simChisq[ii]+=(PSDsim[1,ii,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
   ; print,'adding up simChisq for sim ',ii,' and freq ',i, ' :  ',simChisq[ii]
  endfor
endfor

m=0

for i=0,Nsim-1,1 do begin
  if obsChisq lt simChisq[i] then m++;
  ;print,obsChisq, simChisq[i], m
endfor

print,'success fraction is: ',double(m)/double(Nsim)
;return, double(m)/double(Nsim)

outstruct = {                  $
       suf: double(m)/double(Nsim),       $
       realFreqrebin: 10^(LWAstruct.realFreq),  $
       realPSDrebin: (LWAstruct.realPSD),      $
       drealPSDrebin: LWAstruct.drealPSD,       $
       realFreq:  RealPSD[0,*],         $
       realPSD:   RealPSD[1,*] }
 
return, outstruct

end

