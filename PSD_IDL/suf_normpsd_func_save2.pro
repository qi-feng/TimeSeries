;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real

FUNCTION suf_normpsd_func, realtime, realflux, binwidth, binratio, beta, Nsim

tspan=0L
NV=n_elements(realtime)
tspan=(realtime[NV-1]-realtime[0])/binwidth

;specify a filename base
;LCfname='LWA_LC_beta_'+strcompress(beta,/REMOVE_ALL)+'sim_'
;PSDfname='LWA_PSD_beta_'+strcompress(beta,/REMOVE_ALL)+'sim_'

RESOLVE_ROUTINE, 'LWA_func', /IS_FUNCTION

simChisq=dblarr(Nsim)
obsChisq=0.D
obsPSD=dblarr(Nsim)
 LWAstruct=LWA_func( realtime, realflux, binwidth, binratio,beta )
;;;;  naming convention: PSDsim is simulation, realPSD is from observation
 LCsim=dblarr(2,Nsim,LWAstruct.srN);
 PSDsim=dblarr(2,Nsim,LWAstruct.srN/2)
 realPSD=dblarr(2,LWAstruct.srN/2)
 realPSD[0,*]=LWAstruct.realFreq      
 realPSD[1,*]=LWAstruct.realPSD
 simPSDavg=dblarr(LWAstruct.srN/2)
 simPSDstdev=dblarr(LWAstruct.srN/2)

simnorm=mean(LWAstruct.srPSD)/mean(LWAstruct.realPSD)
simfit0=dblarr(Nsim)
simfit1=dblarr(Nsim)
obsfit=LINFIT(alog10(LWAstruct.realFreq),alog10(LWAstruct.realPSD),CHISQ=obsfitChisq)

for j=0,Nsim-1,1 do begin
    LWAstruct=LWA_func( realtime, realflux, binwidth, binratio,beta )
    LCsim[0,j,*]=LWAstruct.srt
    LCsim[1,j,*]=LWAstruct.srx
    PSDsim[0,j,*]=LWAstruct.srFreq
    PSDsim[1,j,*]=LWAstruct.srPSD
    simfit=LINFIT(alog10(LWAstruct.srFreq),alog10(LWAstruct.srPSD));,CHISQ=simfitChisq[ii])
    simfit0[j]=simfit[0]
    simfit1[j]=simfit[1]
    ;PSDsim[1,j,*]=LWAstruct.srPSD*mean(realPSD[1,*])/mean(LWAstruct.srPSD)
endfor
;note that the sim PSD needs to be renormalized to minimize obsChisq
;so here figure out the norm factor and apply it
meanfreq=mean(LWAstruct.srFreq)
simnorm=(10^((obsfit[0])+(obsfit[1])*alog10(meanfreq)))/(10^((mean(simfit0))+(mean(simfit1))*alog10(meanfreq)))
PSDsim[1,*,*]=PSDsim[1,*,*]*simnorm
print,'simnorm: ',simnorm

;RESOLVE_ROUTINE, 'binner', /IS_FUNCTION
;realPSD_log=alog10(LWAstruct.realFreq)

;realPSDrebin=binner(realPSD_log,LWAstruct.realPSD,binsize=0.1)

; set binning factor
  binfactor = 20
  m = long(1.0* n_elements(LWAstruct.realFreq) / binfactor)  
  nf = m*binfactor
; remove underfilled bin from end
  per = LWAstruct.realPSD[0:nf-1]
  f = LWAstruct.realFreq[0:nf-1]
; rebin data to m bins
  rebinpsd = rebin(per,m)
  rebinfreq = rebin(f,m)
  drebinpsd = sqrt(1.0/binfactor)*rebinpsd

;now iterate through each frequency i to calculate the goodnes of the fit
;note that the sim PSD needs to be renormalized to minimize obsChisq
 for i=0,LWAstruct.srN/2-1,1 do begin
   simPSDavg[i]=mean(PSDsim[1,*,i])
   simPSDstdev[i]=stddev(PSDsim[1,*,i])
   obsChisq+=(realPSD[1,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
   ;print,'adding up obsChisq for freq ',i,' : ',obsChisq
   for ii=0,Nsim-1,1 do begin
     simChisq[ii]+=(PSDsim[1,ii,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
     ;print,'adding up simChisq for sim ',ii,' and freq ',i, ' :  ',simChisq[ii]
   endfor
 endfor

m=0

for i=0,Nsim-1,1 do begin
  if obsChisq lt simChisq[i] then m++;
  print,obsChisq, simChisq[i], m
endfor

print,'success fraction is: ',double(m)/double(Nsim)
;return, double(m)/double(Nsim)

;obsChisq=0.

outstruct = {                  $
       suf: double(m)/double(Nsim),       $
       realFreqrebin: rebinfreq,  $
       realPSDrebin:  rebinpsd,      $
       drealPSDrebin: drebinpsd,       $
       realFreq:  realPSD[0,*],         $
       realPSD:   realPSD[1,*] }
 
return, outstruct

end

