;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real by calculating success fraction (SuF),
;following the PSRESP method described in Chatterjee 2008
COMPILE_OPT idl2, HIDDEN

;provide lc file name in command line, 
;e.g. idl -e ".run psd_suf.pro" -args LC.dat
;LC.dat has 3 columns: time, flux, flux error
;change the value of binwidth below correspondingly 

;args[0] is the data Hspec file e.g. simLC_N1024_beta1.00m300u60_noise0.01_Hspec_linear.dat
;args[1] is the base of simlated processed Hspec file
;        e.g. simLC_realizationsimLC_N1024_beta1.00m300u60
;        then there are ..._beta0.00_trial1_noise0.01_Hspec_linear.dat to ..._beta2.00_trialNsim_noise0.01_Hspec_linear.dat
;args[2] is the Nsim

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
arg_arr = STRMID(args[0],0,strlen(args)-4)

set_plot, 'ps'
device,filename='simhht_'+arg_arr[0]+'_suf.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

device, helvetica=1
device, isolatin1=1
;set_plot,'x'

device, decomposed=1
index=100;

Red=GETCOLOR('red',/true)
Blue=GETCOLOR('blue',/true)
Green=GETCOLOR('green',/true)
Cyan=GETCOLOR('cyan',/true)
Magenta=GETCOLOR('magenta',/true)
Yellow=GETCOLOR('yellow',/true)
Pink=GETCOLOR('pink',/true)
Orange=GETCOLOR('orange',/true)
Gray=GETCOLOR('gray',/true)
Navy=GETCOLOR('navy',/true)

psize=0.6
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;

; have to specify bin width in the same unit as LC time , use days!
;binwidth=10.D/(24.*60.)
;binwidth=1.D
binwidth=50.D/(86400.D)

;1. Set a range of power-law index beta for PSD
;   from 0 to 2.5 in 0.1 steps
;Nbeta=26
;beta=(findgen(Nbeta))/10.D

Nbeta=20
beta=(findgen(Nbeta)+1)/10.D

;0.5 to 1.5 in 0.05 steps
;Nbeta=20
;beta=(findgen(Nbeta)+10)/20.D

;   from 0.3 to 1.5 in 0.01 steps
;Nbeta=121
;beta=(findgen(Nbeta)+30)/100.D

;2. Read the Hspec from real data
readcol,args[0],obsfreq,obshs1,obshs2,obshs3,obshs4,obshs5,obshs6,obshs7,obshs8
obshstotal=obshs1+obshs2+obshs3+obshs4+obshs5+obshs6+obshs7+obshs8
Nobs=n_elements(obshstotal)
obsfreq=obsfreq[1:Nobs-2]
obshstotal=obshstotal[1:Nobs-2]

; specify the number of simulated LCs
Nsim=args[2]

SuFfname='simhht_suf_'+arg_arr[0]

SuF=dblarr(2,Nbeta);
SuF[0,*]=beta

hspecsim=dblarr(2,Nsim,Nobs-2)

simChisq=dblarr(Nsim)
obsChisq=0.D

simhspecavg=dblarr(Nobs-2)
simhspecstdev=dblarr(Nobs-2)

for j=0,Nbeta-1,1 do begin
  for jj=1,Nsim,1 do begin
    readcol,args[1]+'_beta'+string(beta[j],format='(f0.2)')+'_trial'+strtrim(string(jj,format='(i)'),1)+'_noise0.01_Hspec_linear.dat',freq,hs1,hs2,hs3,hs4,hs5,hs6,hs7,hs8
    hstotal=hs1+hs2+hs3+hs4+hs5+hs6+hs7+hs8
    N=n_elements(hstotal)
    ;freq=freq[1:N-2]
    ;hstotal=hstotal[1:N-2]
    freq=freq[1:Nobs-2]
    hstotal=hstotal[1:Nobs-2]
    print,'N is ',N, ' and Nobs is ',Nobs
    ;print,'freq is ', freq, ' and freqobs is ',obsfreq
    ;print,'hspec total is ', hstotal, ' and hspecobs is ',obshstotal
    hspecsim[0,jj-1,*]=freq
    hspecsim[1,jj-1,*]=hstotal
  endfor

  RESOLVE_ROUTINE, 'OBSSIMDISTFUNC2', /IS_FUNCTION
  meanfreq=10^mean(alog10(freq))
  simnorm=mean(hspecsim[1,*,*])/mean(obshstotal)
  X=obssimdistfunc2(hspecsim, obshstotal,simnorm,Nobs-2)
  simnorm=X.simnorm


  for ii=1,N-2,1 do begin
    simhspecavg[ii]=mean(hspecsim[1,*,ii])
    simhspecstdev[ii]=stddev(hspecsim[1,*,ii])
    obsChisq+=(obshstotal[1,ii]-simhspecavg[ii]*simnorm)^2/(simhspecstdev[ii]*simnorm*simhspecstdev[ii]*simnorm)
    for jj=0,Nsim-1,1 do begin     
      simChisq[jj]+=(hspecsim[1,jj,ii]-simhspecavg[ii])^2/(simhspecstdev[ii]*simhspecstdev[ii])
      ;print,'adding up simChisq for sim ',ii,' and freq ',i, ' :  ',simChisq[ii]
    endfor
    ;print,'adding up simChisq for freq ',ii,' : ',simChisq[0]
  endfor

  m=0
   
  for i=0,Nsim-1,1 do begin
     if obsChisq lt simChisq[i] then m++;
     print,obsChisq, simChisq[i], m
  endfor
  SuF[1,j]=double(m)/double(Nsim)
  print,'success fraction for beta ',beta[j],' is: ',SuF[1,j]

endfor

plot,SuF[0,*],SuF[1,*],xtitle='Power-law index of Hilbert spectrum', ytitle='Success Fraction',font=0,charsize=1.2
oplot,SuF[0,*],SuF[1,*],color=Red

openw,lun,SuFfname+'.dat',/get_lun
printf,lun,SuF
Free_lun,lun

device, /close_file

end

