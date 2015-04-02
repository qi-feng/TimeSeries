function obssimdistfunc2,PSDsim, realPSD, simnorm,numberofsimfreq
COMPILE_OPT idl2, HIDDEN

   nsim=50
   ;search from 1/range*simnorm to range*simnorm
   range=2
   simnormlist=dblarr(2,nsim)
   for j=0,nsim-1,1 do begin
      simnormlist[1,j]=simnorm*(1./range+(range-1./range)/nsim*j)
      simPSDavg=dblarr(numberofsimfreq)
      simPSDstdev=dblarr(numberofsimfreq)
      for i=0,numberofsimfreq-1,1 do begin
        simPSDavg[i]=mean(PSDsim[1,*,i]*simnormlist[1,j])
        simPSDstdev[i]=stddev(PSDsim[1,*,i]*simnormlist[1,j])
        simnormlist[0,j]+=(realPSD[i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
        ;simnormlist[0,j]+=(realPSD[1,i]-simPSDavg[i])^2/(simPSDstdev[i]*simPSDstdev[i])
      endfor
   endfor
   x=min(simnormlist[0,*],index)
   outstruct = {                  $
          simnorm:  simnormlist[1,index],         $
          minobschisq: simnormlist[0,index]  };,         $
          ;dpsd:  err }
   return, outstruct
end

