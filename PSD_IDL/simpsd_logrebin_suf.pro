;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real by calculating success fraction (SuF),
;following the PSRESP method described in Chatterjee 2008

;provide lc file name in command line, 
;e.g. idl -e ".run psd_suf.pro" -args LC.dat
;LC.dat has 3 columns: time, flux, flux error
;change the value of binwidth below correspondingly 
args = command_line_args()
arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
lcfile=strtrim(args[0],1)
print, args
print,lcfile

set_plot, 'ps'
device,filename='sim_'+arg_arr[0]+'_suf.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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

;Nbeta=5
;beta=(findgen(Nbeta)+7)/10.D

;0.5 to 1.5 in 0.05 steps
Nbeta=20
beta=(findgen(Nbeta)+10)/20.D


;   from 0.3 to 1.5 in 0.01 steps
;Nbeta=121
;beta=(findgen(Nbeta)+30)/100.D

;2. Read a real LC, careful with unit in LC (assuming days)
;readcol,lcfile,F='D,D,D',VLtime,Vflux,dVflux;
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
readcol,args[0],F='D,D',VLtime,Vflux;

;readcol,'lcCORexploose20min.txt',F='D,D,D',VLtime,Vflux,dVflux

; duration of real LC
tspan=0L
binratio=3; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

; specify the number of simulated LCs
Nsim=100

; specify a filename base to store suf values
SuFfname='sim_suf_'+arg_arr[0]

 SuF=dblarr(2,Nbeta);
 SuF[0,*]=beta

;RESOLVE_ROUTINE, 'SuF_func', /IS_FUNCTION
;RESOLVE_ROUTINE, 'suf_normpsd_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'suf_logrebinpsd_func', /IS_FUNCTION

for j=0,Nbeta-1,1 do begin
 ;SuF[1,j]=SuF_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 ;SuF[1,j]=suf_normpsd_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 sufoutstruct=suf_logrebinpsd_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 SuF[1,j]=sufoutstruct.suf
 print,'success fraction for beta ',beta[j],' is: ',SuF[1,j]
endfor

plot,SuF[0,*],SuF[1,*],xtitle='Power-law index of PSD', ytitle='Success Fraction',font=0,charsize=1.2
oplot,SuF[0,*],SuF[1,*],color=Red

openw,lun,SuFfname+'.dat',/get_lun
printf,lun,SuF
Free_lun,lun

device, /close_file
device,filename='sim_'+arg_arr[0]+'_PSD.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

device, helvetica=1
device, isolatin1=1
;set_plot,'x'

device, decomposed=1
index=100;

Nout=n_elements(sufoutstruct.realFreqrebin)
outPSD=dblarr(3,Nout);
outPSD[0,*]=sufoutstruct.realFreqrebin
outPSD[1,*]=sufoutstruct.realPSDrebin
outPSD[2,*]=sufoutstruct.drealPSDrebin

openw,lun,'sim_'+arg_arr[0]+'_PSD.dat',/get_lun
printf,lun,outPSD
Free_lun,lun

ploterror,sufoutstruct.realFreqrebin , sufoutstruct.realPSDrebin, sufoutstruct.drealPSDrebin, psym=8,yrange=[min(sufoutstruct.realPSDrebin)*0.2,max(sufoutstruct.realPSDrebin)*2],/nohat,/xlog,/ylog,font=0,charsize=1.2; color=Orange,ERRCOLOR=Orange

psdfit=LINFIT(alog10(sufoutstruct.realFreqrebin),alog10(sufoutstruct.realPSDrebin),CHISQ=psdChisq)
oplot,sufoutstruct.realFreqrebin,10^((psdfit[0])+(psdfit[1])*alog10(sufoutstruct.realFreqrebin)),thick=5,linestyle=3,color=Cyan

xyouts,[max(sufoutstruct.realFreqrebin)*0.8,max(sufoutstruct.realFreqrebin)*0.8],[max(sufoutstruct.realPSDrebin)*0.2,max(sufoutstruct.realPSDrebin)*0.2],'index:'+string(psdfit[1])

print,'index is: ',psdfit[1]

device, /close_file

end

