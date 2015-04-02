;Simulate longer light curves with finer sampling rate assuming
;broken power law spectrum density: S(f)~(1/f)^beta1 when f < Fbkn
;or above Fbkn S(f)~(1/f)^beta2 (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real by calculating success fraction (SuF),
;following the PSRESP method described in Chatterjee 2008

COMPILE_OPT idl2, HIDDEN

;provide lc file name in command line, 
;e.g. idl -e ".run psd_suf.pro" -args LC.dat
;LC.dat has 3 columns: time, flux, flux error
;change the value of binwidth below correspondingly 
args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;assume a .xxx file appendix
arg_arr = STRMID(args,0,strlen(args)-4)

lcfile=strtrim(args[0],1)
print, args
print,lcfile

set_plot, 'ps'
device,filename='simbkn_'+arg_arr[0]+'_suf.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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
;beta2 from 0.2 to 3.2 in 0.2 steps
;Nbeta2=15
;beta2=(findgen(Nbeta2)+8)/10.D

;   from 0.5 to 2.5 in 0.1 steps
Nbeta2=20
beta2=(findgen(Nbeta2)+5)/10.D

;beta1 from 0.5 to 1.0 in 0.1 steps
;Nbeta1=9
;beta1=(findgen(Nbeta1)+6)/10.D

Nbeta1=1
beta1=1.0D

;breaking frequency Fbkn
;from 0.5 to 3
;NFbkn=10
;Fbkn=(findgen(NFbkn)+1)*0.5e-6
;Fbkn[0]=1.e-8
;Fbkn[1]=2.e-8
;Fbkn[2]=5.e-8
;Fbkn[3]=1.e-7
;Fbkn[4]=2.e-7
;Fbkn[5]=5.e-7
;Fbkn[6]=1.e-6
;Fbkn[7]=2.e-6
;Fbkn[8]=3.e-6
;Fbkn[9]=5.e-6

;from 0.00025 to 0.005 in steps of 0.00025
NFbkn=20
Fbkn=(findgen(NFbkn)+1)*0.00025
;NFbkn=7
;Fbkn=(findgen(NFbkn)+1)*0.01
;Fbkn[0]=0.0003
;Fbkn[1]=0.0005
;Fbkn[2]=0.001
;Fbkn[3]=0.0015
;Fbkn[4]=0.002
;Fbkn[5]=0.004
;Fbkn[6]=0.007


;2. Read a real LC, careful with unit in LC (assuming days)
;readcol,lcfile,F='D,D,D',VLtime,Vflux,dVflux;
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
readcol,args[0],F='D,D',VLtime,Vflux;

;readcol,'lcCORexploose20min.txt',F='D,D,D',VLtime,Vflux,dVflux

; duration of real LC
tspan=0L
binratio=0; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

; specify the number of simulated LCs
Nsim=100

; specify a filename base to store suf values
SuFfname='simbkn_suf_'+arg_arr[0]

; SuF=dblarr(3,Nbeta1,Nbeta2);
; for i=0,Nbeta2-1,1 do begin
;   SuF[0,*,i]=beta1
; endfor
; for i=0,Nbeta1-1,1 do begin
;   SuF[1,i,*]=beta2
; endfor

SuF=dblarr(3,Nbeta2,NFbkn);
for i=0,NFbkn-1,1 do begin
  SuF[0,*,i]=beta2
endfor
for i=0,Nbeta2-1,1 do begin
  SuF[1,i,*]=Fbkn
endfor

;RESOLVE_ROUTINE, 'SuF_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'suf_brokenpsd_func', /IS_FUNCTION

for i=0,Nbeta1-1,1 do begin
  device,filename='suf2_beta1_'+strtrim(beta1[i],1)+'.eps',/encapsulated,xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8
  device, helvetica=1
  device, isolatin1=1
  openw,lun1,SuFfname+strtrim(beta1[i],1)+'_message.txt',/get_lun

 for I_F=0,NFbkn-1,1 do begin
   for j=0,Nbeta2-1,1 do begin
   sufoutstruct=suf_brokenpsd_func( VLtime, Vflux, binwidth, binratio,beta1[i],beta2[j],Fbkn[I_F], Nsim )
   SuF[2,j,I_F]=sufoutstruct.suf
   
    ;SuF[2,j,I_F]=suf_brokenpsd_func( VLtime, Vflux, binwidth, binratio,beta1[i],beta2[j],Fbkn[I_F], Nsim )
     ;SuF[2,i,j]=suf_brokenpsd_func( VLtime, Vflux, binwidth, binratio,beta1[i],beta2[j],Fbkn[I_F], Nsim ) 

     print,'success fraction for bkn power law at Freq_bkn: ',Fbkn[I_F],' with beta1: ',beta1[i] ,' and beta2: ',beta2[j],' is: ',SuF[2,j,I_F]

     printf,lun1,'success fraction for bkn power law at Freq_bkn: ',Fbkn[I_F],' with beta1: ',beta1[i] ,' and beta2: ',beta2[j],' is: ',SuF[2,j,I_F]
   endfor
  endfor
  Free_lun,lun1

  SURFACE,SuF[2,*,*],SuF[0,*,0],SuF[1,0,*],XST = 1, YST = 1, ZST = 1, xtitle='beta2', ytitle='Fbkn',/LEGO,/ylog;,/NOERASE


  ;oplot,SuF[0,*],SuF[1,*],color=Red

  openw,lun,SuFfname+strtrim(beta1[i],1)+'.dat',/get_lun
  printf,lun,SuF
  Free_lun,lun

  device, /close_file
endfor

device,filename='simpsd_'+arg_arr[0]+'_PSD.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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

ploterror,sufoutstruct.realFreqrebin , sufoutstruct.realPSDrebin, sufoutstruct.drealPSDrebin, psym=8,yrange=[min(sufoutstruct.realPSDrebin)*0.2,max(sufoutstruct.realPSDrebin)*2],xrange=[min(sufoutstruct.realFreqrebin)*0.5,max(sufoutstruct.realFreqrebin)*10],/nohat,/xlog,/ylog,font=0,charsize=1.2; color=Orange,ERRCOLOR=Orange

psdfit=LINFIT(alog10(sufoutstruct.realFreqrebin),alog10(sufoutstruct.realPSDrebin),CHISQ=psdChisq)
oplot,sufoutstruct.realFreqrebin,10^((psdfit[0])+(psdfit[1])*alog10(sufoutstruct.realFreqrebin)),thick=5,linestyle=3,color=Cyan

xyouts,[max(sufoutstruct.realFreqrebin)*0.8,max(sufoutstruct.realFreqrebin)*0.8],[max(sufoutstruct.realPSDrebin)*0.2,max(sufoutstruct.realPSDrebin)*0.2],'input LC index:'+string(psdfit[1],format='(f0.3)'),font=0
print,'index is: ',psdfit[1]

device, /close_file

end

