;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real by calculating success fraction (SuF),
;following the PSRESP method described in Chatterjee 2008

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;assume a .xxx file appendix
arg_arr = STRMID(args,0,strlen(args)-4)

lcfile=strtrim(args)
print, args
print,lcfile

set_plot, 'ps'
device,filename='SuF_'+arg_arr[0]+'_fit_PSDnorm.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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

;1. Set a range of power-law index beta for PSD
;   from 0 to 2.5 in 0.1 steps
Nbeta=26
beta=(findgen(Nbeta))/10.D

;   from 0.3 to 1.5 in 0.01 steps
;Nbeta=121
;beta=(findgen(Nbeta)+30)/100.D

;2. Read a real LC, careful with unit in LC (assuming days)
;readcol,lcfile,F='D,D,D',VLtime,Vflux,dVflux;
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
readcol,args[0],F='D,D',VLtime,Vflux;

;readcol,'lcCORexploose20min.txt',F='D,D,D',VLtime,Vflux,dVflux

; have to specify bin width in the same unit as LC time
binwidth=1.
; duration of real LC
tspan=0L
binratio=3; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

; specify the number of simulated LCs
Nsim=100

; specify a filename base
SuFfname='SuF_normpsd_'+arg_arr[0]

 SuF=dblarr(2,Nbeta);
 SuF[0,*]=beta

;RESOLVE_ROUTINE, 'SuF_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'suf_normpsd_func', /IS_FUNCTION

for j=0,Nbeta-1,1 do begin
 ;SuF[1,j]=SuF_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 ;SuF[1,j]=suf_normpsd_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 sufoutstruct=suf_normpsd_func( VLtime, Vflux, binwidth, binratio,beta[j], Nsim )
 SuF[1,j]=sufoutstruct.suf
 print,'success fraction for beta ',beta[j],' is: ',SuF[1,j]
endfor

plot,SuF[0,*],SuF[1,*],xtitle='Power-law index of PSD', ytitle='Success Fraction'
oplot,SuF[0,*],SuF[1,*],color=Red

openw,lun,SuFfname+'.dat',/get_lun
printf,lun,SuF
Free_lun,lun

device, /close_file

end

