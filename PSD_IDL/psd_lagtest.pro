;This program plots a rebinned psd from input LC
;provide lc file name and psd index in command line, 
;LC.dat has 2 columns: time, flux;;;, flux error
;change the value of binwidth below correspondingly 

COMPILE_OPT idl2, HIDDEN

;arg 0: power law index of PSD
;arg 1: lag time

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

;fname = strtrim(STRMID(args[0],0,strlen(args)-4),1)
beta = double(args[0])
lagtime= uint(args[1])

m=300.
rms=60.
N=1024
; have to specify bin width in the same unit as LC time
binwidth=(600.D)/(86400.D)

funcSimLCbool=simLC_func(N,m,rms,beta,binwidth,x)
Vflux=real_part(x)
VLtime=findgen(N)*binwidth

lagflux=

data=dblarr(2,N)
data[0,*]=VLtime
data[1,*]=Vflux
openw,lun,'simLC_N'+strtrim(string(N,format='(i)'),1)+'_'+strtrim(string(Ngaps,format='(i)'),1)+'gaps_beta'+strtrim(string(beta,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'.dat',/get_lun
printf,lun,data
Free_lun,lun


;print,fname,beta
;fname='simLC_N1024_beta1.00m300u60'
;fname='simLC_N1024_beta1.10m300u30'

;readcol,args[0],F='D,D',VLtime,Vflux;
;beta=2.

set_plot,'ps'
;device,filename=fname+'PSDtest.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

;device,filename='simLC_N'+strtrim(string(N,format='(i)'),1)+'_'+strtrim(string(Ngaps,format='(i)'),1)+'gaps_beta'+strtrim(string(beta,format='(f0.2)'),1)+'m'+strtrim(string(m,format='(i)'),1)+'u'+strtrim(string(rms,format='(i)'),1)+'.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8


;device,filename='powspecPSD_'+arg_arr[0]+'.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

;device, helvetica=1
;device, isolatin1=1
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

;1. Read a real LC, careful with unit in LC (assuming days)
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
;readcol,args[0],F='D,D',VLtime,Vflux;
;readcol,'simLC_beta1_m300u60_6.dat',F='D,D',VLtime,Vflux;
;readcol,'pn_all_0.5_1keV_50s.txt',F='D,D,D',VLtime,Vflux,dVflux;
;readcol,'out_xmm0501_0.5_10keV.dat',F='D,D,D',VLtime,Vflux,dVflux;

;readcol,'NightlyVERLC.txt',F='D,D,D',VLtime,Vflux,dVflux;

;VLtime=VLtime/86400.D
;binwidth=1.
; duration of real LC
tspan=0L
binratio=0; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

RESOLVE_ROUTINE, 'suf_plotpsd_func', /IS_FUNCTION
RESOLVE_ROUTINE, 'powspec', /IS_FUNCTION

sufoutstruct=suf_plotpsd_func(VLtime, Vflux, binwidth, binratio,beta,3 )
powspecoutstruct=powspec(VLtime*86400.D, Vflux)


device, /close_file

end

