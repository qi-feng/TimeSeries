;This program plots a rebinned psd from input LC
;provide lc file name in command line, 
;e.g. idl -e ".run psd.pro" -args LC.dat
;LC.dat has 3 columns: time, flux, flux error
;change the value of binwidth below correspondingly 
args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;assume a .xxx file appendix
arg_arr = STRMID(args,0,strlen(args)-4)

lcfile=strtrim(args[0],1)
print, args
print,lcfile

set_plot,'ps'
device,filename='PSD_'+arg_arr[0]+'.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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

;1. Read a real LC, careful with unit in LC (assuming days)
;readcol,args[0],F='D,D,D',VLtime,Vflux,dVflux;
readcol,args[0],F='D,D',VLtime,Vflux;

; have to specify bin width in the same unit as LC time
binwidth=(50.D)/(86400.D)
; duration of real LC
tspan=0L
binratio=3; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

RESOLVE_ROUTINE, 'suf_normpsd_func', /IS_FUNCTION

sufoutstruct=suf_normpsd_func( VLtime, Vflux, binwidth, binratio,1,1 )


Nout=n_elements(sufoutstruct.realFreqrebin)
outPSD=dblarr(3,Nout);
outPSD[0,*]=sufoutstruct.realFreqrebin
outPSD[1,*]=sufoutstruct.realPSDrebin
outPSD[2,*]=sufoutstruct.drealPSDrebin

openw,lun,arg_arr[0]+'_PSD.dat',/get_lun
printf,lun,outPSD
Free_lun,lun

ploterror,sufoutstruct.realFreqrebin , sufoutstruct.realPSDrebin, sufoutstruct.drealPSDrebin, psym=8,yrange=[min(sufoutstruct.realPSDrebin)*0.2,max(sufoutstruct.realPSDrebin)*2],/nohat,/xlog,/ylog,font=0,charsize=1.2,xtitle='!17 Frequency (Hz)',ytitle='!17 RMS normalized power'; color=Orange,ERRCOLOR=Orange


psdfit=LINFIT(alog10(sufoutstruct.realFreqrebin),alog10(sufoutstruct.realPSDrebin),CHISQ=psdChisq)
oplot,sufoutstruct.realFreqrebin,10^((psdfit[0])+(psdfit[1])*alog10(sufoutstruct.realFreqrebin)),thick=5,linestyle=3,color=Cyan

xyouts,[max(sufoutstruct.realFreqrebin)*0.1,max(sufoutstruct.realFreqrebin)*0.1],[max(sufoutstruct.realPSDrebin),max(sufoutstruct.realPSDrebin)],'index:'+string(psdfit[1]),font=0,charsize=1.2

print,'index is: ',psdfit[1]

device, /close_file

end

