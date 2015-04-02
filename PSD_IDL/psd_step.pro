;This program plots a rebinned psd from input LC
;provide lc file name and psd index in command line, 
;e.g. idl -e ".run psd_test.pro" -args LC.dat 1.0
;LC.dat has 2 columns: time, flux;;;, flux error
;change the value of binwidth below correspondingly 

;arg 0: LC filename
;arg 1: power law index of PSD

args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;lcfile=strtrim(args[0],1)
;print, args
;print,lcfile

fname = strtrim(STRMID(args[0],0,strlen(args)-4),1)

;end
;fname='simLC_N1024_beta1.00m300u60'
;fname='simLC_N1024_beta1.10m300u30'

readcol,args[0],F='D,D',VLtime,Vflux;
;beta=2.

N=n_elements(VLtime)

set_plot,'ps'
device,filename=fname+'PSDtest.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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


;VLtime=VLtime/86400.D
; have to specify bin width in the same unit as LC time
binwidth=(50.D)/(86400.D)
;binwidth=1.
; duration of real LC
tspan=0L
binratio=0; so that sim LC has smaller binwidth/2^3 
NV=n_elements(VLtime)
tspan=(VLtime[NV-1]-VLtime[0])/binwidth

!p.multi=[0,1,2]
!p.position=[0.15,0.57,0.95,0.97]

plot,VLtime*86.400,Vflux,ystyle=3,/xstyle,xtitle='!17Time (ks)',ytitle='!17Flux (arbitary unit) ',font=0;

!p.position=[0.15,0.08,0.95,0.48]

RESOLVE_ROUTINE, 'PSD_fft_func', /IS_FUNCTION
freq=dblarr(N/2)
psd1=dblarr(N/2)

psdbool=PSD_fft_func(VLtime,Vflux,sqrt(Vflux),binwidth,freq,psd1)

;RESOLVE_ROUTINE, 'suf_plotpsd_func', /IS_FUNCTION
;RESOLVE_ROUTINE, 'powspec', /IS_FUNCTION

;sufoutstruct=suf_plotpsd_func(VLtime, Vflux, binwidth, binratio,beta,3 )
;powspecoutstruct=powspec(VLtime*86400.D, Vflux)


plot,freq,psd1,/xstyle,/xlog,/ystyle,/ylog,xrange=[min(freq)*0.5,max(freq)*2],yrange=[min(psd1)*5e-2,max(psd1)*2.],xtitle='!17Frequency (Hz)',ytitle='!17RMS normalized power',font=0,xtickformat='exponent',ytickformat='exponent'

device, /close_file

end

