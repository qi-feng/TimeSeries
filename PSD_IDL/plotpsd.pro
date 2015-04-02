
args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;assume a .xxx file appendix
arg_arr = STRMID(args,0,strlen(args)-4)

infile=strtrim(args)
print, args[0]
print,infile

set_plot, 'ps'
device,filename='PSD_'+arg_arr[0]+'replot.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

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

readcol,args[0],F='D,D,D',freq,psd,dpsd;

ploterror,freq,psd,dpsd,yrange=[min(psd)*0.2,max(psd)*2],psym=8,/nohat,/xlog,/ylog,xstyle=3,ystyle=3; color=Orange,ERRCOLOR=Orange
;yrange=[1.e2,1.e8]

psdfit=LINFIT(alog10(freq),alog10(psd),CHISQ=psdChisq)
oplot,freq,10^((psdfit[0])+(psdfit[1])*alog10(freq)),thick=5,linestyle=3,color=Cyan

xyouts,[max(freq)*0.1,max(freq)*0.1],[max(psd),max(psd)],'index:'+string(psdfit[1])

print,'index is: ',psdfit[1]
device, /close_file

end

