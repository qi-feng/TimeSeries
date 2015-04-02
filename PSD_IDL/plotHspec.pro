
;args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;infile=strtrim(args)
;print, args[0]
;print,infile

fname='simLC_N512_beta1.00000_noise0.01_Hspec_log'
set_plot, 'ps'
device,filename=fname+'replot.eps', /encapsulated, xsize=8 ,ysize=5,/inches,/color,bits_per_pixel=8

device, helvetica=1
device, isolatin1=1
set_plot,'x'

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

;readcol,args[0],F='D,D,D',freq,psd,dpsd;
readcol,fname+'.dat',freq,hs1,hs2,hs3,hs4,hs5,hs6,hs7
hstotal=hs1+hs2+hs3+hs4+hs5+hs6+hs7

N=n_elements(hstotal)
binfactor = 20.
m = long(1.0D* (N/2L) / binfactor)
nbinned = m*binfactor
hscut = hstotal[0:nbinned-1]
hsbinned = rebin(hscut,m)
dhsbinned = sqrt(1.0/binfactor)*hsbinned

ploterror,freq,hsbinned,dhsbinned,yrange=[min(hstotal)*0.2,max(hstotal)*2],psym=8,/nohat,/xlog,/ylog,xstyle=3,ystyle=3; color=Orange,ERRCOLOR=Orange

;ploterror,freq,hstotal,yrange=[min(hstotal)*0.2,max(hstotal)*2],psym=8,/nohat,/xlog,/ylog,xstyle=3,ystyle=3; color=Orange,ERRCOLOR=Orange
;yrange=[1.e2,1.e8]

end

psdfit=LINFIT(alog10(freq),alog10(psd),CHISQ=psdChisq)
oplot,freq,10^((psdfit[0])+(psdfit[1])*alog10(freq)),thick=5,linestyle=3,color=Cyan

xyouts,[max(freq)*0.1,max(freq)*0.1],[max(psd),max(psd)],'index:'+string(psdfit[1])

print,'index is: ',psdfit[1]
device, /close_file

end

