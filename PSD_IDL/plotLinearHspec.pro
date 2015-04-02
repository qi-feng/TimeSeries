
;args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
;infile=strtrim(args)
;print, args[0]
;print,infile

fname='simLC_N1024_beta1.00m300u60_noise0.01_Hspec_linear'
;fname='out_xmm0501_05_10keV_noise0.01_Hspec_linear'

set_plot, 'ps'
device,filename=fname+'replotSmall.eps', /encapsulated, xsize=5 ,ysize=3,/inches,/color,bits_per_pixel=8

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

;readcol,args[0],F='D,D,D',freq,psd,dpsd;
readcol,fname+'.dat',freq,hs1,hs2,hs3,hs4,hs5,hs6,hs7,hs8
hstotal=hs1+hs2+hs3+hs4+hs5+hs6+hs7+hs8

N=n_elements(hstotal)
freq=freq[1:N-2]
hstotal=hstotal[1:N-2]

binfactor = 10.
m = long(1.0D* (N/2L) / binfactor)
nbinned = m*binfactor
hscut = hstotal[0:nbinned-1]
hsbinned = rebin(hscut,m)
dhsbinned = sqrt(1.0/binfactor)*hsbinned

plot,freq,hstotal,xrange=[min(freq)*0.5,max(freq)*2] ,yrange=[min(hstotal)*0.2,max(hstotal)*2],psym=8,/xlog,/ylog,xstyle=3,ystyle=3,thick=5,font=0,xtitle='!17Freq (Hz)',ytitle='!17Marginal Hilbert spectrum',xtickformat='exponent',ytickformat='exponent'; color=Orange,ERRCOLOR=Orange
oplot,freq,hstotal,psym=8, color=Orange,thick=5

;ploterror,freq,hsbinned,dhsbinned,xrange=[min(freq)*0.5,max(freq)*2] ,yrange=[min(hstotal)*0.2,max(hstotal)*2],psym=8,/nohat,/xlog,/ylog,xstyle=3,ystyle=3; color=Orange,ERRCOLOR=Orange

;ploterror,freq,hstotal,yrange=[min(hstotal)*0.2,max(hstotal)*2],psym=8,/nohat,/xlog,/ylog,xstyle=3,ystyle=3; color=Orange,ERRCOLOR=Orange
;yrange=[1.e2,1.e8]

psdfit=LINFIT(alog10(freq),alog10(hstotal),CHISQ=psdChisq)
oplot,freq,10^((psdfit[0])+(psdfit[1])*alog10(freq)),thick=5,linestyle=3,color=Cyan

xyouts,[max(freq)*0.1,max(freq)*0.1],[max(hstotal),max(hstotal)],'index:'+string(psdfit[1]),font=0

print,'index is: ',psdfit[1]
device, /close_file

end

