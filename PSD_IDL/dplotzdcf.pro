args = command_line_args()
arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
dcffile=args[0]

print,arg_arr[0]
print,dcffile

;set_plot,'x'
set_plot, 'ps'
device,filename=arg_arr[0]+'_ZDCF.eps', /encapsulated, xsize=7,ysize=5,/inches,/color,bits_per_pixel=8

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;

;Red=GETCOLOR('red',/true)
;Blue=GETCOLOR('blue',/true)
;Green=GETCOLOR('green',/true)
;device, decomposed=1
;index=100;

psize=1.
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

;readcol,'VX.dcf',VXTau,VXNsigTau,VXPsigTau,VXdcf,NdVXdcf,PdVXdcf,VXbin
readcol,dcffile,VXTau,VXNsigTau,VXPsigTau,VXdcf,NdVXdcf,PdVXdcf,VXbin

VXTau*=86400.
VXNsigTau*=86400.
VXPsigTau*=86400.

;readcol,'SMA_Met.dcf',smTau,smNsigTau,smPsigTau,smdcf,Ndsmdcf,Pdsmdcf,smbin
;readcol,'XO.dcf',XOTau,XONsigTau,XOPsigTau,XOdcf,NdXOdcf,PdXOdcf,XObin
;readcol,'BX.dcf',BXTau,BXNsigTau,BXPsigTau,BXdcf,NdBXdcf,PdBXdcf,BXbin

VXNsigTau=-VXNsigTau
;XONsigTau=-XONsigTau
NdVXdcf=-NdVXdcf
;NdXOdcf=-NdXOdcf
;BXNsigTau=-BXNsigTau
;NdBXdcf=-NdBXdcf

PnumX=1
PnumY=3
x1=0.1
x2=0.95
xm=x2
;xm=(x1+x2)/2.
lwrat=3.

margin=0.05
;xmarg=0.04
xmarg=0.0
y11=1-(x2-x1)/lwrat-margin
y12=1-margin
y21=1-2*(x2-x1)/lwrat-margin
y22=y11
y31=1-3*(x2-x1)/lwrat-margin
y32=y21
y41=1-4*(x2-x1)/lwrat-margin
y42=y31
y51=1-5*(x2-x1)/lwrat-margin
y52=y41
y61=1-6*(x2-x1)/lwrat-margin
y62=y51
;y71=1-7*(x2-x1)/lwrat-margin
;y72=y61

;!p.multi=[0,PnumX,PnumY]
;!p.position=[x1,y11,x2,y12] 
;!p.position=[x1,y11,xm-xmarg,y12]

;!x.tickname=' '

plot,VXTau,VXdcf,psym=8,ytitle='!17ZDCF',charsize=1.2,thick=5,font=0,xrange=[-15e3,15e3],yrange=[-1.,1.0],xstyle=1,ystyle=1,xtitle='!17 time lag (s)'

;oplot,VXTau,VXdcf,psym=8;,color=Red,thick=5
oploterror,VXTau,VXdcf,VXNsigTau,NdVXdcf,/nohat,/lobar;,color=Red,errcolor=Red,thick=5,errthick=3,psym=8
oploterror,VXTau,VXdcf,VXPsigTau,PdVXdcf,/nohat,/hibar;,color=Red,errcolor=Red,thick=5,errthick=3,psym=8
;xyouts,14,0.86,font=0,charsize=0.6,'Veritas XRT '

device, /close_file

end
