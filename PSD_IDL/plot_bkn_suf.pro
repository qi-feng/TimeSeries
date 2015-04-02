 set_plot, 'ps'
  device,filename='suf2_beta1_1.0_beta2_1.5_Fbkn0.002_replot.eps',/encapsulated,xsize=10 ,ysize=8,/inches,/color,bits_per_pixel=8
  device, helvetica=1
  device, isolatin1=1
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

!p.position = [0.1,0.1,0.9,0.9]
  NFbkn=20
  Nbeta2=20
  SuF=dblarr(3,20,20);
  b2=dblarr(20,20)
  fb=dblarr(20,20)
  suf2=dblarr(20,20)
  readcol,'simbkn_suf_simLCbkn_N4096_beta11.00_beta21.50m300u30_Fbkn0.0021.0000000.dat',F='D,D,D',b2,fb,suf1;

;;for i=0,NFbkn-1,1 do begin
;;  SuF[0,*,i]=fb[i]
;;endfor
;;for i=0,Nbeta2-1,1 do begin
;;  SuF[1,i,*]=Fbkn[i]
;;endfor
;print, "b2", b2,"fb", fb,"suf",suf1
 for I_F=0,NFbkn-1,1 do begin
   for j=0,Nbeta2-1,1 do begin
     SuF[2,j,I_F]=suf1[j+20*I_F]
     SuF[0,j,I_F]=b2[j+20*I_F]
     SuF[1,j,I_F]=fb[j+20*I_F]
     print, 'success fraction for bkn power law at Freq_bkn: ',SuF[1,j,I_F],' with beta1: 1.0 and beta2: ',SuF[0,j,I_F],' is: ',SuF[2,j,I_F]
   endfor
  endfor

  ;SURFACE,SuF[2,*,*],SuF[0,*,*],SuF[1,*,*],XST = 1, YST = 1, ZST = 1, xtitle='beta2', ytitle='Fbkn',/LEGO,/ylog,font=0,charsize=1.4, ztitle='SuF';,/NOERASE

  SURFACE,SuF[2,*,*],SuF[0,*,*],SuF[1,*,*],$
    XST = 2, YST = 2, ZST = 1, xtitle='!9b !d!17 2!C ',$
    ytitle='!9n !d!17b !N!17(Hz)!C',zrange=[0.01,1.0], $
    /LEGO,font=0,charsize=1.5, ztitle='!17SuF';,$
    ;XTEXT_ORIENTATION=45, YTEXT_ORIENTATION=45
    ;XMARGIN=[8, 8], YMARGIN=[4, 4];,/NOERASE


end
