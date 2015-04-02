;set_plot,'x'
set_plot, 'ps'
device,filename='M4_PSD.eps', /encapsulated, xsize=8 ,ysize=5,/inches

psize=0.6
a=findgen(16)*(!pi*2/16.)
usersym, psize*cos(a),psize*sin(a), /fill

!p.thick = 5;
!x.thick = 5;
!y.thick = 5;
!z.thick = 5;

readcol,'XRTLC_nightly_all.txt',F='D,D,D,D,D',xrtMJD,dxrtMJD_pos,dxrtMJD_neg,xrtR,dxrtR
readcol,'NightlyVERLC.txt',F='D,D,D',VLtime,Vflux,dVflux;
readcol,'whipple95-09.txt',F='D,D,D',Wmjd,Wrate,dWrate;

readcol,'Mrk421_opticalPol_all.txt',F='D,D,D,D,D,D,D,D,D',opT,op2,op3,op4,op5,opFrac,dopFrac,opAngle,dopAngle
readcol,'Mrk421_opticalFlux_all.txt',F='D,D,D,D,D',oT,o2,o3,oMagn,doMagn
readcol,'test.dat',F='D,D',testF,dtestF

onemjd=86400.
METMJDREF=51910.000742870368

readcol,'OVRO.dat',F='D,D,D',Rmjd,Rflux,dRflux
readcol,'LATweekly.dat',F='D,D,D,D,D',Lstart,Lstop,LATf,dLATf,LFrac;
LATmjd=(Lstart+Lstop)/(2.0*onemjd)+METMJDREF

;RESOLVE_ROUTINE, 'NVA_func', /IS_FUNCTION

mL=mean(LATf)
   NL=n_elements(LATf)
resL=LATf-mL
fftL=fft(resL)

print,fftL

print,'num of points: ',NL

psdL=ABS(fftL)^2
print,psdL

toffset=LATmjd-LATmjd[0]
sampleT=7.*24.*3600.
freqL=(toffset/7.)/(NL*sampleT)

freqL[0]=freqL[1]
psdL[0]=psdL[1]
plot,freqL,psdL,psym=8,/xstyle,/xlog,/ystyle,/ylog,yrange=[1e-22,1e-14],xrange=[1e-9,1e-5],xtitle='!17Frequency (Hz)',ytitle='!17Power',charsize=1.2,font=0


device, /close_file

end

