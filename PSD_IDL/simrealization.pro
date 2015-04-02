;Simulate longer light curves with finer sampling rate assuming
;power law spectrum density: S(f)~(1/f)^beta (Timmer & Konig 1995)
;then rebin, apply window function in the fashion of an imput real LC
;then get PSDs and compare sim and real by calculating success fraction (SuF),
;following the PSRESP method described in Chatterjee 2008
COMPILE_OPT idl2, HIDDEN

;provide lc file name in command line, 
;e.g. idl -e ".run psd_suf.pro" -args LC.dat
;LC.dat has 3 columns: time, flux, flux error
;change the value of binwidth below correspondingly 

;args[0] is the LC file name
;args[1] is beta, power-law index
args = command_line_args()
;arg_arr = STRSPLIT(args, '.', ESCAPE='\', /EXTRACT)
arg_arr = STRMID(args,0,strlen(args)-4)
lcfile=strtrim(args[0],1)
print, args
print,lcfile

; have to specify bin width in the same unit as LC time , use days!
;binwidth=10.D/(24.*60.)
;binwidth=1.D
binwidth=50.D/(86400.D)

beta=double(args[1])

readcol,args[0],F='D,D',realtime,realflux;

; duration of real LC
tspan=0L
binratio=0; so that sim LC has smaller binwidth/2^3 

RESOLVE_ROUTINE, 'LWA_func', /IS_FUNCTION

LWAstruct=LWA_func( realtime, realflux, binwidth, binratio,beta )
LCsim=dblarr(2,LWAstruct.srN);
LCsim[0,*]=LWAstruct.srt
LCsim[1,*]=LWAstruct.srx

openw,lun,'simLC_realization'+arg_arr[0]+'_beta'+string(beta,format='(f0.2)')+'.dat',/get_lun
printf,lun,LCsim
Free_lun,lun

end

