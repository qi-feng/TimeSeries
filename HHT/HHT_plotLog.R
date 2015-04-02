library(hht);

args <- commandArgs(TRUE)
#filename="LC1-3_00130_short.dat"
#filename="data/simLC_beta1.00000.dat"
#filename="hht/simLC_beta1_m300u10.dat"
#filename="hht/simLC_N512_beta1.00000.dat"
#filename="simLC_N512_beta1_10m300u30.dat"
#filename="simLC_N1024_beta1.00m300u60.dat"
#filename="stepLC_N512.dat"

filename=args[1]

#my_hht_func <- function(filename ){
nimf <- 10
trials <- 10
#imf.list <- 1:9
imf.list <- 1:6
#imf.list <- 1:4
#noise.amp <- 0
#noise.amp <- 1.e-2
#noise.amp <- 1.e-8
noise.amp <- 1

trials.dir <- "hht_trials"
#freq.span <- c(1.e-9, 1.e2)

#filesplit<-unlist(strsplit(filename,"[.]"))
#filebase=filesplit[1]
filebase=sub("[.][^.]*$", "",filename)
filebase=paste0(filebase,"_noise",noise.amp)

print(filename)
print(filebase)

#testsimlc1<-scan(filename)
#sig=testsimlc1

strs <- readLines(filename)
dat <- read.table(text=strs,skip=0,header=FALSE)
sig=dat[,2]
tn=dat[,1]*86400.
N=length(sig)

#sig=dat[,1]
#tn=seq(0,N-1)

set.seed(628)

#EEMD(testsimlc1,tn , noise.amp, trials, nimf, trials.dir = trials.dir)
EEMD(sig,tn , noise.amp, trials, nimf, trials.dir = trials.dir)

EEMD.result <- EEMDCompile(trials.dir, trials, nimf)     #Plot the IMFs

dt=mean(diff(tn))
time.span <- c(tn[1], tn[N])
freq.span <- c(1./(tn[N]-tn[1]), 1./dt)
#freq.span <- c(log10(1./(tn[N]-tn[1])), log10(1./dt))
#freq.span <- NULL

os <- TRUE
res <- TRUE

IMFname=paste0(filebase,"_IMF_Rmed.pdf")
pdf(IMFname,width=10,height=8)
PlotIMFs(EEMD.result, time.span,imf.list, os, res)
#PlotIMFs(EEMD.result, time.span, os, res)
dev.off();

#ft <- list()
#ft$nfft <- N
#ft$ns <- 30
#ft$nov <- 29

amp.span=NULL

##fgram <- FTGramImage(sig, dt, ft, time.span = time.span, freq.span = freq.span,amp.span = amp.span, pretty = TRUE, img.x.format = "%.1f",img.y.format = "%.0f",main = "Port Foster Event - Fourier Spectrogram")

hres=EEMD.result
dfreq=2./(tn[N]-tn[1])

#hgram <- HHRender(hres, dt, dfreq, time.span = NULL, freq.span = NULL, scaling = "log", combine.imfs = TRUE, verbose = TRUE)
hgram <- HHRender(hres, dt, dfreq, time.span = time.span, freq.span = freq.span, scaling = "none", combine.imfs = TRUE, verbose = TRUE)

HSname=paste0(filebase,"_HS_log_Rmed_whiteBack.pdf")
pdf(HSname,width=10,height=8)
HHGramImage(hgram, backcol=c(255,255,255),scaling="log")
dev.off();


hspec <- HHSpectrum(hres,dfreq,scaling = "none")

Hspecname=paste0(filebase,"_Hspec_log_Rmed.pdf")
pdf(Hspecname,width=10,height=8)
HHSpecPlot(hspec, scaling="log", show.imfs = TRUE,show.total = TRUE,freq.span=freq.span)
#, scaling = "log")
dev.off();

Hspecdat=paste0(filebase,"_Hspec_log.dat")
hspecdata<-data.frame(freq=hspec$freq,amplitude=hspec$amplitude)
write.table(hspecdata,Hspecdat,row.names = FALSE,col.names = TRUE)

#return(hgram)
#}


