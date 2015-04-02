library(hht);

args <- commandArgs(TRUE)
#filename="data/testSins.dat"
filename=args[1]

#my_hht_func <- function(filename ){
nimf <- 10
trials <- 10
imf.list <- 1:3

#filesplit<-unlist(strsplit(filename,"[.]"))
#filebase=filesplit[1]
filebase=sub("[.][^.]*$", "",filename)

print(filename)
print(filebase)

testsimlc1<-scan(filename)
sig=testsimlc1
N=length(sig)

tn=seq(0,N-1)
     noise.amp <- 1e-06
     trials.dir <- "test"

     set.seed(628)

  EEMD(testsimlc1,tn , noise.amp, trials, nimf, trials.dir = trials.dir)

  EEMD.result <- EEMDCompile(trials.dir, trials, nimf)     #Plot the IMFs

    time.span <- c(0, N)
    #freq.span <- c(1e-9, 0.5)
     os <- TRUE
     res <- TRUE

     IMFname=paste0(filebase,"_IMF_R.pdf")
pdf(IMFname,width=10,height=14)
PlotIMFs(EEMD.result, time.span,imf.list, os, res)
dev.off();

dt=1

ft <- list()
ft$nfft <- N
ft$ns <- 30
ft$nov <- 29

amp.span=NULL

##fgram <- FTGramImage(sig, dt, ft, time.span = time.span, freq.span = freq.span,amp.span = amp.span, pretty = TRUE, img.x.format = "%.1f",img.y.format = "%.0f",main = "Port Foster Event - Fourier Spectrogram")



hres=EEMD.result
dfreq=1e-3



hgram <- HHRender(hres, dt, dfreq, time.span = NULL, freq.span = NULL, scaling = "none", combine.imfs = TRUE, verbose = TRUE)

     HSname=paste0(filebase,"_HS_R.pdf")
pdf(HSname,width=14,height=10)

HHGramImage(hgram)

dev.off();


hspec <- HHSpectrum(hres,dfreq,scaling = "none")

     Hspecname=paste0(filebase,"_Hspec_R.pdf")
pdf(Hspecname,width=14,height=10)

HHSpecPlot(hspec, show.imfs = TRUE,show.total = TRUE, scaling = "none")

dev.off();

print(range(c(hspec$amplitude[,1],hspec$amplitude[,2],hspec$amplitude[,3],hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3])))

LogSpecname=paste0(filebase,"_LogHSpec_R.pdf")
pdf(LogSpecname,width=14,height=10)

plot(hspec$frequency,hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3],col="magenta",pch=16,ylim=range(c(hspec$amplitude[,1],hspec$amplitude[,2],hspec$amplitude[,3],hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3])),log = "xy",main="Marginal Spectrum",xlab = "Frequency", ylab = "Power")

par(new=T)

plot(hspec$frequency,hspec$amplitude[,1],col="blue",pch=16,ylim=range(c(hspec$amplitude[,1],hspec$amplitude[,2],hspec$amplitude[,3],hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3])),log = "xy")

par(new=T)

plot(hspec$frequency,hspec$amplitude[,2],col="red",pch=16,ylim=range(c(hspec$amplitude[,1],hspec$amplitude[,2],hspec$amplitude[,3],hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3])),log = "xy")

par(new=T)                                                                    

plot(hspec$frequency,hspec$amplitude[,3],col="green",pch=16,ylim=range(c(hspec$amplitude[,1],hspec$amplitude[,2],hspec$amplitude[,3],hspec$amplitude[,1]+hspec$amplitude[,2]+hspec$amplitude[,3])),log = "xy")

dev.off()
#return(hgram)
#}


