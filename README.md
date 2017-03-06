# TimeSeries
Basic time series processing with Python under development. 
######
## SimPointProcesses:
Python scripts to simulate a sequence of arriving times of events. A fake energy can be assigned to each event. 
#### To simulate a list of TTEs following a Poisson noise power spectrum, use "massproduce_two_sim_poisson_seq.py". 
Edit these parameters in the file:
```
n_start
n_stop
total
total2
mean_rate
mean_rate2
Efake
Efake2
```
You can get the total number of events and the mean rate for lower and higher energies, the duration is then determined.
######
######
## simRedNoises:
####  To simulate a list of time tagged events (TTE) following a red nosie power spectrum, provide the TTE list file from real data so that the energies from the data can be used. 
Edit these parameters in the file:
```
the above file name (energy_file), 
the number of simulations (n_lc), 
the mean rate (mean_rate), 
and the duration (T) in the file "mass_rednoise_readE_2pp.py" and run it. 
```
The method generate_pl_pp in sim_rednoise_pointprocess is then called. 
The output file names are 'simRedNoiseTTE_duration'+str(T)+'_rate'+str(mean_rate)+'_withEreadFromM4Obs_trial'+str(i)+'.txt'
######
######
## PSD_IDL: 
IDL scripts to simulate 1/f noises (following Timmer & Konig 1995), calculated power spectral densities (PSD), and success fractions (following Uttley et al. 2002 and Chatterjee et al. 2008). 

Light curve file name is provided as command line arguments, e.g.: 
idl -e ".run psd_suf.pro" -args LC.dat
The file (e.g. LC.dat) has 3 columns: time, flux, flux error
######
######
## HHT (see also directory ./pedvar): 
R scripts using the hht R package implimented by Daniel C. Bowman, see 
http://cran.r-project.org/web/packages/hht/index.html

The input light curve should contain two columns, time and flux (or rate etc). The script performs Ensemble Empirical Mode Decomposition (EEMD), calculates the Hilbert spectrogram, and the marginal Hilbert spectrum, and saves the above plots into pdf files. 

An example to run the script:
Rscript HHT_plotLog.R LC.dat

######
######
## pedvar: 
Jupyter notebook of examples to use continuous wavelet transform (CWT) and Hilbert-Huang Transform (HHT) to get scalogram/spectrogram from a time series. We use CWT implemented in ObsPy; and HHT in R (see above), which is run in python. 

Simple examples of composite sinusoidal waves, a chirp signal, and LIGO gravitational signal from GW150914 are shown. The VERITAS pedestal variance curve for a few pixels in one or two runs are shown. 

######
######
## VERITAS_KDE: 
process_s6_root.py:
A python class implementation to read a VEGAS stage 6 root file, 
and produces light curves in the choice of format.

Available formats: histogram, kernel density estimation, bayesian blocks

Need scikit-learn, astroML packages. 

After importing the function above, use the methods like the example here:

process_s6_root.process_one_run(f, 57432, binwidth_min=1. ,kernel_bandwidth_min=1., compfile=fname, use_mjd=False, ea_method='x', plot_hist=True, doplot=True, bb=True, p0=0.01) 

where f is the stage 6 root file name, 
57432 is the run number, 
binwidth_min and kernel_bandwidth_min are the bin width for histogram and bandwidth for kernel density estimation, respectively, 
compfile is an optional text file containing a light curve that was made by other VERITAS software, which is used for results comparison, 
if use_mjd is True the output plot will use MJD as unit for time, otherwise will use seconds since the beginning of the observations, 
if ea_method is 'avg' it uses the average effective area over each run to convert count rate into flux, otherwise it uses the effective area for each event (not implemented for KDE so far), 
plot_hist determines if plotting histograms or not, 
doplot determines if do plot at all, 
bb determines if doing Baysian blocks or not, 
p0 is the prior probability that governs the number of change points (therefore the number of bins) in Baysian blocks. 

######
