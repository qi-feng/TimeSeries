__author__ = 'qfeng'
import process_s6_root
from matplotlib import pyplot as plt
import numpy as np

#f = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/results_Std20140429med0s6.root'
#fname = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/M4_Std_0429_10min_560GeV.txt'

#f = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/results_20140503Std_2s6.root'
#fname = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/M4_Std_0503_10min_316GeV.txt'
#f = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/results_Amy2013s6.root'
f = '/Users/qfeng/Data/Blazars/Mrk421/2014flare/results_20100217exps6.root'

fname = None
print process_s6_root.get_run_list(f, verbose=True)

#t_list, c_list, \
#th_list, ch_list, \
#t_bb, h_bb, bb_hw_list = process_s6_root.process_one_run(
#    f,50100, binwidth_min=4.,kernel_bandwidth_min=2.,
#    compfile=fname, use_mjd=True, plotname=None, #'M420100217_KDE_BB.eps',
#    ea_method='avg',plot_hist=True, doplot=True, bb=True, spectrum=False, ene_bw=100.
#)

#t, f = process_s6_root.process_s6_root(f,energy_GeV=560., binwidth_min=10.,kernel_bandwidth_min=10.,use_mjd=False,compfile=fname)
#t, f = process_s6_root.process_s6_root(f,energy_GeV=316., binwidth_min=10.,kernel_bandwidth_min=10.,use_mjd=False,compfile=fname)
#t, f = process_s6_root.process_s6_root(f,energy_GeV=800., binwidth_min=30.,kernel_bandwidth_min=10.,use_mjd=True,compfile=fname)

#process_s6_root.process_one_run(f, 67956, energy_GeV=800., binwidth_min=30.,kernel_bandwidth_min=10.,use_mjd=True,compfile=fname,doplot=True, bb=True)

elo = 420.
process_s6_root.process_s6_list(f, energy_GeV=elo, binwidth_min=1.,
                                kernel_bandwidth_min=1.,use_mjd=True,compfile=fname,doplot=True)

#t, f = process_s6_root.process_one_run(f,73196,energy_GeV=560., binwidth_min=10.,kernel_bandwidth_min=10.)
#t, f = process_s6_root.process_one_run(f,73285,energy_GeV=316., binwidth_min=10.,kernel_bandwidth_min=10., compfile=fname, use_mjd=False,
#                                       ea_method='avg',plot_hist=True)
plt.show()
