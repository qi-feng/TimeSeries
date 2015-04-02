__author__ = 'qfeng'
import process_s6_root
from matplotlib import pyplot as plt
import numpy as np

f = '/Users/qfeng/Data/veritas/BLLac20110628s6.root'
fname = '/Users/qfeng/Data/Blazars/BLLac/LightCurves/flareLC'

process_s6_root.get_run_list(f)
#t, f = process_s6_root.process_s6_root(f,binwidth_min=4.,kernel_bandwidth_min=2., compfile=fname, use_mjd=False,
#                                       plot_hist=True, plotname='/Users/qfeng/Dropbox/Astro/thesis/QiDissertation_HEAD_skeleton/images/BLLacKDE.eps')
t, f, th, fh, tt, ff, ww = process_s6_root.process_one_run(f,57432, binwidth_min=1.,kernel_bandwidth_min=1., compfile=fname, use_mjd=False,
                                       ea_method='x',plot_hist=True, doplot=True, bb=True, p0=0.01)
#t_list, c_list, \
#th_list, ch_list, \
#t_bb, h_bb, bb_hw_list = process_s6_root.process_one_run(
#    f,57432, binwidth_min=4.,kernel_bandwidth_min=2., p0=0.01,
#    compfile=fname, use_mjd=True, plotname='BLLac_KDE_BB_p0_0_01.eps',
#    ea_method='avg',plot_hist=True, doplot=False, bb=True, spectrum=False #, ene_bw=100.
#)

#t_list1, c_list1, \
#th_list1, ch_list1, \
#t_bb1, h_bb1, bb_hw_list1 = process_s6_root.process_one_run(
#    f,57433, binwidth_min=1.,kernel_bandwidth_min=1.,
#    compfile=fname, use_mjd=True, t_start=
#    ea_method='avg',plot_hist=True, doplot=False, bb=True
#)
#
#t_list = np.concatenate((t_list, t_list1))
#c_list = np.concatenate((c_list, c_list1))
#th_list = np.concatenate((th_list, th_list1))
#ch_list = np.concatenate((ch_list, ch_list1))
#t_bb = np.concatenate((t_bb, t_bb1))
#h_bb = np.concatenate((h_bb, h_bb1))
##bb_hw_list = np.concatenate((bb_hw_list + bb_hw_list1))

#process_s6_root.plot_list(t_list, c_list, th_list, ch_list, np.diff(th_list), t_start=None, t_stop=None,
#              ax = None, use_mjd = False, compfile = fname, plotname='BLLac_KDE_BB.eps',
#              plot_hist = True, kde_label = None, hist_label = None, bb=True, t_bb=t_bb, h_bb=h_bb)

#process_s6_root.process_s6_list(f, binwidth_min=1.,kernel_bandwidth_min=1., compfile=fname, use_mjd=False,
#                                       ea_method='avg',plot_hist=True, doplot=True)

#process_s6_root.process_s6_list(f, energy_GeV=200., binwidth_min=4.,kernel_bandwidth_min=2.,use_mjd=False,
#                                compfile=fname,doplot=True,plot_hist=True)


t_2min = []
f_2min = []
df_2min = []
with open(fname) as fp:
    for line in fp:
        t_2min.append(float(line.split()[0]))
        f_2min.append(float(line.split()[1]))
        df_2min.append(float(line.split()[2]))

t_2min.append(55740.453072)
f_2min.append(2.70743e-7)
df_2min.append(1.046446e-7)

#subplot = 111
#fig = plt.figure(figsize=(8, 5))
#fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
#ax = fig.add_subplot(subplot)

#ax.plot(t,f)
#plt.errorbar(np.array(t_2min), np.array(f_2min), xerr=120./86400., yerr=np.array(df_2min), color='blue', fmt='o',label="Official VERITAS")
#ax.errorbar(np.array(t_2min[-1]), np.array(f_2min[-1]), xerr=480./86400., yerr=np.array(df_2min[-1]), color='blue', fmt='o')

#ax.set_xlabel('t (MJD)')
#ax.set_ylabel('flux (m$^{-2}$ s$^{-1}$)')
#ax.legend(loc='upper right')

#fig.show()
plt.show()
