__author__ = 'qfeng'

import numpy as np
from matplotlib import pyplot as plt

import lc
import bispectrum
import rebin_2d_array
#from mpl_toolkits.mplot3d.axes3d import Axes3D
from sim_poisson_seq import tri_poisson_seq
from mpl_toolkits.mplot3d import axes3d
from read_lc_file import read_tfdf

total = 10000
t_0 = 400.0
t_exp = 10.0
a_exp = 20.0
mean_rate = 30.0
t_1 = 700.0
t_exp1 = 50.0
a_exp1 = 25.0
t_2 = 900.0
t_exp2 = 20.0
a_exp2 = 10.0

bin_num = 1024
#bin_num = 1024*4

# Control the type of simulated time signal
coupled_sine_wave = False
exp_rise = False
exp_fall = False
exp_sym = False
#saw_tooth = True
saw_tooth = False
saw_tooth_rise = True
#saw_tooth_rise = False

#seg_freq, seg_bisp, seg_bic, seg_biph, seg_biph2 = bispectrum.seg_bispec(sub_t, sub_x)
len_seg = 512
n_seg = 100

save_plot = True
save_plot = False
do_rebin = True
do_rebin = False
rebin_scale = 5
if do_rebin == False:
    rebin_scale = 0


t_seq = tri_poisson_seq(total, mean_rate, t_exp, a_exp)

t_hist = lc.make_hist(t_seq, bin_width=1., t_start=0.0)
t_sim = t_hist.bin_center[:bin_num]
x_sim = t_hist.hist[:bin_num]

print type(t_sim), type(x_sim)
print t_sim.shape, x_sim.shape
print len(t_sim), len(x_sim)

fig = plt.figure(figsize=(14, 8))
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)

x_list = []
#t_list = np.arange(0,180.,180./bin_num)
t_list = range(0,bin_num)

# sawtooth signal
if saw_tooth == True:
    for t in t_list:
        x_list.append(-(a_exp/t_exp)*(t-int(t/t_exp)*t_exp)+0+0.5*a_exp+60.)

# sawtooth signal rise
if saw_tooth_rise == True:
    for t in t_list:
        x_list.append((a_exp/t_exp)*(t-int(t/t_exp)*t_exp)+0+0.5*a_exp+60.)

# coupled sine waves
if coupled_sine_wave == True:
    f_nyq = 0.5
    pi = 3.14159
    #wb = 2 * 3.14159 * 0.22/f_nyq
    #wc = 2 * 3.14159 * 0.375 / f_nyq
    wb = 2 * 3.14159 * 0.12
    wc = 2 * 3.14159 * 0.18
    wd = wb + wc
    we = wb + wd
    thb = pi / 3.
    thc = pi / 12.
    thd = pi / 4.
    the = pi * 3. / 8.

    for t in t_list:
        #thb = np.random.random_sample()
        #thc = np.random.random_sample()
        #thd = thb + thc
        #x_list.append( np.cos(wb * t + thb) + np.cos(wc * t + thc) + 0.5 * np.cos(wd * t + thd)
        #             + np.cos(wb * t + thb) * np.cos(wc * t + thc)
        #)
        x_list.append( np.cos(wb * t + thb) + np.cos(wc * t + thc) + 0.5 * np.cos(wd * t + thd)
                     + np.cos(we * t + the) + np.random.normal()*0.2
        )


#below: exp decay flares
if exp_fall == True:
    for t in t_list:
        #x_list.append(-(a_exp/t_exp)*(t-int(t/t_exp)*t_exp)+0+0.5*a_exp+60.)
        if t>t_2:
            x_list.append(mean_rate +
                          a_exp * np.exp(-(t-t_0)/t_exp) +
                          a_exp1 * np.exp(-(t-t_1)/t_exp1) +
                          a_exp2 * np.exp(-(t-t_2)/t_exp2))
        elif t>t_1:
            x_list.append(mean_rate + a_exp * np.exp(-(t-t_0)/t_exp) + a_exp1 * np.exp(-(t-t_1)/t_exp1))
        elif t>t_0:
            x_list.append(mean_rate + a_exp * np.exp(-(t-t_0)/t_exp))
        else:
            x_list.append(mean_rate)

#below: symmetric exp flares
if exp_sym == True:
    for t in t_list:
        #x_list.append(-(a_exp/t_exp)*(t-int(t/t_exp)*t_exp)+0+0.5*a_exp+60.)
        if t>t_2:
            x_list.append(mean_rate +
                          a_exp * np.exp(-abs(t-t_0)/t_exp) +
                          a_exp1 * np.exp(-abs(t-t_1)/t_exp1) +
                          a_exp2 * np.exp(-abs(t-t_2)/t_exp2))
        elif t>t_1 and t<=t_2:
            x_list.append(mean_rate +
                          a_exp * np.exp(-abs(t-t_0)/t_exp) +
                          a_exp1 * np.exp(-abs(t-t_1)/t_exp1) +
                          a_exp2 * np.exp(-abs(t-t_2)/t_exp2))
        elif t>t_0 and t<=t_1:
            x_list.append(mean_rate + a_exp * np.exp(-abs(t-t_0)/t_exp) + a_exp1 * np.exp(-abs(t-t_1)/t_exp1))
        elif t<=t_0:
            x_list.append(mean_rate + a_exp * np.exp(-abs(t-t_0)/t_exp))
        else:
            x_list.append(mean_rate)

#below: exp rise
if exp_rise == True:
    for t in t_list:
        #x_list.append(-(a_exp/t_exp)*(t-int(t/t_exp)*t_exp)+0+0.5*a_exp+60.)
        if t>t_1 and t<=t_2:
            x_list.append(mean_rate +
                          a_exp2 * np.exp((t-t_2)/t_exp2)
            )
        elif t>t_0 and t<=t_1:
            x_list.append(mean_rate +
                          a_exp1 * np.exp(-abs(t-t_1)/t_exp1)+
                          a_exp2 * np.exp((t-t_2)/t_exp2)
            )
        elif t<=t_0:
            x_list.append(mean_rate + a_exp * np.exp((t-t_0)/t_exp)+
                          a_exp2 * np.exp((t-t_2)/t_exp) +
                          a_exp1 * np.exp((t-t_1)/t_exp1)
            )
        else:
            x_list.append(mean_rate)

slc = lc.lc(t_list, x_list)
#file = 'M4_1-day-bin_powerlaw2_LC_trimCol.txt'
file = 'M4_1-week-bin_powerlaw2_LC_trimCol.txt'
t,f,df = read_tfdf(file)

t = (t-t[0])*86400.

print t, f
print len(t), len(f)

#slc = lc.lc(t, f)
#slc.interpolate()

#print "check interp",np.diff(slc.t)
#difflist = np.diff(slc.t)
#for i,dt in enumerate(difflist):
#    if dt != 86400.:
#        print i, dt
#
#print np.count_nonzero(np.diff(slc.t) - 86400.)
#
#t = slc.t
#t_all = np.arange(t[0], t[-1] + (t[1]-t[0]), (t[1]- t[0]))
#t_set = set(t_all)
#t_missing = sorted( t_set.difference(t))
#print t_missing



ax = fig.add_subplot(221)
#ax.plot(t_sim, x_sim-np.mean(x_sim),'r-')
#ax.plot(t_sim, x_sim, 'b-')
ax.plot(slc.t, slc.x, 'b-')

Poisson_lc = lc.lc(np.array(t_sim), np.array(x_sim))

#sub_t, sub_x = Poisson_lc.segmentation(n_seg=8)
sub_t, sub_x = slc.segmentation(n_seg=8)

ax.legend(prop=dict(size=12))
ax.set_xlabel('t (s)')
#ax.set_ylabel('rate')
ax.set_ylabel('flux (cm$^{-2}$ s$^{-1}$)')

#freq, sp = FT_continuous(slc.t, slc.x, method=2)

#sp = np.fft.fftshift(np.fft.fft(slc.x))
#freq = np.fft.fftshift(np.fft.fftfreq(slc.t.shape[-1]))
#freq = freq*len(slc.t)/(float(slc.t[-1])-float(slc.t[0]))
#freq = freq[len(freq)/2:]\
#freq, bis, bic, biph = bispectrum.bispec(slc.t, slc.x, bicoh=True, biphase=True)


assert np.array(sub_t).shape == np.array(sub_x).shape
print np.array(sub_t).shape, np.array(sub_x).shape


seg_freq, seg_bisp, seg_bic, seg_biph, seg_biph2 = bispectrum.ensemble_bispec(slc.t,slc.x,len_seg=len_seg, n_seg=n_seg)

if do_rebin == True:
    seg_bisp = rebin_2d_array.rebin(seg_bisp, scale=rebin_scale)
    seg_bic = rebin_2d_array.rebin(seg_bic, scale=rebin_scale)
    seg_biph2 = rebin_2d_array.rebin(seg_biph2, scale=rebin_scale)


print "RRRRResults:",seg_freq, seg_bisp, seg_bic, seg_biph, seg_biph2
print "RRRRResults len:",(seg_freq.shape), (seg_bisp.shape), (seg_bic.shape), (seg_biph.shape), (seg_biph2.shape)
#freq_sim, bis_sim, bic_sim, biph_sim = bispectrum.bispec(t_sim, x_sim, bicoh=True, biphase=True)


#ax2 = fig.add_subplot(312)
#H, xedges,yedges = np.histogram2d(freq,freq,bins=10, weights=bic)
#ax2.imshow(H, extent=[freq[0], freq[-1],freq[0], freq[-1]],
#           interpolation='nearest', origin='lower',aspect=1)

ax2 = fig.add_subplot(222, projection='3d')

#X, Y = np.meshgrid(freq, freq)
#p = ax2.plot_surface(X, Y, bic, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#Xs, Ys = np.meshgrid(freq_sim, freq_sim)
Xseg, Yseg = np.meshgrid(seg_freq, seg_freq)

if do_rebin == True:
    Xseg = rebin_2d_array.rebin(Xseg, scale=rebin_scale)
    Yseg = rebin_2d_array.rebin(Yseg, scale=rebin_scale)

#p = ax2.plot_surface(Xs, Ys, bic_sim, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)

#p = ax2.plot_surface(np.log(Xseg), np.log(Yseg), seg_bisp.imag, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#cset = ax2.contourf(np.log(Xseg), np.log(Yseg), seg_bisp.imag, zdir='z', offset=np.min(seg_bisp.real), cmap=plt.cm.jet)
#ax2.set_xlabel('log(f1) (Hz)')
#ax2.set_ylabel('log(f2) (Hz)')

p = ax2.plot_surface(Xseg, Yseg, seg_bisp, rstride=8, cstride=8, alpha=0.5,
#p = ax2.plot_surface(Xseg, Yseg, seg_bisp, rstride=16, cstride=16, alpha=0.5,
                     cmap=plt.cm.jet, linewidth=0, antialiased=False)

cb = fig.colorbar(p, shrink=0.9)
cset = ax2.contourf(Xseg, Yseg, seg_bisp, zdir='z', offset=np.min(seg_bisp.real), cmap=plt.cm.jet)

xset = ax2.contourf(Xseg, Yseg, seg_bisp.real, zdir='x', offset=np.min(Xseg), cmap=plt.cm.jet)
yset = ax2.contourf(Xseg, Yseg, seg_bisp.real, zdir='y', offset=np.max(Yseg), cmap=plt.cm.jet)

ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax2.legend(prop=dict(size=12))
ax2.set_xlabel('f1 (Hz)')
ax2.set_ylabel('f2 (Hz)')
ax2.set_zlabel('bispectrum')
#ax2.set_zlim(0.8,1.0)
#ax2.set_xlim(1e-4,1.0)
#ax2.set_ylim(1e-4,1.0)

#ax2.xaxis.set_scale('log')
#ax2.yaxis.set_scale('log')

ax3 = fig.add_subplot(223, projection='3d')

#fac=2
#bic_reb = rebin(seg_bic, scale=fac)
#biph2_reb = rebin(seg_biph2, scale=fac)
#Xreb, Yreb = np.meshgrid(seg_freq[0::fac], seg_freq[0::fac])

#p = ax3.plot_surface(X, Y, biph, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#p = ax3.plot_surface(Xs, Ys, biph_sim, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#p = ax3.scatter(Xseg, Yseg, seg_bic, #rstride=8, cstride=8, alpha=0.5,
p = ax3.plot_surface(Xseg, Yseg, seg_bic, rstride=8, cstride=8, alpha=0.5,
                     cmap=plt.cm.jet, linewidth=0, antialiased=False)

#p = ax3.plot_surface(Xreb, Yreb, bic_reb, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)

cb = fig.colorbar(p, shrink=0.9)

#cset = ax3.contourf(Xseg, Yseg, seg_bic, zdir='z', offset=np.min(seg_bic), cmap=plt.cm.jet)
cset = ax3.contourf(Xseg, Yseg, seg_bic, zdir='z', offset=-0.2, cmap=plt.cm.jet)
xset = ax3.contourf(Xseg, Yseg, seg_bic, zdir='x', offset=np.min(Xseg), cmap=plt.cm.jet)
yset = ax3.contourf(Xseg, Yseg, seg_bic, zdir='y', offset=np.max(Yseg), cmap=plt.cm.jet)

ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax3.legend(prop=dict(size=12))
ax3.set_xlabel('f1 (Hz)')
ax3.set_ylabel('f2 (Hz)')
ax3.set_zlim(-0.2,1.02)
ax3.set_zlabel('bicoherence')
#ax3.set_xscale('log')
#ax3.set_yscale('log')

print "max bicoh", np.max(seg_bic)
print "min bicoh", np.min(seg_bic)

ax4 = fig.add_subplot(224, projection='3d')

#p = ax3.plot_surface(X, Y, biph, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#p = ax3.plot_surface(Xs, Ys, biph_sim, rstride=8, cstride=8, alpha=0.5,
#                     cmap=plt.cm.jet, linewidth=0, antialiased=False)
#seg_biph2 = seg_biph

p = ax4.plot_surface(Xseg, Yseg, seg_biph2, rstride=8, cstride=8, alpha=0.5,
#p = ax4.plot_surface(Xreb, Yreb, biph2_reb, rstride=8, cstride=8, alpha=0.5,
                     cmap=plt.cm.jet, linewidth=0, antialiased=False)


cb = fig.colorbar(p, shrink=0.9)
cset = ax4.contourf(Xseg, Yseg, seg_biph2, zdir='z', offset=np.min(seg_biph2), cmap=plt.cm.jet)
xset = ax4.contourf(Xseg, Yseg, seg_biph2, zdir='x', offset=np.min(Xseg), cmap=plt.cm.jet)
yset = ax4.contourf(Xseg, Yseg, seg_biph2, zdir='y', offset=np.max(Yseg), cmap=plt.cm.jet)

ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax4.legend(prop=dict(size=12))
ax4.set_xlabel('f1 (Hz)')
ax4.set_ylabel('f2 (Hz)')
ax4.set_zlabel('biphase')



fig.tight_layout()
#plt.savefig('exp_rise_bispectrum.eps',format='eps', dpi=1000)
#plt.savefig('sawtooth_decay_bispectrum.eps',format='eps', dpi=1000)
#plt.savefig('exp_rise_and_decay_bispectrum.eps',format='eps', dpi=1000)

if save_plot==False:
    plt.show()
else:
    #plt.savefig('bispectrum_multi_sinewave'+str(n_seg)+'segments_'+str(len_seg)+'len_seg.eps',format='eps', dpi=1000)
    #plt.savefig('bispectrum_FermiLAT_daily_'+str(n_seg)+'segments_'+str(len_seg)+'len_seg_rebin'+str(rebin_scale)+'.eps',
    #plt.savefig('bispectrum_FermiLAT_weekly_'+str(n_seg)+'segments_'+str(len_seg)+'len_seg_rebin'+str(rebin_scale)+'.eps',
    #            format='eps', dpi=1000)
    plt.savefig('bispectrum_exp_decay'+str(n_seg)+'segments_'+str(len_seg)+'len_seg.eps',format='eps', dpi=1000)
