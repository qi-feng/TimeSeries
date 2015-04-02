__author__ = 'qfeng'
# This function processes a VEGAS stage 6 file
# and produces light curves in the choice of format
# available formats: histogram, kernel density estimation, bayesian blocks

from matplotlib import pyplot as plt
import numpy as np
from astroML.plotting import hist
from sklearn.neighbors import KernelDensity
from scipy.interpolate import interp1d
from astroML.density_estimation.bayesian_blocks import bayesian_blocks

import ROOT

def norm_ea_kde(t, kde, t_evt, weight_evt, live_t_factor, bws):
    # make histogram of effective area and interpolate in between
    """
    :param t: MJD time list for plotting kde, ** should be longer than t_evt!! **
    :param kde: un-normed kde list, ** same len as t **
    :param t_evt: actual list of event arrival times in MJD
    :param weight_evt: actual list of weight for each evt (1/effective area for on, and alpha/effective area for off)
    :param live_t_factor: LT/dur
    :param bws: bandwidth in second, actually used for smoothing EA
    :return:
    """
    assert len(t_evt) == len(weight_evt)
    assert len(t) == len(kde)
    norm_kde = np.zeros(len(kde))
    bw_hist = bws/86400.
    tt_hist = np.arange(t[0], t[-1] + bw_hist, bw_hist)
    wt_hist, t_hist = np.histogram(t_evt, tt_hist, weights=weight_evt, normed=False)
    #t_evt = np.append(t[0],t_evt)
    #t_evt = np.append(t_evt,t[-1])
    #weight_evt = np.append(weight_evt[0], weight_evt)
    #weight_evt = np.append(weight_evt, weight_evt[-1])
    #f_interp = interp1d(t_evt, weight_evt, kind='cubic')
    t_hist = np.append(t[0],t_hist+0.5*np.mean(np.diff(t)))
    #t_hist = np.append(t_hist,t[-1])
    wt_hist = np.append(wt_hist[0], wt_hist)
    wt_hist = np.append(wt_hist, wt_hist[-1])
    f_interp = interp1d(t_hist, wt_hist, kind='cubic')
    print "interp mean weight:", np.mean(f_interp(t)),np.mean(wt_hist), "sum", wt_hist.sum()
    print "original mean weight:", np.mean(weight_evt), "sum", np.sum(weight_evt)
    print "factor of", np.mean(weight_evt)/np.mean(f_interp(t)),"off"
    subplot = 111
    fig = plt.figure(figsize=(8, 5))
    fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
    ax = fig.add_subplot(subplot)
    ax.plot((t-t[0])*1440.,1/f_interp(t), label = 'interpolated effective area')
    ax.set_xlabel('t (minutes since the start of observations)')
    ax.set_ylabel('interpolated effective area (m$^{-2}$)')
    fig.tight_layout()
    ax.legend(loc='best')
    ax.set_ylim([0,20000])


    norm_kde = f_interp(t)*kde * len(t_evt) / (live_t_factor * 86400.)*np.mean(weight_evt)/np.mean(f_interp(t))
    return norm_kde

def norm_spec_kde(ene, kde, ene_evt, weight_evt, bws):
    """
    :param t: MJD time list for plotting kde, ** should be longer than t_evt!! **
    :param kde: un-normed kde list, ** same len as t **
    :param t_evt: actual list of event arrival times in MJD
    :param weight_evt: actual list of weight for each evt (1/effective area for on, and alpha/effective area for off)
    :param live_t_factor: LT/dur
    :param bws: bandwidth in second, actually used for smoothing EA
    :return:
    """
    assert len(ene_evt) == len(weight_evt)
    assert len(ene) == len(kde)
    norm_kde = np.zeros(len(kde))
    bw_hist = bws
    ene_edges= np.arange(ene[0], ene[-1] + bw_hist, bw_hist)
    wt_hist, ene_hist = np.histogram(ene_evt, ene_edges, weights=weight_evt, normed=False)
    ene_hist = np.append(ene[0],ene_hist+0.5*np.mean(np.diff(ene)))
    #t_hist = np.append(t_hist,t[-1])
    wt_hist = np.append(wt_hist[0], wt_hist)
    wt_hist = np.append(wt_hist, wt_hist[-1])
    f_interp = interp1d(ene_hist, wt_hist, kind='cubic')
    print "interp mean weight:", np.mean(f_interp(ene)),np.mean(wt_hist), "sum", wt_hist.sum()
    print "original mean weight:", np.mean(weight_evt), "sum", np.sum(weight_evt)
    print "factor of", np.mean(weight_evt)/np.mean(f_interp(ene)),"off"
    subplot = 111
    fig = plt.figure(figsize=(8, 5))
    fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
    ax = fig.add_subplot(subplot)
    ax.plot(ene,1/f_interp(ene), label = 'interpolated effective area')
    ax.set_xlabel('energy (GeV)')
    ax.set_ylabel('interpolated effective area (m$^{-2}$)')
    fig.tight_layout()
    ax.legend(loc='best')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_ylim([0,20000])


    norm_kde = f_interp(ene)*kde * len(ene_evt) / np.mean(weight_evt)/np.mean(f_interp(ene))
    return norm_kde

class run_data:

        def __init__(self, run_num):
            self.RunNumber = run_num
            self.RunStartMJD = 0.0
            self.RunStopMJD = 0.0
            self.Alpha = 0.0
            self.LiveTime = 0.0
            self.LiveTimeFactor = 0.0
            self.MinSafeEnergyGeV = 0.0
            self.MaxSafeEnergyGeV = 0.0
            self.tonList = []
            self.toffList = []
            #wt lists are 1/effective area
            self.wtonList = []
            self.wtoffList = []
            self.eonList = []
            self.eoffList = []
            self.ewtonList = []
            self.ewtoffList = []
            self.hasRST = False

        def get_rst(self, s6_filename):

            f = ROOT.TFile(s6_filename, "read")
            tRST = f.Get('RunStatsTree')
            ROOT.gROOT.ProcessLine(
                "struct run_t {\
               Double_t        fAlpha;\
               Int_t           faRunNumber;\
               Double_t        faRunStartMJD;\
               Double_t        faRunEndMJD;\
               Double_t        faLiveTime;\
               Double_t        fMinSafeEnergyGeV;\
               Double_t        fMaxSafeEnergyGeV;\
            }")

            rst = ROOT.run_t()

            tRST.SetBranchAddress("faRunNumber", ROOT.AddressOf(rst, 'faRunNumber'))
            tRST.SetBranchAddress("faRunStartMJD", ROOT.AddressOf(rst, 'faRunStartMJD'))
            tRST.SetBranchAddress("faRunEndMJD", ROOT.AddressOf(rst, 'faRunEndMJD'))
            tRST.SetBranchAddress("fAlpha", ROOT.AddressOf(rst, 'fAlpha'))
            tRST.SetBranchAddress("faLiveTime", ROOT.AddressOf(rst, 'faLiveTime'))
            tRST.SetBranchAddress("fMinSafeEnergyGeV", ROOT.AddressOf(rst, 'fMinSafeEnergyGeV'))
            tRST.SetBranchAddress("fMaxSafeEnergyGeV", ROOT.AddressOf(rst, 'fMaxSafeEnergyGeV'))

            totalLT = 0.0
            alpha = 0.0
            minSafeE = 100000.0
            maxSafeE = -10000.0

            found_run = False
            for i in range(0, tRST.GetEntries()):
                tRST.GetEntry(i)
                print rst.faRunNumber
                if self.RunNumber == rst.faRunNumber:
                    self.RunStartMJD = rst.faRunStartMJD
                    self.RunStopMJD = rst.faRunEndMJD
                    self.Alpha = rst.fAlpha
                    self.LiveTime = rst.faLiveTime
                    self.MinSafeEnergyGeV = rst.fMinSafeEnergyGeV
                    self.MaxSafeEnergyGeV = rst.fMaxSafeEnergyGeV

                    print "Run number:", rst.faRunNumber, "Alpha:", rst.fAlpha, "Live time:", rst.faLiveTime, \
                        "fMinSafeEnergyGeV", rst.fMinSafeEnergyGeV, "fMaxSafeEnergyGeV", rst.fMaxSafeEnergyGeV
                    found_run = True
                    break

            if found_run == False:
                print "Cannot find run stats tree for this run", self.RunNumber

            self.LiveTimeFactor = self.LiveTime / ((self.RunStopMJD - self.RunStartMJD) * 86400. )
            self.hasRST = True

        def get_est(self, s6_filename, energy_GeV = 200., energy_hi_GeV= 30000., spectrum=False):

            f = ROOT.TFile(s6_filename, "read")
            tEST = f.Get('EventStatsTree')

            ROOT.gROOT.ProcessLine(
                "struct evt_t {\
               Double_t        MJDDbl;\
               Float_t         EnergyGeV;\
               Float_t         EffectiveArea;\
               Bool_t          OnEvent;\
               Bool_t          OffEvent;\
               Float_t         Alpha;\
               UInt_t          RunNum;\
            }")

            evt = ROOT.evt_t()

            tEST.SetBranchAddress("EnergyGeV", ROOT.AddressOf(evt, 'EnergyGeV'))
            tEST.SetBranchAddress("MJDDbl", ROOT.AddressOf(evt, 'MJDDbl'))
            tEST.SetBranchAddress("EffectiveArea", ROOT.AddressOf(evt, 'EffectiveArea'))
            tEST.SetBranchAddress("OnEvent", ROOT.AddressOf(evt, 'OnEvent'))
            tEST.SetBranchAddress("OffEvent", ROOT.AddressOf(evt, 'OffEvent'))

            if self.hasRST == False:
                self.get_rst(s6_filename)

            on_count = 0
            off_count = 0
            skip_on_count = 0
            skip_off_count = 0

            for i in range(1, tEST.GetEntries()):
                tEST.GetEntry(i)
                if evt.OnEvent and evt.EnergyGeV >= energy_GeV and energy_hi_GeV < energy_hi_GeV \
                        and evt.MJDDbl >= self.RunStartMJD and evt.MJDDbl <= self.RunStopMJD:
                    on_count += 1
                    if evt.EffectiveArea <= 0:
                        print "ERROR: found an on event with a negative or null effective area",\
                        evt.EffectiveArea, ", skipping event"
                        skip_on_count += 1
                        continue
                    else:
                        self.wtonList.append(1. / float(evt.EffectiveArea))
                        self.tonList.append(evt.MJDDbl)

                elif evt.OffEvent and evt.EnergyGeV >= energy_GeV and energy_hi_GeV < energy_hi_GeV and \
                                evt.MJDDbl >= self.RunStartMJD and evt.MJDDbl <= self.RunStopMJD:
                    #f_off_evt.append(1. / float(evt.EffectiveArea) * alpha)
                    off_count += 1
                    if evt.EffectiveArea <= 0:
                        print "ERROR: found an off event with a negative or null effective area",\
                        evt.EffectiveArea, ", skipping event"
                        skip_off_count += 1
                        continue
                    else:
                        self.wtoffList.append(1. / float(evt.EffectiveArea))
                        self.toffList.append(evt.MJDDbl)
                if spectrum == True:
                    if evt.OnEvent and evt.EffectiveArea > 0.1:
                        self.eonList.append(evt.EnergyGeV)
                        self.ewtonList.append(1. / float(evt.EffectiveArea))
                    elif evt.OffEvent and evt.EffectiveArea > 0.1:
                        self.eoffList.append(evt.EnergyGeV)
                        self.ewtoffList.append(1. / float(evt.EffectiveArea))

            if spectrum==True:
                #print "Getting spectrum"
                self.eonList, self.ewtonList = zip(*sorted(zip(self.eonList, self.ewtonList)))
                self.eoffList, self.ewtoffList = zip(*sorted(zip(self.eoffList, self.ewtoffList)))
                self.eonList = np.array(self.eonList)
                self.eoffList = np.array(self.eoffList)
                self.ewtonList = np.array(self.ewtonList)
                self.ewtoffList = np.array(self.ewtoffList)


            if skip_off_count + skip_on_count > 0:
                print "In file", s6_filename
                print skip_on_count, "out of", on_count, "on events skipped"
                print "ratio of", float(skip_on_count)/float(on_count)
                print skip_off_count, "out of", off_count, "off events skipped"
                print "ratio of", float(skip_off_count)/float(off_count)

            #t_evt = np.array(t_evt)[:, np.newaxis]
            #t_off_evt = np.array(t_off_evt)[:, np.newaxis]

def get_run_list(s6_filename, verbose=False):
    import ROOT

    f = ROOT.TFile(s6_filename, "read")
    tRST = f.Get('RunStatsTree')
    ROOT.gROOT.ProcessLine(
        "struct run_num {\
       Int_t           faRunNumber;\
       Double_t        faRunStartMJD;\
       Double_t        faRunEndMJD;\
       Double_t        faLiveTime;\
       Double_t        fMinSafeEnergyGeV;\
       Double_t        fMaxSafeEnergyGeV;\
    }")

    rst = ROOT.run_num()

    tRST.SetBranchAddress("faRunNumber", ROOT.AddressOf(rst, 'faRunNumber'))
    tRST.SetBranchAddress("faRunStartMJD", ROOT.AddressOf(rst, 'faRunStartMJD'))
    tRST.SetBranchAddress("faRunEndMJD", ROOT.AddressOf(rst, 'faRunEndMJD'))
    tRST.SetBranchAddress("faLiveTime", ROOT.AddressOf(rst, 'faLiveTime'))
    tRST.SetBranchAddress("fMinSafeEnergyGeV", ROOT.AddressOf(rst, 'fMinSafeEnergyGeV'))
    tRST.SetBranchAddress("fMaxSafeEnergyGeV", ROOT.AddressOf(rst, 'fMaxSafeEnergyGeV'))

    totalLT = 0.0
    minSafeE = 100000.0
    maxSafeE = -10000.0

    run_num_list = []
    for i in range(0, tRST.GetEntries()):
        tRST.GetEntry(i)
        totalLT += rst.faLiveTime
        if minSafeE > rst.fMinSafeEnergyGeV:
            minSafeE = rst.fMinSafeEnergyGeV
        if maxSafeE < rst.fMaxSafeEnergyGeV:
            maxSafeE = rst.fMaxSafeEnergyGeV

        run_num_list.append(rst.faRunNumber)
        #print "Run number:", rst.faRunNumber, "Alpha:", rst.fAlpha, "Live time:", rst.faLiveTime, \
        #    "fMinSafeEnergyGeV", rst.fMinSafeEnergyGeV, "fMaxSafeEnergyGeV", rst.fMaxSafeEnergyGeV

        if verbose == True:
            print "This run:", rst.faRunNumber
            print "Live time:", rst.faLiveTime, \
            "Duration:", (rst.faRunEndMJD - rst.faRunStartMJD)*86400.,"s"
            #"fMinSafeEnergyGeV", rst.fMinSafeEnergyGeV, "fMaxSafeEnergyGeV", rst.fMaxSafeEnergyGeV

    return run_num_list

def plot_list(t, kde, t_hist, f_hist, w_hist, t_start = None, t_stop = None, ax = None,
              use_mjd = True, plotname=None, compfile = None, plot_hist = True,
              kde_label = None, hist_label = None, bb=None, t_bb=None, h_bb=None, bb_label="Bayesian blocks"):

    if t_start == None:
        t_start = t[0]
    if t_stop == None:
        t_stop = t[-1]

    doplot = False
    MJD_flag = use_mjd
    if ax == None:
        subplot = 111
        fig = plt.figure(figsize=(8, 5))
        fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
        ax = fig.add_subplot(subplot)
        doplot = True
    # plot TTE
    ##ax.plot(t_evt, min(min(norm_dens_kde),min(norm_dens_kde_epanechnikov)) * 0.992 * np.ones(len(t_evt)), '|k')

    #c_bb, t_bb, p_bb = hist((t_evt[:,0]-t_start)*1440., bins='blocks', ax=ax, color='black',histtype='step', normed=True)
    #ax.cla()
    #weight_seq = c_bb.sum()/(np.diff(t_bb)*live_t_factor)

    #ax.bar(t_bb[:-1], c_bb*weight_seq, width=np.diff(t_bb), alpha=0.6, color='green')

    if kde_label == None:
        kde_label = " Kernel Density Estimation"
    if MJD_flag == False:
        ax.plot((t - t_start) * 1440., kde, '-', color='red', zorder=3,
                label=kde_label)
        #ax.plot((t - t_start) * 1440., norm_off_dens_kde, '--', color='red',
        #        zorder=3)  #, label="Off events Epanechnikov Kernel Density (h="+str(bws)+"s)")
    else:
        ax.plot((t), kde, '-', color='red', zorder=3,
                label=kde_label)
        #ax.plot((t), norm_off_dens_kde, '--', color='red',
        #        zorder=3)  #, label="Off events Epanechnikov Kernel Density (h="+str(bws)+"s)")

    if hist_label == None:
        hist_label = "Flux histogram"
    if plot_hist==True:
        if MJD_flag == True:
            ax.bar(t_hist[:-1], f_hist, width=w_hist, alpha=0.5,
                   color='gray', label=hist_label)
        else:
            ax.bar((t_hist[:-1] - t_start) * 1440., f_hist,
                   width=np.array(w_hist) * 1440., alpha=0.5, color='gray', label=hist_label)

    if bb==True:
        if MJD_flag == True:
            ax.bar(t_bb[:-1], h_bb, width=np.diff(t_bb), alpha=0.5,
                   color='green', label="Bayesian blocks")
        else:
            ax.bar((t_bb[:-1] - t_start) * 1440., h_bb,
                   width=np.array(np.diff(t_bb)) * 1440., alpha=0.5, color='green', label="Bayesian blocks")

    if compfile != None:
        t_min = []
        f_min = []
        df_min = []
        with open(compfile) as fp:
            for line in fp:
                t_min.append(float(line.split()[0]))
                f_min.append(float(line.split()[1]))
                df_min.append(float(line.split()[2]))
        if MJD_flag == True:
            ax.errorbar(np.array(t_min), np.array(f_min), xerr=np.mean(np.diff(t_min))/2., yerr=np.array(df_min),
                        color='blue', fmt='o',label="Rigorous flux histogram")
            #            color='blue', fmt='o',label="Official VERITAS")
            #ax.bar(t2_hist[:-1], (c2_hist * hist_norm - c2_off_hist * hist_off_norm), width=np.diff(t_bins), alpha=0.6,
            #   color='gray', label="Flux histogram")
        else:
            ax.errorbar((np.array(t_min) - t_start) * 1440., np.array(f_min), xerr=np.mean(np.diff(t_min))/2.* 1440., yerr=np.array(df_min),
                        color='blue', fmt='o',label="Rigorous flux histogram")

    if doplot == True:
        if MJD_flag == True:
            ax.set_xlabel('t (MJD)')
        else:
            ax.set_xlabel('t (minutes since the start of observations)')

        ax.set_ylabel('flux (m$^{-2}$ s$^{-1}$)')
        ax.legend(loc='best')

        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        fig.tight_layout()

        if plotname == None:
            plt.show()
        else:
            fig.savefig(plotname,format='eps', dpi=1000)

def process_one_run(s6_filename, run_num, energy_GeV = 200., binwidth_min = 20., kernel = 'epanechnikov',
                    kernel_bandwidth_min = 1., use_mjd = True, fixstart=None, plotname=None, compfile = None,
                    plot_hist = True, spectrum=False, ene_bw=50., energy_hi_GeV = 30000.,
                    ea_method = 'avg', doplot = False, ax = None, bb=True, p0=0.05):

    MJD_flag = use_mjd

    rd = run_data(run_num)
    rd.get_rst(s6_filename)
    rd.get_est(s6_filename, energy_GeV = energy_GeV, energy_hi_GeV= energy_hi_GeV, spectrum=spectrum)
    t_evt = np.array(rd.tonList)[:, np.newaxis]
    t_off_evt = np.array(rd.toffList)[:, np.newaxis]

    bws = 60.* kernel_bandwidth_min
    bw = bws / 86400.
    # t for plotting KDE
    t = np.arange(rd.RunStartMJD, rd.RunStopMJD + bw, bw / 10.)[:, np.newaxis]

    kde = KernelDensity(bw, kernel=kernel).fit(t_evt)
    dens_kde = np.exp(kde.score_samples(t))

    off_kde = KernelDensity(bw, kernel=kernel).fit(t_off_evt)
    off_dens_kde = np.exp(off_kde.score_samples(t))

    if spectrum == True:
        print rd.eonList
        #ene = np.arange(rd.MinSafeEnergyGeV, rd.MaxSafeEnergyGeV + ene_bw, ene_bw / 10.)[:, np.newaxis]
        ene = np.arange(100., 10000. + ene_bw, ene_bw / 10.)[:, np.newaxis]
        ene_off = ene
        spec_kde = KernelDensity(ene_bw, kernel=kernel).fit(np.array(rd.eonList)[:, np.newaxis])
        spec_kde_dens = np.exp(spec_kde.score_samples(ene))
        norm_spec_kde_dens = norm_spec_kde(ene[:,0], spec_kde_dens, np.array(rd.eonList), np.array(rd.ewtonList), 100.)
        norm_spec_kde_dens = norm_spec_kde_dens/rd.LiveTime
        #ene_off = np.arange(rd.MinSafeEnergyGeV, rd.MaxSafeEnergyGeV + ene_bw, ene_bw / 10.)[:, np.newaxis]
        spec_kde_off = KernelDensity(ene_bw, kernel=kernel).fit(np.array(rd.eoffList)[:, np.newaxis])
        spec_kde_dens_off = np.exp(spec_kde_off.score_samples(ene_off))
        norm_spec_kde_dens_off = norm_spec_kde(ene_off[:,0], spec_kde_dens_off, np.array(rd.eoffList), np.array(rd.ewtoffList), 100.)
        norm_spec_kde_dens_off = norm_spec_kde_dens_off * rd.Alpha / rd.LiveTime
        fig1 = plt.figure(figsize=(10, 8))
        fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
        ax1= fig1.add_subplot(111)
        ax1.plot(ene, norm_spec_kde_dens,'r-')
        ax1.plot(ene_off, norm_spec_kde_dens_off,'b--')
        ax1.plot(ene, norm_spec_kde_dens - norm_spec_kde_dens_off,'g-')
        ax1.set_xscale('log')
        ax1.set_yscale('log')



    if ea_method == 'avg':
        norm_on = bw * len(t_evt) / (rd.LiveTimeFactor * bws) * np.mean(rd.wtonList)
        norm_off = bw * len(t_off_evt) / (rd.LiveTimeFactor * bws) * np.mean(rd.wtoffList) * rd.Alpha
        norm_dens_kde = dens_kde * norm_on
        norm_off_dens_kde = off_dens_kde * norm_off
    else:
        norm_dens_kde = norm_ea_kde(t[:,0], dens_kde,t_evt[:,0],rd.wtonList, rd.LiveTimeFactor, bws*2)
        norm_off_dens_kde = norm_ea_kde(t[:,0], off_dens_kde,t_off_evt[:,0],rd.wtoffList, rd.LiveTimeFactor, bws*2)*rd.Alpha

    print "check", norm_dens_kde.sum()
    print "check off KDE", norm_off_dens_kde.sum(), norm_off_dens_kde.sum() / np.mean(rd.wtonList)

    bws_hist = binwidth_min * 60.
    bw_hist = bws_hist / 86400.
    # generate histogram bin edges
    t_bins = np.arange(rd.RunStartMJD, rd.RunStopMJD + bw_hist, bw_hist)

    #ea_hist, t_hist = np.histogram(t_evt[:,0], t_bins, weights=f_evt, normed=False)
    #ax.bar(t_hist[:-1], c_off_hist*alpha/(live_t_factor*120.), width=np.diff(t_bins), alpha=0.6, color='gray',label="Off flux")
    c_hist, t_hist = np.histogram(t_evt[:, 0], t_bins, weights=rd.wtonList, normed=False)
    c_off_hist, t_off_hist = np.histogram(t_off_evt[:, 0], t_bins, weights=rd.wtoffList, normed=False)

    c2_hist, t2_hist = np.histogram(t_evt[:, 0], t_bins, weights=rd.wtonList, density=True)
    c2_off_hist, t2_off_hist = np.histogram(t_off_evt[:, 0], t_bins, weights=rd.wtoffList, density=True)

    hist_norm = bw_hist * c_hist.sum() / (rd.LiveTimeFactor * bws_hist)
    hist_off_norm = bw_hist * c_off_hist.sum() / (rd.LiveTimeFactor * bws_hist) * rd.Alpha

    print "kernel:", norm_dens_kde.sum()
    print "raw sum histo:", c_hist.sum()
    print "length of on event list:", len(t_evt[:, 0])
    print "histo not normed:", c_hist.sum(), "approx count:", c_hist.sum() / np.mean(rd.wtonList), "rate:", c_hist.sum() / (
        rd.LiveTimeFactor * bws)
    print "histo density=True approx count :", c2_hist.sum() * hist_norm * (rd.LiveTimeFactor * bws_hist) / np.mean(rd.wtonList)

    #print "bb normed:", c_bb.sum()
    print "off histo not normed:", c_off_hist.sum(), "approx count:", c_off_hist.sum() / np.mean(
        rd.wtoffList), "rate:", c_hist.sum() / (rd.LiveTimeFactor * bws)
    print c_hist.sum() / (rd.LiveTimeFactor * bws), c2_hist.sum() * hist_norm
    print "off", c2_off_hist.sum() * hist_off_norm
    print "total on events:", len(t_evt)
    print "total off events (before applying alpha):", len(t_off_evt)
    print "total off events (after applying alpha):", len(t_off_evt)*rd.Alpha
    print "histo density=True approx count :", c2_off_hist.sum() * hist_off_norm * (rd.LiveTimeFactor * bws_hist) / np.mean(rd.wtoffList)

    flux_hist = (c2_hist * hist_norm - c2_off_hist * hist_off_norm)

    if MJD_flag == False and fixstart == None:
        t_bb = bayesian_blocks((t_evt[:, 0]-rd.RunStartMJD)*1440., p0=p0)
        c_bb, t_bb = np.histogram((t_evt[:, 0]-rd.RunStartMJD)*1440., t_bb, normed=False)
        weight_seq = 1./(np.diff(t_bb)*60.*rd.LiveTimeFactor)*np.mean(rd.wtonList)
        c_off_bb, t_off_bb = np.histogram((t_off_evt[:, 0]-rd.RunStartMJD)*1440., t_bb, normed=False)
        weight_off_seq = 1./(np.diff(t_bb)*60.*rd.LiveTimeFactor)*np.mean(rd.wtoffList)*rd.Alpha
    elif MJD_flag == False and fixstart != None:
        t_bb = bayesian_blocks((t_evt[:, 0]-fixstart)*1440., p0=p0)
        c_bb, t_bb = np.histogram((t_evt[:, 0]-fixstart)*1440., t_bb, normed=False)
        weight_seq = 1./(np.diff(t_bb)*60.*rd.LiveTimeFactor)*np.mean(rd.wtonList)
        c_off_bb, t_off_bb = np.histogram((t_off_evt[:, 0]-fixstart)*1440., t_bb, normed=False)
        weight_off_seq = 1./(np.diff(t_bb)*60.*rd.LiveTimeFactor)*np.mean(rd.wtoffList)*rd.Alpha
    else:
        t_bb = bayesian_blocks(t_evt[:, 0], p0=p0)
        c_bb, t_bb = np.histogram(t_evt[:, 0], t_bb, normed=False)
        weight_seq = 1./(np.diff(t_bb)*86400.*rd.LiveTimeFactor)*np.mean(rd.wtonList)
        c_off_bb, t_off_bb = np.histogram(t_off_evt[:, 0], t_bb, normed=False)
        weight_off_seq = 1./(np.diff(t_bb)*86400.*rd.LiveTimeFactor)*np.mean(rd.wtoffList)*rd.Alpha
    c_bb = c_bb*weight_seq - c_off_bb*weight_off_seq
    width_list = np.diff(t_bb)

    if doplot == True:
        if bb == True:
            fig = plt.figure(figsize=(10, 8))
            fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
            ax = fig.add_subplot(111)
            ax.bar(t_bb[:-1], c_bb*weight_seq, width=np.diff(t_bb), alpha=0.6, color='green', label = "Bayesian Blocks")
        if fixstart == None:
            plot_list(t[:, 0], norm_dens_kde - norm_off_dens_kde, t2_hist, flux_hist, np.diff(t2_hist), t_start=rd.RunStartMJD, t_stop=rd.RunStopMJD,
                  ax = ax, use_mjd = use_mjd, plotname=plotname, compfile = compfile,
                  plot_hist = plot_hist, kde_label = str(kernel)+" Kernel Density (h=" + str(bws) + "s)", hist_label = None)
        else:
            plot_list(t[:, 0], norm_dens_kde - norm_off_dens_kde, t2_hist, flux_hist, np.diff(t2_hist), t_start=fixstart, t_stop=rd.RunStopMJD,
                  ax = ax, use_mjd = use_mjd, plotname=plotname, compfile = compfile,
                  plot_hist = plot_hist, kde_label = str(kernel)+" Kernel Density (h=" + str(bws) + "s)", hist_label = None)

    if MJD_flag == True:
        return t[:, 0], norm_dens_kde - norm_off_dens_kde, t2_hist, flux_hist, t_bb, c_bb, width_list
    elif MJD_flag == False and fixstart != None:
        return (t[:, 0] - fixstart) * 1440., norm_dens_kde - norm_off_dens_kde, (t2_hist - fixstart) * 1440., flux_hist, (t_bb - fixstart) * 1440., c_bb, width_list
    elif MJD_flag == False and fixstart == None:
        return (t[:, 0] - rd.RunStartMJD) * 1440., norm_dens_kde - norm_off_dens_kde, (t2_hist - rd.RunStartMJD) * 1440., flux_hist, (t_bb - rd.RunStartMJD) * 1440., c_bb, width_list


def process_s6_list(s6_filename, energy_GeV = 200., binwidth_min = 20., kernel = 'epanechnikov', ea_method = 'avg',
                    kernel_bandwidth_min = 1., use_mjd = True, plotname=None, compfile=None, plot_hist=True,
                    doplot=False, bb=True, p0=0.05):
    run_list = get_run_list(s6_filename)
    t_list = []
    c_list = []
    th_list = []
    ch_list = []
    hw_list = []
    tbb_list = []
    cbb_list = []
    bbw_list = []
    if use_mjd == False:
        rd = run_data(run_list[0])
        rd.get_rst(s6_filename)
        fixstart = rd.RunStartMJD
    else:
        fixstart = None

    for run in run_list:
        print "processing run ", run

        t, c , t_h, c_h, t_bb, c_bb, bb_w = process_one_run(s6_filename, run, energy_GeV = energy_GeV,
                    energy_hi_GeV = 30000., binwidth_min = binwidth_min, kernel = kernel,
                    kernel_bandwidth_min = kernel_bandwidth_min, use_mjd = True, fixstart=fixstart, plotname=None, compfile = None, plot_hist = False,
                    ea_method = ea_method, doplot=False, p0=p0)

        t_list = t_list + list(t)
        c_list = c_list + list(c)
        th_list = th_list + list(t_h[:-1])
        ch_list = ch_list + list(c_h)
        hw_list = hw_list + list(np.diff(t_h))
        tbb_list = tbb_list + list(t_bb[:-1])
        cbb_list = cbb_list + list(c_bb)
        bbw_list = bbw_list + list(bb_w)

    th_list.append(t_h[-1])
    tbb_list.append(t_bb[-1])
    print len(t_list), len(c_list)
    print len(th_list), len(ch_list), len(hw_list)
    print len(tbb_list), len(cbb_list), len(bbw_list)

    print "start", fixstart
    rd2 = run_data(run_list[-1])
    rd2.get_rst(s6_filename)
    print "stop", rd2.RunStopMJD

    if doplot ==True:
        if bb == True:
            from astroML.density_estimation.bayesian_blocks import bayesian_blocks
            fig = plt.figure(figsize=(10, 8))
            fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
            ax = fig.add_subplot(111)
            if use_mjd == True:
                ax.bar(tbb_list[:-1], cbb_list, width=bbw_list, alpha=0.6, color='green', label = "Bayesian Blocks")
            else:
                ax.bar((np.array(tbb_list[:-1])-fixstart)*1440., cbb_list, width=np.array(bbw_list)*1440., alpha=0.6, color='green', label = "Bayesian Blocks")

        plot_list(t_list, c_list, th_list, ch_list, hw_list, t_start=None, t_stop=None,
              ax = ax, use_mjd = use_mjd, plotname=plotname, compfile = compfile,
              plot_hist = plot_hist, kde_label = None, hist_label = None)

        if use_mjd == True:
            ax.set_xlabel('t (MJD)')
        else:
            ax.set_xlabel('t (minutes since the start of observations)')

        #ax.set_xlim(fixstart, )
        ax.set_ylabel('flux (m$^{-2}$ s$^{-1}$)')
        ax.legend(loc='best')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_ylim(0,max(c_list)*1.2)

    return t_list, c_list, th_list, ch_list, tbb_list, cbb_list

def process_s6_root(s6_filename, energy_GeV = 200., binwidth_min = 20., kernel = 'epanechnikov',
                    kernel_bandwidth_min = 1., use_mjd = True, plotname=None, compfile=None, plot_hist=True):

    f = ROOT.TFile(s6_filename, "read")
    #f = ROOT.TFile("/Users/qfeng/Data/Blazars/RXJ0648.7+1516/results_UA2014updateRebin5s6.root", "read")
    # fname = 'M4_Std_20100217exp_420GeVto1TeV_2min.txt'

    MJD_flag = use_mjd
    tEST = f.Get('EventStatsTree')

    ROOT.gROOT.ProcessLine(
        "struct evt_t {\
       Double_t        MJDDbl;\
       Float_t         EnergyGeV;\
       Float_t         EffectiveArea;\
       Bool_t          OnEvent;\
       Bool_t          OffEvent;\
       Float_t         Alpha;\
       UInt_t          RunNum;\
    }")

    evt = ROOT.evt_t()

    tRST = f.Get('RunStatsTree')
    ROOT.gROOT.ProcessLine(
        "struct run_t {\
       Double_t        fAlpha;\
       Int_t           faRunNumber;\
       Double_t        faRunStartMJD;\
       Double_t        faRunEndMJD;\
       Double_t        faLiveTime;\
       Double_t        fMinSafeEnergyGeV;\
       Double_t        fMaxSafeEnergyGeV;\
    }")

    rst = ROOT.run_t()

    tRST.SetBranchAddress("faRunNumber", ROOT.AddressOf(rst, 'faRunNumber'))
    tRST.SetBranchAddress("faRunStartMJD", ROOT.AddressOf(rst, 'faRunStartMJD'))
    tRST.SetBranchAddress("faRunEndMJD", ROOT.AddressOf(rst, 'faRunEndMJD'))
    tRST.SetBranchAddress("fAlpha", ROOT.AddressOf(rst, 'fAlpha'))
    tRST.SetBranchAddress("faLiveTime", ROOT.AddressOf(rst, 'faLiveTime'))
    tRST.SetBranchAddress("fMinSafeEnergyGeV", ROOT.AddressOf(rst, 'fMinSafeEnergyGeV'))
    tRST.SetBranchAddress("fMaxSafeEnergyGeV", ROOT.AddressOf(rst, 'fMaxSafeEnergyGeV'))


    totalLT = 0.0
    alpha = 0.0
    minSafeE = 100000.0
    maxSafeE = -10000.0
    for i in range(0, tRST.GetEntries()):
        tRST.GetEntry(i)
        totalLT += rst.faLiveTime
        alpha += rst.fAlpha
        if minSafeE > rst.fMinSafeEnergyGeV:
            minSafeE = rst.fMinSafeEnergyGeV
        if maxSafeE < rst.fMaxSafeEnergyGeV:
            maxSafeE = rst.fMaxSafeEnergyGeV

        print "Run number:", rst.faRunNumber, "Alpha:", rst.fAlpha, "Live time:", rst.faLiveTime, \
            "fMinSafeEnergyGeV", rst.fMinSafeEnergyGeV, "fMaxSafeEnergyGeV", rst.fMaxSafeEnergyGeV

    alpha = alpha / tRST.GetEntries()
    print "mean alpha:", alpha, "total live time:", totalLT

    tEST.SetBranchAddress("EnergyGeV", ROOT.AddressOf(evt, 'EnergyGeV'))
    tEST.SetBranchAddress("MJDDbl", ROOT.AddressOf(evt, 'MJDDbl'))
    tEST.SetBranchAddress("EffectiveArea", ROOT.AddressOf(evt, 'EffectiveArea'))
    tEST.SetBranchAddress("OnEvent", ROOT.AddressOf(evt, 'OnEvent'))
    tEST.SetBranchAddress("OffEvent", ROOT.AddressOf(evt, 'OffEvent'))

    #Get the time of the first event
    tEST.GetEntry(0)
    print evt.MJDDbl, evt.EnergyGeV, evt.OnEvent, evt.OffEvent, evt.EffectiveArea
    t_start = evt.MJDDbl

    print tEST.GetEntries(), "entries found"

    #Get the time of the last event
    tEST.GetEntry(tEST.GetEntries() - 1)
    t_stop = evt.MJDDbl
    live_t_factor = totalLT / ((t_stop - t_start) * 86400. )

    t_evt = []
    f_evt = []
    t_off_evt = []
    f_off_evt = []
    #tf_evt = []

    for i in range(1, tEST.GetEntries()):
        tEST.GetEntry(i)
        if evt.OnEvent and evt.EnergyGeV >= energy_GeV:
            f_evt.append(1. / float(evt.EffectiveArea))
            t_evt.append(evt.MJDDbl)
            #tf_evt.append([evt.MJDDbl, 1.])
        elif evt.OffEvent and evt.EnergyGeV >= energy_GeV:
            #f_off_evt.append(1. / float(evt.EffectiveArea) * alpha)
            f_off_evt.append(1. / float(evt.EffectiveArea))
            t_off_evt.append(evt.MJDDbl)

    print "mean effective area", 1./np.mean(np.array(f_evt))
    print "raw duration:", float((t_evt[-1] - t_evt[0])) * 86400., "s"

    t_evt = np.array(t_evt)[:, np.newaxis]
    t_off_evt = np.array(t_off_evt)[:, np.newaxis]
    #tf_evt = np.array(tf_evt)
    # bandwidth

    bws = 60.* kernel_bandwidth_min
    bw = bws / 86400.
    # t for plotting KDE
    #t = np.arange(t_start, t_stop+bw, bw)[:, np.newaxis]
    t = np.arange(t_start, t_stop + bw, bw / 100.)[:, np.newaxis]

    kde = KernelDensity(kernel=kernel, bandwidth=bw).fit(t_evt)
    #kde_exp = KernelDensity(kernel='exponential', bandwidth=bw).fit(t_evt)
    #kde_epanechnikov = KernelDensity(bw, kernel='epanechnikov').fit(t_evt)
    dens_kde = np.exp(kde.score_samples(t))
    #dens_kde_epanechnikov = np.exp(kde_epanechnikov.score_samples(t))
    #dens_kde_exp = np.exp(kde_exp.score_samples(t))

    off_kde = KernelDensity(bw, kernel=kernel).fit(t_off_evt)
    #off_dens_kde = np.exp(off_kde.score_samples(t))*alpha
    off_dens_kde = np.exp(off_kde.score_samples(t))
    #off_kde_epanechnikov = KernelDensity(bw, kernel='epanechnikov').fit(t_off_evt)
    #off_dens_kde_epanechnikov = np.exp(off_kde_epanechnikov.score_samples(t))*alpha
    #off_dens_kde_epanechnikov = np.exp(off_kde_epanechnikov.score_samples(t))
    #off_kde_exp = KernelDensity(bw, kernel='exponential').fit(t_off_evt)
    #off_dens_kde_exp = np.exp(off_kde_exp.score_samples(t))*alpha
    #off_dens_kde_exp = np.exp(off_kde_exp.score_samples(t))

    #print dens_kde.sum() * bw, dens_kde_epanechnikov.sum() * bw, dens_kde_exp.sum() * bw
    #print off_dens_kde.sum() * bw, off_dens_kde_epanechnikov.sum() * bw, off_dens_kde_exp.sum() * bw

    #norm_on = bw*len(t_evt)/(live_t_factor*(float(t_stop)-float(t_start))*86400.)*np.mean(f_evt)
    #norm_off = bw*len(t_off_evt)/(live_t_factor*(float(t_stop)-float(t_start))*86400.)*np.mean(f_off_evt)
    norm_on = bw * len(t_evt) / (live_t_factor * bws) * np.mean(f_evt)
    norm_off = bw * len(t_off_evt) / (live_t_factor * bws) * np.mean(f_off_evt) * alpha

    norm_dens_kde = dens_kde * norm_on
    #norm_dens_kde_epanechnikov = dens_kde_epanechnikov * norm_on
    #norm_dens_kde_exp = dens_kde_exp * norm_on

    norm_off_dens_kde = off_dens_kde * norm_off
    #norm_off_dens_kde_epanechnikov = off_dens_kde_epanechnikov * norm_off
    #norm_off_dens_kde_exp = off_dens_kde_exp * norm_off

    print "check", norm_dens_kde.sum()
    print "check off KDE", norm_off_dens_kde.sum(), off_dens_kde.sum() / np.mean(f_evt)

    subplot = 111
    fig = plt.figure(figsize=(8, 5))
    fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
    ax = fig.add_subplot(subplot)
    # plot TTE
    ##ax.plot(t_evt, min(min(norm_dens_kde),min(norm_dens_kde_epanechnikov)) * 0.992 * np.ones(len(t_evt)), '|k')

    bws_hist = binwidth_min * 60.
    bw_hist = bws_hist / 86400.
    # generate histogram bin edges
    t_bins = np.arange(t_start, t_stop + bw_hist, bw_hist)

    #c_bb, t_bb, p_bb = hist((t_evt[:,0]-t_start)*1440., bins='blocks', ax=ax, color='black',histtype='step', normed=True)
    #ax.cla()
    #weight_seq = c_bb.sum()/(np.diff(t_bb)*live_t_factor)

    #ax.bar(t_bb[:-1], c_bb*weight_seq, width=np.diff(t_bb), alpha=0.6, color='green')

    if MJD_flag == False:
        #ax.plot(t[:, 0], norm_dens_kde, '-', color='black', zorder=3, label="On events Gaussian Kernel Density (h="+str(bws)+"s)")
        ax.plot((t[:, 0] - t_start) * 1440., norm_dens_kde, '--', color='red', zorder=3,
                label="On/Off events "+str(kernel)+" Kernel Density (h=" + str(bws) + "s)")
        #ax.plot(t[:, 0], norm_dens_kde_exp, '-', color='blue', zorder=3, label="On events Exponential Kernel Density (h="+str(bws)+"s)")
        #ax.plot(t[:, 0], norm_off_dens_kde, '--', color='black', zorder=3, label="Off events Gaussian Kernel Density (h="+str(bws)+"s)")
        ax.plot((t[:, 0] - t_start) * 1440., norm_off_dens_kde, '--', color='red',
                zorder=3)  #, label="Off events Epanechnikov Kernel Density (h="+str(bws)+"s)")
        #ax.plot(t[:, 0], norm_off_dens_kde_exp, '--', color='blue', zorder=3, label="Off events exponential Kernel Density (h="+str(bws)+"s)")
    else:
        ax.plot((t[:, 0]), norm_dens_kde, '--', color='red', zorder=3,
                label="On/Off events "+str(kernel)+" Kernel Density (h=" + str(bws) + "s)")
        ax.plot((t[:, 0]), norm_off_dens_kde, '--', color='red',
                zorder=3)  #, label="Off events Epanechnikov Kernel Density (h="+str(bws)+"s)")

    if MJD_flag == False:
        ax.plot((t[:, 0] - t_start) * 1440., norm_dens_kde - norm_off_dens_kde, '-', color='red',
                zorder=3, label="Flux "+str(kernel)+" Kernel Density (h=" + str(bws) + "s)")
    else:
        ax.plot((t[:, 0]), norm_dens_kde - norm_off_dens_kde, '-', color='red', zorder=3,
                label="Flux "+str(kernel)+" Kernel Density (h=" + str(bws) + "s)")

    #ea_hist, t_hist = np.histogram(t_evt[:,0], t_bins, weights=f_evt, normed=False)
    #ax.bar(t_hist[:-1], c_off_hist*alpha/(live_t_factor*120.), width=np.diff(t_bins), alpha=0.6, color='gray',label="Off flux")
    c_hist, t_hist = np.histogram(t_evt[:, 0], t_bins, weights=f_evt, normed=False)
    c_off_hist, t_off_hist = np.histogram(t_off_evt[:, 0], t_bins, weights=f_off_evt, normed=False)

    c2_hist, t2_hist = np.histogram(t_evt[:, 0], t_bins, weights=f_evt, density=True)
    c2_off_hist, t2_off_hist = np.histogram(t_off_evt[:, 0], t_bins, weights=f_off_evt, density=True)

    hist_norm = bw_hist * c_hist.sum() / (live_t_factor * bws_hist)
    hist_off_norm = bw_hist * c_off_hist.sum() / (live_t_factor * bws_hist) * alpha

    print "kernel:", norm_dens_kde.sum()
    print "raw sum histo:", c_hist.sum()
    print "length of on event list:", len(t_evt[:, 0])
    print "histo not normed:", c_hist.sum(), "approx count:", c_hist.sum() / np.mean(f_evt), "rate:", c_hist.sum() / (
        live_t_factor * bws)
    print "histo density:", c2_hist.sum() * hist_norm
    #print "bb normed:", c_bb.sum()
    print "off histo not normed:", c_off_hist.sum(), "approx count:", c_off_hist.sum() / np.mean(
        f_off_evt), "rate:", c_hist.sum() / (live_t_factor * bws)
    print c_hist.sum() / (live_t_factor * bws), c2_hist.sum() * hist_norm
    print "off", c2_off_hist.sum() * hist_off_norm
    print "total on events:", len(t_evt)
    print "total off events (before applying alpha):", len(t_off_evt)

    #ax.bar(t_hist[:-1], (c_hist-c_off_hist*alpha)/(live_t_factor*bws_hist), width=np.diff(t_bins), alpha=0.6, color='magenta',label="Flux histogram")
    #ax.bar(t_hist[:-1], (c_hist)/(live_t_factor*bws_hist), width=np.diff(t_bins), alpha=0.6, color='magenta',label="On Flux")
    if plot_hist==True:
        if MJD_flag == True:
            ax.bar(t2_hist[:-1], (c2_hist * hist_norm - c2_off_hist * hist_off_norm), width=np.diff(t_bins), alpha=0.5,
                   color='gray', label="Flux histogram")
        else:
            ax.bar((t2_hist[:-1] - t_start) * 1440., (c2_hist * hist_norm - c2_off_hist * hist_off_norm),
                   width=np.diff(t_bins) * 1440., alpha=0.5, color='gray', label="Flux histogram")

    if compfile != None:
        t_min = []
        f_min = []
        df_min = []
        with open(compfile) as fp:
            for line in fp:
                t_min.append(float(line.split()[0]))
                f_min.append(float(line.split()[1]))
                df_min.append(float(line.split()[2]))
        if MJD_flag == True:
            ax.errorbar(np.array(t_min), np.array(f_min), xerr=np.mean(np.diff(t_min))/2., yerr=np.array(df_min),
                        color='blue', fmt='o',label="Rigorous flux histogram")
            #            color='blue', fmt='o',label="Official VERITAS")
            #ax.bar(t2_hist[:-1], (c2_hist * hist_norm - c2_off_hist * hist_off_norm), width=np.diff(t_bins), alpha=0.6,
            #   color='gray', label="Flux histogram")
        else:
            ax.errorbar((np.array(t_min) - t_start) * 1440., np.array(f_min), xerr=np.mean(np.diff(t_min))/2.* 1440., yerr=np.array(df_min),
                        color='blue', fmt='o',label="Rigorous flux histogram")
                        #color='blue', fmt='o',label="Official VERITAS")
            #ax.errorbar((np.array(55740.453072) - t_start) * 1440., np.array(2.70743e-7), xerr=8, yerr=np.array(1.046446e-7), color='blue', fmt='o')

            #ax.bar((t2_hist[:-1] - t_start) * 1440., (c2_hist * hist_norm - c2_off_hist * hist_off_norm),
            #   width=np.diff(t_bins) * 1440., alpha=0.6, color='gray', label="Flux histogram")




    #ax.bar(t2_hist[:-1], (c2_hist*hist_norm), width=np.diff(t_bins), alpha=0.6, color='gray',label="On Flux histogram")

    #ax.errorbar((np.array(t_2min)-t_start)*1440., np.array(f_2min), xerr=120./86400*1440., yerr=np.array(df_2min), color='blue', fmt='o',label="Official VERITAS")
    #ax.errorbar((np.array(t_2min[-1])-t_start)*1440., np.array(f_2min[-1]), xerr=480./86400*1440., yerr=np.array(df_2min[-1]), color='blue', fmt='o')
    #ax.errorbar(t_2min, f_2min, xerr=np.array([60.]*len(t_2min)), yerr=np.array(df_2min), fmt='o')

    #ax.plot(t_evt[:,0],scale_f_evt,'-', color='gray', label="1/effective area")

    #print max(dens_kde)*1.1
    #ax.set_ylim(min(min(norm_dens_kde),min(norm_dens_kde_epanechnikov))*0.8, max(max(norm_dens_kde),max(norm_dens_kde_epanechnikov))*1.5)

    if MJD_flag == True:
        ax.set_xlabel('t (MJD)')
    else:
        ax.set_xlabel('t (minutes since the start of observations)')

    ax.set_ylabel('flux (m$^{-2}$ s$^{-1}$)')
    ax.legend(loc='best')

    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    fig.tight_layout()
    #plt.tight_layout()

    if plotname == None:
        #fig.show()
        plt.show()
    else:
        fig.savefig(plotname,format='eps', dpi=1000)

    if MJD_flag == True:
        return t[:, 0], norm_dens_kde - norm_off_dens_kde
    else:
        return (t[:, 0] - t_start) * 1440., norm_dens_kde - norm_off_dens_kde



