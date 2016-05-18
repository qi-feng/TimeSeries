__author__ = 'qfeng'
import numpy as np
from sim_poisson_seq import *
import matplotlib.pyplot as plt
#from scipy import special
from scipy.stats import entropy
import pyfits
#import seaborn as sns
#sns.set_style("ticks")
#sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
from matplotlib.ticker import NullFormatter
import pandas as pd

try:
    from astroML.density_estimation.bayesian_blocks import bayesian_blocks
except:
    print "astroML is not found, can't do bayesian_blocks"



def get_lat_ea(plot=False):
    ea_hdu = pyfits.open('/Users/qfeng/HEAsoft/caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_FB.fits')
    #print ea_hdu.info()
    """
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU      10   ()
    1    EFFECTIVE AREA_FRONT  BinTableHDU     58   1R x 5C      ['74E', '74E', '32E', '32E', '2368E']
    2    PHI_DEPENDENCE_FRONT  BinTableHDU     61   1R x 6C      [23E, 23E, 8E, 8E, 184E, 184E]
    3    EFFICIENCY_PARAMS_FRONT  BinTableHDU     36   2R x 1C      [6E]
    4    EFFECTIVE AREA_BACK  BinTableHDU     58   1R x 5C      [74E, 74E, 32E, 32E, 2368E]
    5    PHI_DEPENDENCE_BACK  BinTableHDU     61   1R x 6C      [23E, 23E, 8E, 8E, 184E, 184E]
    6    EFFICIENCY_PARAMS_BACK  BinTableHDU     36   2R x 1C      [6E]
    """

    #print ea_hdu[0].header

    #front ea
    ea_data = ea_hdu[1].data
    """
        name = 'ENERG_LO'; format = '74E'; unit = 'MeV'
        name = 'ENERG_HI'; format = '74E'; unit = 'MeV'
        name = 'CTHETA_LO'; format = '32E'; unit = ''
        name = 'CTHETA_HI'; format = '32E'; unit = ''
        name = 'EFFAREA'; format = '2368E'; unit = 'm2'; dim = '(74, 32)'
    """

    ea_elo = ea_data[0][0]/1000.
    ea_ehi = ea_data[0][1]/1000.
    ea_coslo = ea_data[0][2]
    ea_coshi = ea_data[0][3]
    ea_ea = ea_data[0][4]

    ea_onaxis = ea_ea[-1]

    #back ea
    eab_data = ea_hdu[4].data
    eab_elo = eab_data[0][0]/1000.
    eab_ehi = eab_data[0][1]/1000.
    eab_coslo = eab_data[0][2]
    eab_coshi = eab_data[0][3]
    eab_ea = eab_data[0][4]

    eab_onaxis = eab_ea[-1]

    if plot:
        plt.errorbar((ea_elo+ea_ehi)/2., ea_onaxis, xerr=(ea_ehi-ea_elo)/2., fmt='.', color='r')
        plt.errorbar((eab_elo+eab_ehi)/2., eab_onaxis, xerr=(eab_ehi-eab_elo)/2., fmt='.', color='b')
        plt.errorbar((eab_elo+eab_ehi)/2., ea_onaxis+eab_onaxis, xerr=(eab_ehi-eab_elo)/2., fmt='.', color='k')

        plt.xscale('log')
        plt.show()

    return ea_onaxis+eab_onaxis

class powerlaw:
    # Class can calculate a power-law pdf, cdf from x_min to x_max,
    # always normalized so that the integral pdf is 1
    # (so is the diff in cdf between x_max and xmin)
    ############################################################################
    #                       Initialization of the object.                      #
    ############################################################################
    def __init__(self, a, x_min, x_max):
        """
        :param a: photon index, dN/dE = E^a
        :param x_min: e lo
        :param x_max: e hi
        :return:
        """
        self.a = a
        self.x_min = x_min
        self.x_max = x_max

    def pdf(self, x):
        if self.a==-1:
            return -1
        self.pdf_norm = (self.a+1)/(self.x_max**(self.a+1.0) - self.x_min**(self.a+1.0))
        return self.pdf_norm*x**self.a

    def cdf(self, x):
        if self.a==-1:
            return -1
        self.cdf_norm = 1./(self.x_max**(self.a+1.0) - self.x_min**(self.a+1.0))
        return self.cdf_norm*(x**(self.a+1.0) - self.x_min**(self.a+1.0))

    def ppf(self, q):
        if self.a==-1:
            return -1
        norm = (self.x_max**(self.a+1.0) - self.x_min**(self.a+1.0))
        return (q*norm*1.0+self.x_min**(self.a+1.0))**(1.0/(self.a+1.0))

    def random(self, n):
        r_uniform = np.random.random_sample(n)
        return self.ppf(r_uniform)


class discan:
    ############################################################################
    #                       Initialization of the object.                      #
    ############################################################################
    def __init__(self, theta1=0, theta2=0):
        self.theta1 = theta1
        self.theta2 = theta2

    def read_data_tte(self, ts, es):
        assert len(ts)==len(es), "Please provide arrival times and energies of the same dimension"
        self.tdata = np.array(ts)
        self.edata = np.array(es)
        self.theta1_data = 0
        self.theta2_data = 0

    def gen_pulse(self, t, A, tzero, tau1, tau2):
        """
        :param t: a positive 1-D np array of time
        "param A: peak intensity
        :param tau1: see Norris 2005 (in Scargle 07)
        :param tau2:
        :return:
        """
        # lists are pulse specific, t and It are for the entire LC
        self.pulse_tau1 = [tau1]
        self.pulse_tau2 = [tau2]
        self.pulse_mu = [np.sqrt(tau1/tau2*1.0)]
        self.pulse_lam = [np.exp(2.*self.pulse_mu[0])]
        self.t = t
        self.pulse_A = [A]
        self.It = np.zeros(len(t))
        self.It[np.where(t<tzero)] = 0.0
        t_offset = t[np.where(t>tzero)] - tzero
        self.pulse_It = [A*self.pulse_lam[0]*np.exp(-(tau1/t_offset)-(t_offset/tau2))]
        self.It[np.where(t>tzero)] = self.pulse_It[0]
        self.pulse_w = [tau2*np.sqrt(1.0+4.0*self.pulse_mu[0])]
        self.pulse_k = [tau2/self.pulse_w[0]]
        self.pulse_peak = [np.sqrt(tau1*tau2)]
        self.pulse_tzero = [tzero]

    def add_pulse(self, A, tzero, tau1, tau2):
        if not hasattr(self, 't'):
            print "Have to run gen_pulse first"
            return
        self.pulse_tau1.append(tau1)
        self.pulse_tau2.append(tau2)
        self.pulse_mu.append(np.sqrt(tau1/tau2*1.0))
        self.pulse_lam.append(np.exp(2.*self.pulse_mu[-1]))
        self.pulse_A.append(A)
        t_offset = self.t[np.where(self.t>tzero)] - tzero
        self.pulse_It.append(A*self.pulse_lam[-1]*np.exp(-(tau1/t_offset)-(t_offset/tau2)))
        self.It[np.where(self.t>tzero)] += self.pulse_It[-1]
        #self.It += A*self.pulse_lam[-1]*np.exp(-(tau1/self.t)-(self.t/tau2))
        self.pulse_w.append(tau2*np.sqrt(1.0+4.0*self.pulse_mu[-1]))
        self.pulse_k.append(tau2/self.pulse_w[-1])
        self.pulse_peak.append(np.sqrt(tau1*tau2))
        self.pulse_tzero.append(tzero)

    def throw_arrival_times(self, Nsim=4000, sort=True):
        if not hasattr(self, 't'):
            print "Have to run gen_pulse first, or fill self.t and self.It"
            return
        if not hasattr(self, 'Nsim'):
            self.Nsim = Nsim
        else:
            assert self.Nsim == Nsim, "You are throwing inconsistent number of arrival times and energies"
        print "Throwing Nsim="+str(Nsim)+" arrival times"
        self.tsims = throw_arrival_times(self.t, self.It, Nsim, sort=sort)
        if sort:
            self.tsims = np.sort(self.tsims)

    def throw_energies(self, nu, E_min, E_max, Nsim=4000, add_gaussian_sigma=None, add_gaussian_sigma_hi=None):
        """
        :param nu:
        :param E_min:
        :param E_max:
        :param Nsim:
        :param add_gaussian_sigma: the fraction of sigma_E to add to E
        :param add_gaussian_sigma_hi: if specified, the added sigma_E changes from add_gaussian_sigma at E_min to
                                        add_gaussian_sigma_hi at E_max
        :return:
        """
        ### note esim_true is the true power-law energy thrown,
        ### while esim is the "detected" energy with an E_rms added
        if not hasattr(self, 't'):
            print "Have to run gen_pulse first, or fill self.t and self.It"
            return
        if not hasattr(self, 'Nsim'):
            self.Nsim = Nsim
        else:
            assert self.Nsim == Nsim, "You are throwing inconsistent number of arrival times and energies"
        print "Throwing Nsim="+str(Nsim)+" energies"
        self.Emin = E_min
        self.Emax = E_max
        pl_nu = powerlaw(nu, E_min, E_max)
        self.esims = pl_nu.random(Nsim)
        self.esims_true = self.esims
        if add_gaussian_sigma is not None:
            if add_gaussian_sigma_hi is not None:
                for i, e in enumerate(self.esims):
                    sigma_e = add_gaussian_sigma+1.0*(e-E_min)*(add_gaussian_sigma_hi-add_gaussian_sigma)/(E_max-E_min)
                    self.esims[i] += np.random.random()*sigma_e*e
            else:
                for i, e in enumerate(self.esims):
                    self.esims[i] += np.random.random()*add_gaussian_sigma*e

    def add_dispersion(self, theta1=None, theta2=None, Emin=None, sort=True):
        ### dispersions are added to the arrival times based on esims_true, not esims (the detected energies)
        ### however, the energies detected still need to be esims
        if not hasattr(self, 'tsims'):
            print "Have to run throw_arrival_times and throw_energies first, or fill self.tsims and self.esims"
            return
        if Emin is not None:
            self.Emin = Emin
        if not hasattr(self, 'Emin'):
            print "Have to specify an Emin"
            return
        if theta1 is not None:
            self.theta1 = theta1
        if theta2 is not None:
            self.theta2 = theta2
        self.tsims_shifted = np.zeros(self.tsims.shape)
        for i, e in enumerate(self.esims_true):
            self.tsims_shifted[i] = self.tsims[i] + self.theta1*(e-self.Emin) + self.theta2**(e-self.Emin)**(e-self.Emin)
        if sort:
            self.esims_shifted = self.esims[np.argsort(self.tsims_shifted)]
            self.tsims_shifted = np.sort(self.tsims_shifted)
        else:
            self.esims_shifted = self.esims
        return self.tsims_shifted

    def sub_dispersion(self, theta1=None, theta2=None, Emin=None, sort=True):
        ### dispersions are subtracted based on edata (or esims in the case of simulation test),
        ###  as we can't know esims_true in the real world
        if not hasattr(self, 'tdata'):
            print "Have to run read_data_tte first, or fill self.tdata and self.edata"
            return
        if Emin is not None:
            self.Emin = Emin
        if not hasattr(self, 'Emin'):
            print "Have to specify an Emin"
            return
        if theta1 is not None:
            self.theta1_data = theta1
        if theta2 is not None:
            self.theta2_data = theta2
        self.tdata_shifted = np.zeros(self.tdata.shape)
        for i, e in enumerate(self.edata):
            self.tdata_shifted[i] = self.tdata[i] - self.theta1_data*(e-self.Emin) - self.theta2_data*(e-self.Emin)*(e-self.Emin)
        if sort:
            self.edata_shifted = self.edata[np.argsort(self.tdata_shifted)]
            self.tdata_shifted = np.sort(self.tdata_shifted)
        else:
            self.edata_shifted = self.edata
        return self.tdata_shifted

    def get_bayesian_blocks_bins(self, ts, p0=0.05):
        return bayesian_blocks(ts, p0=p0)

    def tte2hist(self, ts, bins):
        """
        return rates and bin centers
        """
        count, bins = np.histogram(ts, bins, normed=False)
        return np.array(count)/np.array((bins[1:]-bins[:-1])), np.array(bins[1:]+bins[:-1])/2.

    def make_cells(self, ts, zero_thresh=1.e-5, get_time=False):
        if np.sum(np.diff(ts)<0):
            print "input ts not sorted!"
            ts = np.sort(ts)
        cell_rates = []
        cell_t = []
        cell_widths = []
        this_interval = 0
        this_count = 0
        intervals = np.diff(ts)
        bin_widths = np.zeros(ts.shape)
        bin_widths[1:-1] = (intervals[:-1]+intervals[1:])/2.
        bin_widths[0] = intervals[0]
        bin_widths[-1] = intervals[-1]
        #cell_rates = 1./bin_widths
        small_bin = False
        bin_counter = -1
        for i in bin_widths:
            bin_counter +=1
            if i>=zero_thresh and not small_bin:
                cell_rates.append(1./i)
                if get_time:
                    cell_widths.append(i)
                    cell_t.append(ts[bin_counter])
            else:
                this_count += 1
                this_interval += i
                small_bin = True
                if this_interval >= zero_thresh:
                    cell_rates.append(1.0*this_count/this_interval)
                    if get_time:
                        cell_widths.append(this_interval)
                        cell_t.append(ts[bin_counter]+i/2.-this_interval/2.)
                    #print this_count, this_interval, 1.0*this_count/this_interval
                    this_interval = 0
                    this_count = 0
                    small_bin = False
        if get_time:
            return np.array(cell_rates), np.array(cell_t), np.array(cell_widths)
        return np.array(cell_rates)


    def shannon_entropy(self, rates):
        rates = np.array(rates)
        sum_rate = np.sum(rates)
        I_shannon = 0.0
        one = 0.0
        for rate in rates:
            if rate == 0:
                print "A rate is zero"
                #raise RuntimeError
                continue
            I_shannon += -1.0*rate/sum_rate * np.log2(1.0*rate/sum_rate)
        return I_shannon

    def total_variation(self, rates):
        rates = np.array(rates)
        sum_rate = np.sum(rates)
        if np.sum(rates)<=0:
            print "the sum of rates provided is negative or zero"
            return np.inf
        tv = 0.0
        ps = rates*1.0/sum_rate
        #for rate in rates:
        tv = np.sum(abs(np.diff(ps)))
        return tv


    def cost_function(self, rates, method='shannon'):
        if method=='shannon':
            #print "Calculating home-brewed Shannon entropy"
            return self.shannon_entropy(rates)
        if method=='entropy':
            #p_rates = np.array(rates)/np.sum(np.array(rates))
            #return entropy(p_rates)
            #print "Calculating scipty.stats.entropy"
            return entropy(rates)
        if method=='totalVariation':
            #flip sign to minimize
            return -self.total_variation(rates)


    def recover_dispersion1(self, theta1s, bins=None, t_start=None, t_stop=None, cost_function='entropy',
                            hist_method='bb', p0=0.05, sort=True, norm_entropy_Nbins=False):
        """
        returns best_theta1, shannon_data, and 1-D array of shannons as a function of theta1s
        """
        if hist_method=='bb':
            bins_data = bayesian_blocks(self.tdata, p0=p0)
            rate_data, bins_data = self.tte2hist(self.tdata, bins_data)
            cost_function_data = self.cost_function(rate_data, method=cost_function)
            print len(rate_data), "Bayesian blocks were made for data"
            bins_sims = bayesian_blocks(self.tsims, p0=p0)
            rate_sims, bins_sims = self.tte2hist(self.tsims, bins_sims)
            cost_function_sims = self.cost_function(rate_sims, method=cost_function)
            print len(rate_sims), "Bayesian blocks were made for sims with no dispersion"
        elif hist_method=='cell':
            cell_rates = self.make_cells(self.tdata)
            cost_function_data = self.cost_function(cell_rates, method=cost_function)
        else:
            rate_data, bins_data = self.tte2hist(self.tdata, bins)
            cost_function_data = self.cost_function(rate_data, method=cost_function)
        cost_functions = []
        for theta1 in theta1s:
            tdata_sub_dispersed = self.sub_dispersion(theta1=theta1, sort=sort)
            if hist_method=='bb':
                bins_bb = bayesian_blocks(tdata_sub_dispersed, p0=p0)
                data_sub_dispersed_rate, data_shifted_bins = self.tte2hist(tdata_sub_dispersed, bins_bb)
                print len(data_sub_dispersed_rate), "Bayesian blocks were made when testing theta1=", theta1
                #plt.plot(data_shifted_bins, data_sub_dispersed_rate)
                #plt.show()
            elif hist_method=='cell':
                #data_sub_dispersed_rate = self.make_cells(tdata_sub_dispersed)
                data_sub_dispersed_rate, cell_t, cell_widths = self.make_cells(tdata_sub_dispersed, get_time=True)
                #cost_function_data = self.cost_function(cell_rates, method=cost_function)
                #plt.errorbar(cell_t, data_sub_dispersed_rate, xerr=cell_widths, label="Theta1="+str(theta1))
                #plt.title(self.cost_function(data_sub_dispersed_rate, method=cost_function))
                #plt.legend(loc='best')
                #plt.show()
            else:
                data_sub_dispersed_rate, data_shifted_bins = self.tte2hist(tdata_sub_dispersed, bins)
                #plt.plot(data_shifted_bins, data_sub_dispersed_rate, label="Theta1="+str(theta1))
                #plt.title(self.cost_function(data_sub_dispersed_rate, method=cost_function))
                #plt.legend(loc='best')
                #plt.show()
            #print data_sub_dispersed_rate
            cost_function_data_sub_dispersed = self.cost_function(data_sub_dispersed_rate, method=cost_function)
            if norm_entropy_Nbins:
                cost_function_data_sub_dispersed = cost_function_data_sub_dispersed*1.0/np.log2(1.0*len(data_sub_dispersed_rate))
            cost_functions.append(cost_function_data_sub_dispersed)
            #print "theta1:", theta1, "cost", cost_function_data_sub_dispersed
        cost_functions = np.array(cost_functions)
        #print "Looped over theta1s, cost function values", cost_functions
        min_cost = np.min(cost_functions[np.isfinite(cost_functions)])
        #print "Min cost function value:", min_cost
        #print "index:", np.where(cost_functions==min_cost)
        best_theta1 = theta1s[np.where(cost_functions==min_cost)]
        #print "Theta 1 for min cost function value",
        return best_theta1, cost_function_data, cost_functions


def throw_arrival_times(t, It, n, t_start=None, t_stop=None, sort=True):
    """
    :param t: 1-D array of LC sampling time
    :param It: 1-D array of intensities of the same dimension of t
    :param n: number of events to simulate
    :param t_start: optional boundaries
    :param t_stop:
    :return: 1-D array of n arrival times that follows the pdf of It
    """
    t = np.array(t)
    It = np.array(It)
    assert len(t) == len(It), "Input time and intensities are of different sizes"
    if t_start is None:
        t_start = t[0]
    if t_stop is None:
        t_stop = t[-1]
    Imax = np.max(It)
    ts = []
    i_printed = -1
    while len(ts)<n:
        #if len(ts) % 1000 == 0 and len(ts) != i_printed:
        #    i_printed = len(ts)
        #    print str(i_printed)+" arrival times simulated so far..."
        r1 = np.random.random()
        r2 = np.random.random()
        # uniform random from t_start to t_stop
        t_random = (1.0*r1*(t_stop-t_start)+t_start)
        # if ti < t_random < t_i+1, use I(ti)
        ti = t[np.argmax(t>t_random)-1]
        #print "t_sim", t_random, "It sim", r2*Imax, "It real", It[np.argmax(t>t_random)-1]
        # r2*Imax is uniform random from 0 to Imax
        if r2*Imax <= It[np.argmax(t>t_random)-1]:
            ts.append(t_random)

    ts = np.array(ts)
    if sort:
        ts = np.sort(ts)
    return ts

def plot_t_e(ts, es, bins=100, pcolor='k', hcolor='gray', ptitle=None, ax=None, show=True):
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    if ax is None:
        plt.figure(1, figsize=(8, 8))
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
    else:
        plt.axes(ax)
        ax.axes = plt.axes(rect_scatter)
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)


    # the scatter plot:
    axScatter.scatter(ts, es, s=5, color=pcolor)

    xymax = np.max([np.max(np.fabs(ts)), np.max(np.fabs(es))])

    axScatter.set_xlim((np.min(ts), np.max(ts)))
    axScatter.set_ylim((np.min(es), np.max(es)))
    plt.title(ptitle)

    #bins = np.arange(-lim, lim + binwidth, binwidth)
    axHistx.hist(ts, bins=bins, histtype='step', color=hcolor, linewidth=2)
    axHisty.hist(es, bins=np.logspace(np.log10(np.min(es)), np.log10(np.max(es)), bins), orientation='horizontal', histtype='step', color='gray', linewidth=2)

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.set_xscale('log')
    axHisty.set_yscale('log')
    axScatter.set_yscale('log')
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    axScatter.set_ylabel('log E(GeV)')
    axScatter.set_xlabel('t (ms)')

    #ax.plot(tsims, esims, 'k.')
    #plt.tight_layout()
    if show:
        plt.show()


def visualize_powerlaw_example():
    fig, ax = plt.subplots(1, 1)
    a = -2
    b = -4
    x_max=1000
    x_min=10
    pla = powerlaw(a, x_min, x_max)
    plb = powerlaw(b, x_min, x_max)

    x = np.linspace(x_min, x_max, 100)

    N=100000
    r_pl_a = pla.random(N)
    r_pl_b = plb.random(N)

    #print r_pl_a, r_pl_b
    #print min(r_pl_a), max(r_pl_a)
    #print min(r_pl_b), max(r_pl_b)


    ax.plot(x, pla.pdf(x), 'r-', lw=5, alpha=0.6, label=('powerlaw pdf a=%.1f'%a))
    ax.plot(x, plb.pdf(x), 'b-', lw=5, alpha=0.6, label=('powerlaw pdf a=%.1f'%b))
    #ax.plot(x, pla.cdf(x), 'r--', lw=5, alpha=0.6, label=('powerlaw cdf a=%.1f'%a))
    #ax.plot(x, plb.cdf(x), 'b--', lw=5, alpha=0.6, label=('powerlaw cdf a=%.1f'%a))
    print pla.pdf_norm, plb.pdf_norm

    ax.hist(r_pl_a, color='r', histtype='step', bins=100, normed=True, linewidth=2)
    ax.hist(r_pl_b, color='b', histtype='step', bins=100, normed=True, linewidth=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-10, 1)
    plt.legend(loc='best')
    plt.show()

def read_bat_lc(f, t_trigger=156822675.776, t_start=-0.17, t_stop=0.43):
    df = pd.read_csv(f, sep=r"\s+", skiprows=[1, 2])
    df.TIME -= t_trigger
    df = df[df.TIME>=t_start]
    df = df[df.TIME<=t_stop]
    return df.TIME[:]-df.TIME.values[0], df.TOTCOUNTS[:]

def test_discan(Ntrial=1, Nsim=4400, plot_intermediate=False, bins=200, p0=0.01, theta1_true=20.,
                norm_entropy_Nbins=False, costfunction='entropy',
                plotname="theta1_entropy.png", histmethod='cell', theta1s=np.arange(0., 40., 0.1)):
    """
    :param Ntrial:   #number of trials to recover a theta1
    :param Nsim:     #number of TTEs to throw
    :param plot_intermediate:
    :param plotname:
    :param histmethod:
    :return:
    """
    #theta1s = np.arange(0., 40., 0.25)

    dc = discan()

    choice=1
    ######################
    # choice 1:
    # make a LC throwing pulses t, It
    if choice==1:
        t = np.arange(0.,600,1.)
        As =[18.8, 19.3, 62.4, 141.6, 105.7, 120.3, 16.7, 102.3, 15.3]
        #teffs are the times that the pulses reach 0.01 of their peak intensities
        #teffs = [111, 172, 204, 256, 290, 314, 332, 394, 460]
        tzeros = [0, 110, 155, 230, 285, 310, 255, 390, 365]
        #BATSE taus
        #tau1s = [5.64E+03, 2.45E+03, 8.57E+02, 3.79E+02, 2.55E+00, 3.42E+00, 2.77E+03, 6.31E+00, 4.65E+03]
        #tau2s = [4.24E+00, 2.60E+00, 6.44E+00, 3.64E+00, 1.43E+01, 1.28E+01, 5.30E+00, 2.12E+01, 3.58E+00]
        #LAT taus should be narrower by a factor of ~0.1
        tau1s = [5.64E+02, 2.45E+02, 8.57E+1, 3.79E+01, 2.55E-1, 3.42E-1, 2.77E+02, 6.31E-1, 4.65E+02]
        tau2s = [4.24E-1, 2.60E-1, 6.44E-1, 3.64E-1, 1.43E+0, 1.28E+0, 5.30E-1, 2.12E+0, 3.58E-1]

        #print t
        count = 0
        for A, tzero, tau1, tau2 in zip(As, tzeros, tau1s, tau2s):
            if count == 0:
                dc.gen_pulse(t, A, tzero, tau1, tau2)
            else:
                dc.add_pulse(A, tzero, tau1, tau2)
            count += 1

        print "A", dc.pulse_A
        print "w", dc.pulse_w
        print "peak", dc.pulse_peak
    ######################

    ######################
    # choice 2: read BAT LC:
    if choice==2:
        t, It = read_bat_lc("/Users/qfeng/Data/SWIFT/GRB051221A/00173780000/bat/event/GRB051221A_15_350keV_1ms_t1.csv")
        #convert to ms
        dc.t = t*1000.
        dc.It = It
    ######################

    ######################
    # choice 3: throw a fat pulse:
    if choice==3:
        bkg_rate = 10.
        bkg_rms = 3.
        t = np.arange(0,600,1)
        dc = discan()
        dc.gen_pulse(t, 50, 80, 30, 80)
        dc.It += np.random.random(size=len(t))*bkg_rms+bkg_rate

    ######################


    #tsims = throw_arrival_times(dc.t, dc.It, Nsim)

    #ea = get_lat_ea()
    tau_QG = 17. #ms/GeV
    E_min = 0.03
    E_max = 10.
    nu = -4./3
    #pl_nu = powerlaw(nu, E_min, E_max)
    #esims = pl_nu.random(Nsim)

    best_theta1_trials = np.zeros(Ntrial)


    allCostFunctions = np.zeros((len(theta1s), Ntrial))
    sort=True

    for trial in range(Ntrial):
        dc.throw_arrival_times(Nsim=Nsim, sort=sort)
        #dc.throw_energies(nu, E_min, E_max, Nsim=Nsim)
        # add a gaussian as energy uncertainties in the detector
        dc.throw_energies(nu, E_min, E_max, Nsim=Nsim, add_gaussian_sigma=0.1, add_gaussian_sigma_hi=0.2)
        #plot_t_e(dc.tsims, dc.esims)
        #sim_rate40, sim_bins_40 = dc.tte2hist(dc.tsims, 40)

        dc.add_dispersion(theta1=theta1_true, sort=sort)
        #plot_t_e(dc.tsims_shifted, dc.esims)

        dc.read_data_tte(dc.tsims_shifted, dc.esims_shifted)
        #plot_t_e(dc.tdata, dc.edata)
        #sim_shifted_rate40, sim_shifted_bins40 = dc.tte2hist(dc.tsims_shifted, 40)

        if histmethod=='cell':
            best_theta1, data_cost_function, cost_functions = dc.recover_dispersion1(theta1s,
                                                                                     cost_function=costfunction, sort=sort,
                                                                                     hist_method='cell')
        if histmethod=='bb' or histmethod=='bayesian_blocks':
            best_theta1, data_cost_function, cost_functions = dc.recover_dispersion1(theta1s,
                                                                                     norm_entropy_Nbins=norm_entropy_Nbins,
                                                                                     cost_function=costfunction, sort=sort,
                                                                                     hist_method='bb', p0=p0)
        if histmethod=='hist':
            best_theta1, data_cost_function, cost_functions = dc.recover_dispersion1(theta1s,
                                                                                     cost_function=costfunction, sort=sort,
                                                                                     bins=bins, hist_method='hist')

        print "best_theta1 that is not nan or inf", best_theta1[0]

        print "data cost function", data_cost_function
        #for t_, s_ in zip(theta1s, cost_functions):
        #    print "theta1:", t_, "cost function:", s_

        dc.sub_dispersion(theta1=best_theta1[0], sort=sort)
        best_theta1_trials[trial] = best_theta1[0]
        #fig, ax = plt.subplots(2,2)
        show=True
        if plot_intermediate:
            plot_t_e(dc.tsims, dc.esims, pcolor='k', hcolor='gray', ptitle="Sim", show=show) #, ax=ax[1,0])
            ##plot_t_e(dc.tsims_shifted, dc.esims_shifted, pcolor='r', hcolor='r') #, ax=ax[1,1], show=False)
            plot_t_e(dc.tdata, dc.edata, pcolor='r', hcolor='r', ptitle="Dispersion added", show=show) #, ax=ax[0,0], show=False)
            plot_t_e(dc.tdata_shifted, dc.edata_shifted, pcolor='b', hcolor='b', ptitle="Dispersion recovered", show=show) #, ax=ax[0,1], show=False)
            plt.plot(theta1s, cost_functions)
            plt.axvline(theta1_true, linewidth=2, color = 'b', linestyle='--')
            plt.xlabel(r'$\theta_{QG}$')
            plt.ylabel('Entropy')
            plt.show()

        allCostFunctions[:,trial] = cost_functions

    #cell_rates, cell_t, cell_w = dc.make_cells(dc.tdata, get_time=True)
    #print cell_rates.shape
    #plt.errorbar(cell_t, cell_rates, xerr=cell_w/2.)

    plt.hist(best_theta1_trials, bins=20)
    plt.axvline(theta1_true, linewidth=2, color = 'b', linestyle='--')
    plt.xlabel(r'$\theta_{QG}$')
    #plt.ylabel('Counts')
    plt.savefig(plotname.split('.')[0]+"_hist.png", fmt='png')
    #plt.show()
    plt.clf()

    cf_hi = np.percentile(allCostFunctions, 95, axis=1)
    cf_lo = np.percentile(allCostFunctions, 5, axis=1)
    plt.fill_between(theta1s, cf_lo, cf_hi, alpha=.4, color='g',label="95% CI")
    plt.plot(theta1s, np.mean(allCostFunctions, axis=1))
    plt.axvline(theta1_true, linewidth=2, color = 'b', linestyle='--')
    plt.xlabel(r'$\theta_{QG}$')
    plt.ylabel('Entropy')
    plt.savefig(plotname, fmt='png')
    plt.clf()
    #plt.savefig("discan/theta_entropy_true_theta20_bayesian_blocks_100trials.png", fmt='png')
    #plt.show()

    return best_theta1_trials

def scan_theta1s():
    Ntrial = 100
    theta1_trues = np.arange(10, 31, 1)
    best_theta1_alltrials = np.zeros((len(theta1_trues), Ntrial))
    for i, theta1_true in enumerate(theta1_trues):
        #Hist:
        #best_theta1_alltrials[i,:] = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False, histmethod='hist',theta1s=np.arange(0., 40., 0.1),
        #        bins=100, plotname="discan/theta_entropy_true_theta"+str(theta1_true)+"_hist_Nbins20_Nphoton200_100trial.png")
        #cell:
        best_theta1_alltrials[i,:] = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False,
                        histmethod='cell',theta1s=np.arange(0., 40., 0.1),
                        plotname="discan/cell/theta_entropy_true_theta"+str(theta1_true)+"_cell_Nphoton200_100trial.png")

    plt.errorbar(theta1_trues, np.mean(best_theta1_alltrials, axis=1), yerr=np.std(best_theta1_alltrials, axis=1), fmt='.')
    plt.plot(theta1_trues, theta1_trues, linewidth=1, color = 'gray', linestyle='--')
    plt.xlabel(r'$\theta_{true}$')
    plt.ylabel(r'$\theta_{recovered}$')
    #plt.savefig("discan/scan_theta1_profile_hist_100bins_Nphoton200.png", fmt='png')
    plt.savefig("discan/scan_theta1_profile_cell_Nphoton200.png", fmt='png')
    plt.show()


    for k in range(len(theta1_trues)):
        print np.array([theta1_trues[k]]*Ntrial), best_theta1_alltrials[k]
        plt.plot(np.array([theta1_trues[k]]*Ntrial), best_theta1_alltrials[k], 'k.', ms=5)

    plt.plot(theta1_trues, theta1_trues, linewidth=1, color = 'gray', linestyle='--')
    plt.xlabel(r'$\theta_{true}$')
    plt.ylabel(r'$\theta_{recovered}$')
    #plt.savefig("discan/scan_theta1_scatterplot_hist_100bins_Nphoton200.png", fmt='png')
    plt.savefig("discan/scan_theta1_scatterplot_cell_Nphoton200.png", fmt='png')
    plt.show()

    #test_discan(Ntrial=100, Nsim=100, plot_intermediate=False, histmethod='cell', plotname="discan/theta_entropy_true_theta20_cell_100trial.png")
    #test_discan(Ntrial=100, Nsim=200, plot_intermediate=False, histmethod='hist',theta1s=np.arange(0., 40., 0.1),
    #            bins=100, plotname="discan/theta_entropy_true_theta20_hist_Nbins20_Nphoton200_100trial.png")
    #test_discan(Ntrial=100, Nsim=200, plot_intermediate=False, histmethod='bb', theta1s=np.arange(0., 40., 0.1),
    #            p0=0.1, plotname="discan/theta_entropy_true_theta20_bb_p0_0p1_Nphoton200_100trial.png")
    #test_discan(Ntrial=1, Nsim=200, plot_intermediate=True, theta1s=np.arange(0., 40., 0.1), histmethod='bb', p0=0.01, plotname="discan/theta_entropy_true_theta20_bb_1trial.png")

if __name__=="__main__":
    Ntrial = 100
    #theta1_trues = np.arange(10, 31, 1)
    theta1_true = 20
    best_theta1_trials = np.zeros(Ntrial)
    #Hist:
    #best_theta1_alltrials[i,:] = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False, histmethod='hist',theta1s=np.arange(0., 40., 0.1),
    #        bins=100, plotname="discan/theta_entropy_true_theta"+str(theta1_true)+"_hist_Nbins20_Nphoton200_100trial.png")
    #total variation
    #best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False, bins=100,
    #                histmethod='hist',theta1s=np.arange(0., 40., 0.1), p0=0.01, norm_entropy_Nbins=False, costfunction='totalVariation',
    #                plotname="discan/hist/theta_totalVariation_true_theta"+str(theta1_true)+"_hist_Nbins100_Nphoton200_"+str(Ntrial)+"trials.png")

    #cell:
    #best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False,
    #                histmethod='cell',theta1s=np.arange(0., 40., 0.1),
    #                plotname="discan/cell/theta_entropy_true_theta"+str(theta1_true)+"_cell_Nphoton200_10trials.png")
    #cell:
    best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=50, theta1_true=theta1_true, plot_intermediate=False,
                    histmethod='cell',theta1s=np.arange(0., 40., 0.1),
                    plotname="discan/cell/theta_entropy_true_theta"+str(theta1_true)+"_cell_Nphoton50_"+str(Ntrial)+"trials.png")
    #total variation
    #best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=100, theta1_true=theta1_true, plot_intermediate=False,
    #                histmethod='cell',theta1s=np.arange(0., 40., 0.1), costfunction='totalVariation',
    #                plotname="discan/cell/theta_totalVariation_true_theta"+str(theta1_true)+"_cell_Nphoton100_"+str(Ntrial)+"trials.png")

    #bb bayesian blocks:
    #best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False,
    #                histmethod='bb',theta1s=np.arange(0., 40., 0.1), p0=0.01, norm_entropy_Nbins=True,
    #                plotname="discan/bb/theta_entropy_true_theta"+str(theta1_true)+"_bb_p0_0p01_Nphoton200_"+str(Ntrial)+"trials.png")

    #total variation
    #best_theta1_trials = test_discan(Ntrial=Ntrial, Nsim=200, theta1_true=theta1_true, plot_intermediate=False,
    #                histmethod='bb',theta1s=np.arange(0., 40., 0.1), p0=0.01, norm_entropy_Nbins=True, costfunction='totalVariation',
    #                plotname="discan/bb/theta_totalVariation_true_theta"+str(theta1_true)+"_bb_p0_0p01_Nphoton200_"+str(Ntrial)+"trials.png")

