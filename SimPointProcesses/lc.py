#!/usr/bin/env python
__author__ = 'qfeng'

import numpy as np
from matplotlib import pyplot as plt

class hist(object):
    def __init__(self, *args):
        if len(args) == 2:
            self.fill_hist(args[0], args[1])

    def fill_hist(self, hist, bin_edges):
        self.hist = np.array(hist)
        self.bin_edges = np.array(bin_edges)
        self.bin_center = (bin_edges[:-1] + bin_edges[1:])/2.
        self.bin_width = np.diff(bin_edges)
        self.bin_num = len(bin_edges) - 1
        self.duration = bin_edges[-1] - bin_edges[0]

    def rebin(self, scale, norm=True, verbose=False):
        # merge n=scale bins into one bin
        # new histogram has a length of self.bin_num//scale or self.bin_num//scale + 1
        """
        :param scale: how many bins to merge into one new bin
        :param norm: if True (default), returns rate; other wise counts
        :param verbose:
        :return: new_hist and new_bin_edges, similar to style of numpy hist
        """
        new_bin_num = self.bin_num//scale
        if verbose==True:
            print "Rebin every",scale,"bins into a new bin"
            print new_bin_num, "full new bins"

        if self.bin_num - self.bin_num//scale*scale > 0:
            print "the last bin is partially filled"
            print "only", self.bin_num - self.bin_num//scale*scale, "old bins in the last new bin"
            new_bin_num = self.bin_num//scale+1

        new_bin_edges = np.zeros(new_bin_num+1)
        new_hist = np.zeros(new_bin_num)
        x = np.array([self.hist[i:self.bin_num//scale*scale:scale] for i in range(0,scale)])

        if len(self.bin_edges[0::scale]) == new_bin_num+1:
        #if self.bin_num//scale+1 == new_bin_num:
            new_bin_edges[:] = self.bin_edges[0::scale]
            new_hist[:] = np.array([x[:,k].sum() for k in range(0,new_bin_num)])

        else:
        #if len(self.bin_edges[0::scale]) == new_bin_num:
            new_bin_edges[:-1] = self.bin_edges[0::scale]
            new_bin_edges[-1] = self.bin_edges[-1]
            new_hist[:-1] = np.array([x[:,k].sum() for k in range(0,new_bin_num-1)])
            new_hist[-1] = np.sum(self.hist[self.bin_num//scale*scale:])

        if norm == True:
            new_hist[:] = new_hist[:]/np.diff(new_bin_edges)

        return new_hist, new_bin_edges

    def segmentation(self, n_seg=None, seg_bin_num=None, verbose=False):
        # divide the histogram to n_seg equal-length shorter ones,
        # each sub histogram has a length of self.bin_num//n_seg
        # there may be a (n_seg + 1)th sub histogram that has a length
        # of self.bin_num - (self.bin_num//n_seg) * n_seg
        """
        :param n_seg: how many segments to divide into
        :param verbose:
        :return: list_of_sub_hist and list_of_sub_bin_edges
        """
        if n_seg==None and seg_bin_num==None:
            print "Please provide the number of equal-length segments n_seg="
            print "or provide the number of bins for each segment seg_bin_num="
        elif seg_bin_num == None:
            seg_bin_num = self.bin_num//n_seg
        elif n_seg == None:
            n_seg = self.bin_num//seg_bin_num

        if verbose==True:
            print "segmentation: Dividing the histogram into",n_seg,"equal-length segments"
            print "each segment has",seg_bin_num,"bins"

        list_of_sub_hist = []
        list_of_sub_bin_edges = []

        for i in range(0,n_seg):
            sub_hist = self.hist[(i*seg_bin_num):((i+1)*seg_bin_num)]
            sub_bin_edges = self.bin_edges[(i*seg_bin_num):((i+1)*seg_bin_num + 1)]
            list_of_sub_hist.append(sub_hist)
            list_of_sub_bin_edges.append(sub_bin_edges)

        if self.bin_num - (self.bin_num//n_seg) * n_seg > 0:
            print "In call of segmentation: the", n_seg+1, "th segment is shorter"
            print "has only", self.bin_num - (self.bin_num//n_seg) * n_seg, "bins"
            last_hist = self.hist[(n_seg*seg_bin_num):]
            last_bin_edges = self.bin_edges[(n_seg*seg_bin_num):]
            list_of_sub_hist.append(last_hist)
            list_of_sub_bin_edges.append(last_bin_edges)


        return list_of_sub_hist, list_of_sub_bin_edges

    def plot(self, ax=None, xlabel="Time", ylabel="Histogram", label=None, color='gray', alpha=0.5, fname=None):

        if ax == None:
            subplot = 111
            fig = plt.figure(figsize=(8, 5))
            fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
            ax = fig.add_subplot(subplot)
            doplot = True
        else:
            doplot = False

        ax.bar(self.bin_edges[:-1], self.hist, width=self.bin_width, alpha=alpha,
                color=color, label=label)

        if doplot == True:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend(loc='best')
            #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            fig.tight_layout()

            if fname == None:
                plt.show()
            else:
                fig.savefig(fname,format='eps', dpi=1000)

class lc(object):
    def __init__(self, *args):
        if len(args) == 2:
            self.fill_lc(args[0], args[1])

    def fill_lc(self, t, x):
        assert len(t) == len(x), "Input t and x have different lengths"
        self.x = np.array(x)
        self.t = np.array(t)
        self.t_intervals = np.diff(t)
        self.t_num = len(t)
        self.duration = t[-1] - t[0]

    def rebin(self, scale, norm=True, verbose=False):
        assert np.count_nonzero(np.diff(self.t) - (self.t[1]-self.t[0])) == 0,\
               "can't rebin unevenly binned light curve as of now"
        new_t_num = self.t_num//scale
        if verbose==True:
            print "Combining every",scale,"sampling times into a new sampling time"
            print new_t_num, "full new sampling times"

        if self.t_num - self.t_num//scale*scale > 0:
            print "the last new time interval is partially filled"
            print "only", self.t_num - self.t_num//scale*scale, "old time intervals in the last new interval"
            new_bin_num = self.t_num//scale+1

        new_t = np.zeros(new_t_num)
        new_x = np.zeros(new_t_num)
        x = np.array([self.x[i:self.t_num//scale*scale:scale] for i in range(0,scale)])
        t = np.array([self.t[i:self.t_num//scale*scale:scale] for i in range(0,scale)])

        if len(self.t[0::scale]) == new_t_num:
            new_t[:] = np.array([t[:,k].mean() for k in range(0,new_t_num)])
            new_x[:] = np.array([x[:,k].sum() for k in range(0,new_t_num)])

        else:
        #if len(self.bin_edges[0::scale]) == new_bin_num:
            new_t[:-1] = np.array([t[:,k].mean() for k in range(0,new_t_num-1)])
            new_t[-1] = np.mean(self.t[self.t_num//scale*scale:])
            new_x[:-1] = np.array([t[:,k].sum() for k in range(0,new_t_num-1)])
            new_x[-1] = np.sum(self.t[self.t_num//scale*scale:])

        if norm == True:
            new_x[:] = new_x[:]/scale

        return new_t, new_x

    def segmentation(self, n_seg=None, seg_sample_num=None, verbose=False):
        # divide the histogram to n_seg equal-length shorter ones,

        """
        :param n_seg: how many segments to divide into
        :param verbose:
        :return: list_of_sub_t and list_of_sub_x
        """
        if n_seg==None and seg_sample_num==None:
            print "Please provide the number of equal-length segments n_seg="
            print "or provide the number of bins for each segment seg_sampe_num="
        elif seg_sample_num == None:
            seg_sample_num = self.t_num//n_seg
        elif n_seg == None:
            n_seg = self.t_num//seg_sample_num

        if verbose==True:
            print "segmentation: Dividing the histogram into",n_seg,"equal-length segments"
            print "each segment has",seg_sample_num,"bins"

        list_of_sub_t = []
        list_of_sub_x = []

        for i in range(0,n_seg):
            sub_x = self.x[(i*seg_sample_num):((i+1)*seg_sample_num)]
            sub_t = self.t[(i*seg_sample_num):((i+1)*seg_sample_num)]
            list_of_sub_x.append(sub_x)
            list_of_sub_t.append(sub_t)

        if self.t_num - (self.t_num//n_seg) * n_seg > 0:
            print "In call of segmentation: the", n_seg+1, "th segment is shorter"
            print "has only", self.t_num - (self.t_num//n_seg) * n_seg, "bins"
            last_x = self.x[(n_seg*seg_sample_num):]
            last_t = self.t[(n_seg*seg_sample_num):]
            list_of_sub_x.append(last_x)
            list_of_sub_t.append(last_t)

        return list_of_sub_t, list_of_sub_x


    def interpolate(self, interval=None):
        ### !! fix me: the first two sampling defines the timer interval
        assert np.count_nonzero(np.diff(self.t) - (self.t[1]-self.t[0])) != 0,\
        "Appears that the light curve is evenly sampled, nothing to interpolate"
        from scipy.interpolate import interp1d
        if interval == None:
            interval = self.t[1]- self.t[0]
        t_all = np.arange(self.t[0], self.t[-1] + interval, interval)
        t_set = set(t_all)
        t_missing = sorted( t_set.difference(self.t))

        for tx in t_missing:
            indices = np.where(t_all == tx)
            in2 = np.where(self.t == t_all[indices[0]-1])
            self.t = np.insert(self.t, in2[0]+1, (t_all[indices[0]]))
            self.x = np.insert(self.x, in2[0]+1, (self.x[in2[0]]+self.x[in2[0]+1])/2.)
            #xt_comp =  np.insert(xt, in2[0]+1, f_interp(t_all[indices[0]-1]))

    def write_to_file(self, filename=None):
        if filename == None:
            print "provide filename=? to write light curve to"
            raise SystemExit
        with open(filename, 'w') as f:
            for t, x in zip(self.t, self.x):
                f.write(str(t)+'\t'+str(x)+'\n')


    def plot(self, ax=None, xlabel="Time", ylabel=None, label=None, color='gray', alpha=0.5, fname=None):

        if ax == None:
            subplot = 111
            fig = plt.figure(figsize=(8, 5))
            fig.subplots_adjust(bottom=0.08, top=0.95, right=0.95, hspace=0.1)
            ax = fig.add_subplot(subplot)
            doplot = True
        else:
            doplot = False

        ax.plot(self.t, self.x, alpha=alpha,
                color=color, label=label)

        if doplot == True:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend(loc='best')
            #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            fig.tight_layout()

            if fname == None:
                plt.show()
            else:
                fig.savefig(fname,format='eps', dpi=1000)

def make_hist(t_tte, t_start=None, bin_width=None, bin_num=None):
    """
    make a histogram
    :param t_tte:
    :param bin_width:
    :param bin_num:
    :param t_start:
    :return:
    """
    if t_start==None:
        t_start = t_tte[0]
    if bin_width==None and bin_num==None:
        print "provide one of the two parameters: bin_width= or bin_num="
    elif bin_num==None:
        bin_num = (t_tte[-1]-t_start)//bin_width+1
    elif bin_width==None:
        bin_width = float(t_tte[-1]-t_start)/float(bin_num)
    t_end = t_start + bin_num*bin_width

    bins = np.arange(t_start,t_end+bin_width,bin_width)
    _hist, _bins = np.histogram(t_tte,bins=bins)
    h = hist(_hist, _bins)
    return h