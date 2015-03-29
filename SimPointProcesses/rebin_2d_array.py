import numpy as np

__author__ = 'unknown'
def rebin(a, shape=None, scale=None):
    if shape == None:
        if a.shape[0]%scale > 0:
            print "padding the last few rows"
            #a = a[:(a.shape[0]//scale*scale)]
            print a[(a.shape[0]%scale-scale),:]
            a = np.concatenate((a, a[(a.shape[0]%scale-scale):,:]), 0 )
            print a
        if a.shape[1]%scale > 0:
            print "padding the last few columns"
            #a = a[:,:(a.shape[1]//scale*scale)]
            a = np.concatenate((a, a[:,(a.shape[1]%scale-scale):]), 1)
            print a
        sh = (a.shape[0]/scale),a.shape[0]//(a.shape[0]/scale),(a.shape[1]/scale),a.shape[1]//(a.shape[1]/scale)
    else:
        if a.shape[0] - a.shape[0]//shape[0]*shape[0] >0:
            print "throwing away the last few elements"
            a = a[:(a.shape[0]//shape[0]*shape[0])]
            print a
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]

    return a.reshape(sh).mean(-1).mean(1)
