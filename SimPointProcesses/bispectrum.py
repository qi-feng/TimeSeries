__author__ = 'qfeng'

import numpy as np
from matplotlib import pyplot as plt
import lc

# function to calculate bispectrum etc
def bispec(t,x, bicoh=True, biphase=True):
    #assert np.mean(np.diff(t) - (t[1]-t[0])) < 1.e-6,\
    assert np.count_nonzero(np.diff(t) - (t[1]-t[0])) == 0,\
        "can't process unevenly binned light curve as of now"
    #important: subtract mean to have a zero-mean series
    x=x-np.mean(x)
    sp1 = np.fft.fftshift(np.fft.fft(x))
    freq1 = np.fft.fftshift(np.fft.fftfreq(t.shape[-1]))
    freq1 = freq1*len(t)/(float(t[-1])-float(t[0]))
    N = len(freq1)/2
    freq = freq1[N:]
    sp = sp1[N:]
    #sp = sp1
    print sp.dtype
    print np.count_nonzero(sp),"out of",len(sp),"Fourier coeffs are non-zero"
    #print freq
    #print sp[:100]
    bispec = np.zeros((N, N), complex)
    #if biphase == True:
    biphase_ = np.zeros((N, N))
    #if bicoh == True:
    #bicoh_ = np.zeros((N-1, N-1))
    # mod_xf1xf2sq = |X_fj * X_fk|^2
    # mod_xf1f2sq = |X_{fj+fk}|^2
    # bicoh^2 = bispec^2 / (mod_xf1xf2sq * mod_xf1f2sq)

    mod_xf1xf2sq = np.zeros((N, N))
    mod_xf1f2sq = np.zeros((N, N))

    for j in range(0,N):
        for i in range(0,j-1):
            if i+j >= N:
                continue
            else:
                bispec[i][j]=np.conjugate(sp[i+j])*sp[i]*sp[j]
                if bicoh == True:
                    #bispec[i][j]=np.conjugate(sp[i+j])*sp[i]*sp[j]/np.sqrt(abs(sp[i+j])**2*abs(sp[i]*sp[j])**2)
                    try:
                        #bicoh_[i][j] = float(abs(bispec[i][j]))/float((abs(sp[i+j])*abs(sp[i]*sp[j])))
                        #bicoh_[i][j] = np.sqrt(
                        #                      float( (abs(bispec[i][j]))**2 ) / float( (abs(sp[i+j]))**2 * (abs(sp[i]*sp[j]))**2 )
                        #                      )
                        mod_xf1xf2sq[i][j] = abs(sp[i]*sp[j])**2
                        mod_xf1f2sq[i][j] = abs(sp[i+j])**2
                    except ZeroDivisionError:
                        print "ERROR: float division by zero"
                        print sp[i], sp[j], sp[i+j]
                        continue
                if biphase == True:
                    biphase_[i][j] = np.arctan(bispec[i][j].imag/bispec[i][j].real)
    print bispec.dtype
    #return freq, bispec, bicoh_, biphase_
    return freq, bispec, mod_xf1xf2sq, mod_xf1f2sq, biphase_

def seg_bispec(list_of_sub_t, list_of_sub_x, bicoh=True, biphase=True):
    assert np.array(list_of_sub_t).shape == np.array(list_of_sub_x).shape,\
        "unequal lengths of input lists t and x for sub time series "
    n_seg = np.array(list_of_sub_t).shape[0]
    #print list_of_sub_t.shape
    freq1 = np.fft.fftshift(np.fft.fftfreq(list_of_sub_t[0].shape[-1]))
    freq1 = freq1*len(list_of_sub_t[0])/(float(list_of_sub_t[0][-1])-float(list_of_sub_t[0][0]))
    N = len(freq1)/2
    print N, np.array(list_of_sub_t).shape[0], np.array(list_of_sub_t).shape[-1]
    bisp = np.zeros((N, N), complex)
    biph = np.zeros((N, N))
    biph2 = np.zeros((N, N))
    bic = np.zeros((N, N))
    mod_xf1f2sq = np.zeros((N, N))
    mod_xf12sq = np.zeros((N, N))
    count = 0
    for tt, xx in zip(list_of_sub_t, list_of_sub_x):
        sub_freq, _bis, _bicf1f2, _bicf12, _biph = bispec(tt, xx, bicoh=bicoh, biphase=bicoh)
        bisp += _bis
        biph += _biph
        #bic += _bic
        mod_xf1f2sq += _bicf1f2
        mod_xf12sq += _bicf12
        count += 1
    assert n_seg == count
    mod_xf1f2sq = mod_xf1f2sq/float(n_seg)
    mod_xf12sq = mod_xf12sq/float(n_seg)
    bisp = bisp/float(n_seg)
    biph = biph/float(n_seg)

    for j in range(0,N):
        for i in range(0,j-1):
            if i+j >= N or mod_xf12sq[i][j]*mod_xf1f2sq[i][j]==0:
                continue
            else:
                bic[i][j]= abs(bisp[i][j])/np.sqrt(mod_xf12sq[i][j]*mod_xf1f2sq[i][j])
                biph2[i][j] = np.arctan(bisp[i][j].imag/bisp[i][j].real)

    return sub_freq, bisp, bic, biph, biph2


def ensemble_bispec(t, x, len_seg = 128, n_seg = 20, bicoh=True, biphase=True):
    assert np.count_nonzero(np.diff(t) - (t[1]-t[0])) == 0,\
        "can't process unevenly binned light curve as of now"
    len_shift = len_seg
    if len_seg + n_seg >= len(t):
        print "input light curve too short for dividing into",\
            n_seg, "segments of", len_seg,"sub intervals"
        raise(RuntimeError)
    elif len_seg * n_seg < len(t)-len_seg:
        print "we can make more than", n_seg, "segments from input light curve"
        print "change number of segments to", len(t)//len_seg
        n_seg = len(t)//len_seg
    elif len_seg * n_seg > len(t):
        len_shift = (len(t) - len_seg)//n_seg
    if len_seg < 32:
        print "do you think a 32-point segment is too short??"

    N = len_seg/2
    bisp = np.zeros((N, N), complex)
    biph = np.zeros((N, N))
    biph2 = np.zeros((N, N))
    bic = np.zeros((N, N))
    mod_xf1f2sq = np.zeros((N, N))
    mod_xf12sq = np.zeros((N, N))
    count = 0
    for i in range(0,n_seg):
        assert len_seg == len(t[i*len_shift:(i*len_shift+len_seg)]),str(i)+"th segment failed"
        sub_freq, _bis, _bicf1f2, _bicf12, _biph = bispec(t[i*len_shift:(i*len_shift+len_seg)], x[i*len_shift:(i*len_shift+len_seg)],
                                                          bicoh=bicoh, biphase=bicoh)
        bisp += _bis
        biph += _biph
        #bic += _bic
        mod_xf1f2sq += _bicf1f2
        mod_xf12sq += _bicf12
        count += 1

    assert n_seg == count
    mod_xf1f2sq = mod_xf1f2sq/float(n_seg)
    mod_xf12sq = mod_xf12sq/float(n_seg)
    bisp = bisp/float(n_seg)
    biph = biph/float(n_seg)

    for j in range(0,N):
        for i in range(0,j-1):
            if i+j >= N or mod_xf12sq[i][j]*mod_xf1f2sq[i][j]==0:
                continue
            else:
                bic[i][j]= abs(bisp[i][j])/np.sqrt(mod_xf12sq[i][j]*mod_xf1f2sq[i][j])
                biph2[i][j] = np.arctan(bisp[i][j].imag/bisp[i][j].real)

    #frac = abs(bisp)
    #denom = mod_xf12sq * mod_xf1f2sq
    #non_zero_indices = np.where(denom!=0)
    #bic[non_zero_indices] = np.sqrt(frac[non_zero_indices]/denom[non_zero_indices])
    ##bic = bic - 1./n_seg
    #biph2[non_zero_indices] = np.arctan(bisp[non_zero_indices].imag/bisp[non_zero_indices].real)

    return sub_freq, bisp, bic, biph/float(n_seg), biph2

