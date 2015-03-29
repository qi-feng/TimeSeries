__author__ = 'qfeng'
import lc
import numpy as np

N = 20
t = np.arange(0,N,1)
data = np.ones(N)
scale = 3
th = lc.lc(t,data)
print "rebin unnorm"
tr, xr = th.rebin(scale, norm=False)
print tr, xr
#print "rebin norm", th.rebin(scale, norm=True)

x = []
print th.t, th.x
#subh, subb = th.segmentation(seg_bin_num=5)
subh, subb = th.segmentation(n_seg=3)

for h, b in zip(subh, subb):
    sub_th = lc.lc(h, b)
    x.append(sub_th)
    print sub_th.t, sub_th.x

