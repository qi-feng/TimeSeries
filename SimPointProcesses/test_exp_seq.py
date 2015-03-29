#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
import scipy
from astroML.plotting import hist

from sim_poisson_seq import exp_poisson_seq

total = 9000
mean_rate = 60.
t_0 = 15.
t_exp = 30.
a_exp = 23.24
#a_exp = 77.5

t_seq = exp_poisson_seq(total, mean_rate, t_0, t_exp, a_exp)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)

nbins = 200
weight_seq=np.zeros(total)
print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))

fname = 'simPoisson+exp_TTE_'+str(total)+'sample_rate'+str(mean_rate)+'exp_t0'+str(t_0)+'P'+str(t_exp)+'expA'+str(a_exp)+'_withE.txt'
#f = open(fname, 'w')
with open(fname) as f:
    for x in t_seq:
        f.write(str(x)+'\t'+str("1.0")+'\n')

#f.close()


ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")
#ax.plot(t_seq,mean_rate+2.*a_exp*np.sin(2*np.pi/t_exp*t_seq))
x_list = []
for t in t_seq:
  if t>t_0:
    x_list.append(mean_rate + a_exp * np.exp(-(t-t_0)/t_exp))
  else:
    x_list.append(mean_rate)

ax.plot(t_seq,x_list)
plt.show()
