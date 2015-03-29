#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
from astroML.time_series import generate_power_law

from sim_poisson_seq import poisson_seq

total = 300000
mean_rate = 30.

t_seq = poisson_seq(total,mean_rate)

nbins = 200
weight_seq=np.zeros(total)
print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))

ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")
#ax.plot(t_seq,1./beta+2.*Asine*np.sin(2*np.pi/Tsine*t_seq))

fname = 'simPoissonTTE_'+str(total)+'sample_rate'+str(mean_rate)+'_withE.txt'
f = open(fname, 'w')
for x in t_seq:
  f.write(str(x)+'\t'+str("1.0")+'\n')

f.close()

plt.show()
