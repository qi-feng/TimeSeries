#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
from astroML.time_series import lomb_scargle
from astroML.time_series import multiterm_periodogram
from astroML.time_series import generate_power_law
from astroML.time_series import ACF_scargle , ACF_EK


from sim_poisson_seq import poisson_seq

#total = 10000
total = 3000
mean_rate = 0.15
Efake = 0.1

t_seq = poisson_seq(total,mean_rate)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)

nbins = 200
weight_seq=np.zeros(total)
print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))

ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")
#ax.plot(t_seq,1./beta+2.*Asine*np.sin(2*np.pi/Tsine*t_seq))

fname = 'simPoissonTTE_'+str(total)+'sample_rate'+str(mean_rate)+'_withE'+str(Efake)+'.txt'
f = open(fname, 'w')
for x in t_seq:
  f.write(str(x)+'\t'+str(Efake)+'\n')

f.close()

plt.show()
