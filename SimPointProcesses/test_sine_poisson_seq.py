#%matplotlib inline

import matplotlib.pyplot as plt
import numpy as np
import scipy
from astroML.plotting import hist

from sim_poisson_seq import sin_poisson_seq

total = 9000
mean_rate = 60.
t_sine = 30.
a_sine = 6.

t_seq = sin_poisson_seq(total, mean_rate, t_sine, a_sine)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)

nbins = 200
weight_seq=np.zeros(total)
print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))

fname = 'simPoisson+sine_TTE_'+str(total)+'sample_rate'+str(mean_rate)+'sineP'+str(t_sine)+'sineA'+str(a_sine)+'_withE.txt'
f = open(fname, 'w')
for x in t_seq:
    f.write(str(x)+'\t'+str("1.0")+'\n')

f.close()


ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")
ax.plot(t_seq,mean_rate+2.*a_sine*np.sin(2*np.pi/t_sine*t_seq))
plt.show()
