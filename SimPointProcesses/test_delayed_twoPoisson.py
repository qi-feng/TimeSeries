#%matplotlib inline
#!!!!! unfinished

import matplotlib.pyplot as plt
import numpy as np
from astroML.time_series import lomb_scargle
from astroML.time_series import multiterm_periodogram
from astroML.time_series import generate_power_law
from astroML.time_series import ACF_scargle , ACF_EK


from sim_poisson_seq import poisson_seq

#total = 10000
#mean_rate = 10.0
total = 3000
mean_rate = 0.15
bg_rate = 0.05
sig_rate = 0.1

Efake = 0.1
Efake2 = 1.0

t_seq = poisson_seq(total,sig_rate)
t_seq_2 = t_seq + 30.0
t_seq_bg = poisson_seq(total,bg_rate)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)

nbins = 200
weight_seq=np.zeros(total)
weight_seq2=np.zeros(total)

print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))
weight_seq2[:]=float((nbins)/(t_seq2[-1]-t_seq2[0]))

ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")
#ax.plot(t_seq,1./beta+2.*Asine*np.sin(2*np.pi/Tsine*t_seq))
ax.hist(t_seq2, bins=nbins, weights=weight_seq2,histtype='stepfilled', alpha=0.2,     normed=False, color='r', label="Standard histogram")

fname = 'simPoissonTTE_'+str(total)+'sample_rate'+str(mean_rate)+'_withE'+str(Efake)+'and'+str(Efake2)+'.txt'

f = open(fname, 'w')
n1 = len(t_seq)
n2 = len(t_seq2)
i=0
j=0

while i<n1 and j<n2:
  if float(t_seq[i])<=t_seq2[j]:
    f.write(str(t_seq[i])+'\t'+str(Efake)+'\n')
    i += 1
  else:
    f.write(str(t_seq2[j])+'\t'+str(Efake2)+'\n')
    j += 1

f.close()

plt.show()
