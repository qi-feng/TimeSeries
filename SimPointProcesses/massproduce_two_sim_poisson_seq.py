import numpy as np

from sim_poisson_seq import poisson_seq

#n_lc = 1000
n_start = 0
n_stop = 100
#total = 10000
total = 3000
mean_rate = 0.15
Efake = 0.1
Efake2 = 1.0

#for i in range(0,n_lc):
for i in range(n_start,n_stop):
	t_seq = poisson_seq(total,mean_rate)
        t_seq2 = poisson_seq(total,mean_rate)
	
	print "Duration comparison: ",t_seq[-1]-t_seq[0],total/mean_rate
	
	fname = 'simTwoPoissonTTE_'+str(total)+'sample_rate'+str(mean_rate)+'_withE'+str(Efake)+'and'+str(Efake2)+'trial'+str(i)+'.txt'
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

#plt.show()
