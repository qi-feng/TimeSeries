#%matplotlib inline

import numpy as np
#from astroML.time_series import lomb_scargle
#from astroML.time_series import multiterm_periodogram
#from astroML.time_series import generate_power_law
#from astroML.time_series import ACF_scargle , ACF_EK

from sim_rednoise_pointprocess import generate_pl_pp
#from sim_poisson_seq import poisson_seq
import pandas as pd

energy_file='M4V20140429_OnEvents_T_E.txt'
df = pd.read_csv(energy_file, header=None, delim_whitespace=True)
df.columns = ['t','e']
e_list = df['e']

n_lc = 100
n_start = 0
n_stop = 10
#total = 10000
mean_rate = 0.15
T = 16000 
dt = 0.1/mean_rate
std_dev = 0.03
beta = 1.0

for i in range(0,n_lc):
#for i in range(n_start,n_stop):
        t_seq = generate_pl_pp(beta=beta, T=T, dt=dt, mean_rate=mean_rate, std_dev=std_dev)
        fname = 'simRedNoiseTTE_duration'+str(T)+'_rate'+str(mean_rate)+'_withEreadFromM4Obs_trial'+str(i)+'.txt'
        f = open(fname, 'w')
        n1 = len(t_seq)
        for i, x in enumerate(t_seq):
            f.write(str(x)+'\t'+str(e_list[i%len(e_list)])+'\n')

        f.close()


