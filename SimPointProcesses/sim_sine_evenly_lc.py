__author__ = 'qfeng'

import numpy as np
import matplotlib.pyplot as plt
import lc
from scipy.interpolate import interp1d

N = 1000
t_start = 0
dt = 1.0
mean_rate = 60.
period = 3.
amplitude = 23.24

t = np.arange(t_start, t_start + (N+1)*dt, dt)
x = amplitude * np.sin(t/period) + mean_rate
f_interp = interp1d(t, x, kind='cubic')
t_x = np.arange(t_start, t_start + (N)*dt, dt/5.)

slc = lc.lc(t, x)
slc.write_to_file(filename='perfect_sine_period'+str(period)+'mean'+str(mean_rate)+'.txt')


#plt.plot(t_x,f_interp(t_x), 'r.')
plt.plot(t,x, 'bo')

plt.show()
