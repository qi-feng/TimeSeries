__author__ = 'qfeng'

import numpy as np
from astroML.time_series import generate_power_law
import matplotlib.pyplot as plt

from astroML.time_series import lomb_scargle
from astroML.time_series import multiterm_periodogram

def generate_pl_pp(beta=1.0, T=1024, dt=1., mean_rate=0.15, std_dev=0.03):
    """
    simulate a point process that has a red noise PSD
    :param beta: slope of the red noise (the power-law index in P~f^-beta)
    :param T: duration of the light curve simulated
    :param dt: time step in the sim LC
    :param mean_rate: mean rate of the sim LC
    :param std_dev: std of the sim LC
    :return: a sequence of the arrival time
    """
    #t_seq = np.arange(0,N*dt,dt)
    assert(dt>0)
    N = int(T/dt)
    print N

    y = generate_power_law(N=N, dt=dt, beta=beta)
    y = y/np.std(y)*std_dev
    y = y + mean_rate
    total_num_points = int(T*mean_rate)
    t_seq = np.zeros(total_num_points)
    for k in range(1, total_num_points):
        try:
            index=int(t_seq[k-1]/dt)
        except:
            print("Funny index for {kk} th bin at t={tt}".format(kk=k,tt=t_seq[k-1]))
        if index>=N-1:
            index=N-1
        t_k = np.random.exponential(1./y[index])
        t_seq[k] = "%.10f" % (float(t_seq[k-1]+t_k))

    return t_seq

if __name__ == "__main__":
    beta = 1.0
    mean_rate=0.15
    std_dev=0.03
    T = 16000
    t_seq = generate_pl_pp(beta=beta, T=T, mean_rate=mean_rate, std_dev=std_dev)

    nbins = 200
    weight_seq=np.zeros(len(t_seq))
    print "Duration comparison: ",t_seq[-1]-t_seq[0],len(t_seq)/mean_rate
    weight_seq[:]=float((nbins)/(t_seq[-1]-t_seq[0]))

    fig = plt.figure(figsize=(5, 8))
    ax = fig.add_subplot(211)

    n, bins, patches = ax.hist(t_seq, bins=nbins, weights=weight_seq,histtype='stepfilled', alpha=0.2, normed=False,label="Standard histogram")

    y = (bins[:-1]+bins[1:])/2.

    #N=1024
    dt=1.
    N = int(T/dt)
    t_seq2 = np.arange(0,N*dt,dt)
    y2 = generate_power_law(N=N, dt=dt, beta=beta)
    y2 = y2/np.std(y2)*std_dev + mean_rate

    ax.plot(y, n, 'ro')
    ax.plot(t_seq2, y2)
    #omega=np.linspace(float(1./(t_seq[-1]-t_seq[0])),float(0.5/t_seq[1]-t_seq[0]),int((t_seq[-1]-t_seq[0])/2))
    #omega=np.linspace(float(1./(t_seq[-1]-t_seq[0])),float(0.5/(t_seq[1]-t_seq[0])),int((t_seq[-1]-t_seq[0])/2))
    omega=np.linspace(float(1./(y[-1]-y[0])),float(0.5/(y[1]-y[0])),int((y[-1]-y[0])/(2*(y[1]-y[0]))))
    omega2 = np.linspace(float(1./(t_seq2[-1]-t_seq2[0])),float(0.5/(t_seq2[1]-t_seq2[0])),int((t_seq2[-1]-t_seq2[0])/2))
    print len(omega)
    n = n - np.mean(n)
    P_LS = lomb_scargle(y, n, np.sqrt(np.abs(n)), omega , generalized=True)
    P_M = multiterm_periodogram(y, n, np.sqrt(np.abs(n)), omega)
    #P_M = multiterm_periodogram(n, y, np.sqrt(np.abs(y)), omega)
    P_LS2 = lomb_scargle(t_seq2, y2, np.sqrt(np.abs(y2)), omega2 , generalized=True)
    P_M2 = multiterm_periodogram(t_seq2, y2, np.sqrt(np.abs(y2)), omega2)

    ax2 = fig.add_subplot(212)
    ax2.plot(omega,P_LS,'r',label="Generalized Lomb-Scargle Periodogram")
    ax2.plot(omega,P_M, 'b',label="Multiterm Periodogram")
    ax2.plot(omega2,P_LS2,'g')
    ax2.set_xlabel(r'Freq (Hz)')
    ax2.set_ylabel(r'Power')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    plt.tight_layout()
    plt.show()





