import numpy as np


def poisson_seq(total, mean_rate):

    """
    :param total: total number of points to simulate
    :param mean_rate: mean rate lambda in Poisson distr
                      Prob(t=wait time for 1 or more photon)=1-exp(-lambda*t)
                      1/mean_rate is the "beta" for numpy exponential random number
                      f(t;1/beta)=1/beta*exp(-t/beta)
    :return: a 1-d array containing time tags for each events
             of a Poisson process with known mean rate
    """

    # duration of the simulated sequence is approx. total/mean_rate
    
    # beta is scale param in exponential distr,
    # f(x;1/beta)=1/beta*exp(-x/beta)
    # 1/beta is the expected rate
    beta = 1./mean_rate

    t_seq = np.zeros(total)

    for t in range(0, total):
        k = np.random.exponential(beta)
        t_seq[t] = "%.10f" % (float(t_seq[t-1]+k))

    return t_seq


def sin_poisson_seq(total, mean_rate, t_sine, a_sine):

    '''
    :param total:       total number of points to simulate
    :param mean_rate:   mean rate lambda in Poisson distr
                        Prob(t=wait time for 1 or more photon)=1-exp(-lambda*t)
                        1/mean_rate is the "beta" for numpy exponential random number
                        f(t;1/beta)=1/beta*exp(-t/beta)
    :param t_sine:      period of the added sine modulation
    :param a_sine:      amplitude (same unit as mean rate) of the added sine modulation

    :return:            a 1-d array containing time tags for each events
                        of a Poisson process with known mean rate
    '''

    # beta is scale param in exponential distr,
    # f(x;1/beta)=1/beta*exp(-x/beta)
    # 1/beta is the expected rate
    beta = 1./mean_rate

    t_seq = np.zeros(total)

    for t in range(0, total):
        x_sine = 2.*a_sine*np.sin(2*np.pi/t_sine*t_seq[t-1])
        k = np.random.exponential(1./(1./beta+x_sine))
        t_seq[t] = "%.10f" % (float(t_seq[t-1]+k))

    return t_seq


def tri_poisson_seq(total, mean_rate, t_sine, a_sine):
    '''
    :param total:       total number of points to simulate
    :param mean_rate:   mean rate lambda in Poisson distr
                        Prob(t=wait time for 1 or more photon)=1-exp(-lambda*t)
                        1/mean_rate is the "beta" for numpy exponential random n
umber
                        f(t;1/beta)=1/beta*exp(-t/beta)
    :param t_sine:      period of the added triangle modulation
    :param a_sine:      amplitude (same unit as mean rate) of the added sine mod
ulation

    :return:            a 1-d array containing time tags for each events
                        of a Poisson process with known mean rate
    '''

    # beta is scale param in exponential distr,
    # f(x;1/beta)=1/beta*exp(-x/beta)
    # 1/beta is the expected rate
    beta = 1./mean_rate

    t_seq = np.zeros(total)

    for t in range(0, total):
        x_sine = -(a_sine/t_sine)*(t_seq[t-1]-int(t_seq[t-1]/t_sine)*t_sine)+mean_rate+0.5*a_sine
        k = np.random.exponential(1./(x_sine))
        t_seq[t] = "%.10f" % (float(t_seq[t-1]+k))

    return t_seq

def exp_poisson_seq(total, mean_rate, t0, t_exp, a_exp):
    '''
    :param total:       total number of points to simulate
    :param mean_rate:   mean rate lambda in Poisson distr
                        Prob(t=wait time for 1 or more photon)=1-exp(-lambda*t)
                        1/mean_rate is the "beta" for numpy exponential random n
umber
                        f(t;1/beta)=1/beta*exp(-t/beta)
    :param t_exp:      if t>t_0, flux(t) = mean_rate + a_exp * exp(-(t-t0)/t_exp)
                        e-folding time of the added pulse
    :param a_exp:      amplitude (same unit as mean rate) of the added pulse
ulation

    :return:            a 1-d array containing time tags for each events
                        of a Poisson process with known mean rate
    '''
    # beta is scale param in exponential distr,
    # f(x;1/beta)=1/beta*exp(-x/beta)
    # 1/beta is the expected rate
    beta = 1./mean_rate

    t_seq = np.zeros(total)

    for t in range(0, total):
        if t_seq[t-1] < t0:
          x = mean_rate
        else:
          x = mean_rate + a_exp * np.exp(-(t_seq[t-1]-t0)/t_exp)

        k = np.random.exponential(1./(x))
        t_seq[t] = "%.10f" % (float(t_seq[t-1]+k))

    return t_seq


