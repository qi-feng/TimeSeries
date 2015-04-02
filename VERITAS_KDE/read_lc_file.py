__author__ = 'qfeng'

import numpy as np


def read_tfdf(filename=None):
    if filename == None:
        print "provide filename=? to read light curve from"
        raise SystemExit
    t = []
    flux = []
    df = []
    with open(filename) as f:
        for line in f:
            try:
                t.append(float(line.strip().split()[0]))
                flux.append(float(line.strip().split()[1]))
                df.append(float(line.strip().split()[2]))
            except IndexError:
                print "A line in the file doesn't have enough entries."
    return np.array(t), np.array(flux), np.array(df)