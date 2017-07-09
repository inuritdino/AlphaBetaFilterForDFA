#!/usr/bin/python3

"""Alpha-beta filter for DFA.  Following the paper Echeverria et
al. "Interpretation of heart rate variability via detrended
fluctuation analysis and ab filter.", Chaos 2003; 13(2):467-75
"""

import numpy as np
import matplotlib.pyplot as plt

def ab_filter(ls,Fs,num_ls=500):
    """
    Usage: me,ll = ab_filter(ls,Fs,num_ls=500)
    Inputs:
    ls - window sizes in log
    Fs - fluctuations in log
    num_ls - a number of window sizes to use in linear interpolation (default: 500).
    Outputs:
    me - estimation of the exponents as a function of window size
    ll - uniformly sampled window sizes on a log scale (ls)
    """
    ### Linear interpolation between the log10(ls) and log10(Fs) points 
    ### with uniformly sampled log(ls) points:
    uni_ls = np.linspace(start=ls[0],stop=ls[-1],num=num_ls)
    Fs_interp = np.interp(uni_ls,ls,Fs)

    ###############################
    ### Alpha-beta filter smoothing
    ###############################
    ## Ge -- log(Fs) estimator, Gp(k) -- log(Fs) predictor based on (k-1) value
    ## me -- d(log(Fs))/d(log(ls)) estimator, i.e. gradient
    Gp = np.zeros(len(Fs_interp))
    Ge = np.zeros(len(Fs_interp))
    me = np.zeros(len(Fs_interp))
    ## log(ls)-step/delta
    d = uni_ls[1] - uni_ls[0]
    ## The first estimators are the same as data
    Ge[0] = Fs_interp[0]
    me[0] = (Fs_interp[1] - Fs_interp[0])/(uni_ls[1] - uni_ls[0])
    ## The loop
    for i in range(1,len(Fs_interp)):
        Gp[i] = Ge[i-1] + d * me[i-1]
        r = Fs_interp[i] - Gp[i]
        Ge[i] = Gp[i] + alpha(i) * r
        me[i] = me[i-1] + (beta(i)/d) * r

    return (me,uni_ls);

def alpha(k):
    """Window size(k) dependent alpha value. 
    With suggestion from Echeverria et al. "Interpretation of
    heart rate variability..."
    """
    Q = 500
    if k > Q:
        return alpha(Q);
    else:
        return 2*(2*k-1)/(k*(k+1))

def beta(k):
    """Window size(k) dependent beta value.
    With suggestion from Echeverria et al. "Interpretation of
    heart rate variability..."
    """
    Q = 500
    if k > Q:
        return beta(Q);
    else:
        return 6/(k*(k+1));
