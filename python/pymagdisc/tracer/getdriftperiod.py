import os
import sys
import scipy as sp
import numpy as np
from scipy.optimize import curve_fit


def getdriftperiod(t, x, omb, a4td):
    """Calculates the drift period using a time series of longitudes.
    omb is the estimate of the bounce frequency and a4td the angle for 1 drift period.
    360 for longitude in degrees, 2*pi in radians."""

    # function [kvg,td,dtd,f]=getdriftperiod(t,x,omb,a4td)

    # fit straight line
    pol = np.polyfit(t, x, 1, rcond=None, full=False, w=None, cov=False)
    # slope is estimate for drift angular frequency
    omd = pol[0]

    # initial parameters
    # [modulation amplitude, modulation_period, phase, drift_frequency]
    pin = [max(abs((x - np.polyval(pol, x)) / omd)), 2 * omb, 0, omd]

    # what to fit
    dp = [1, 1, 1, 1]
    # fixed omb
    # dp = [1, 0, 1, 1]

    # Parameters
    stol = 1e-6
    niter = 200
    minstep = [0, 0, 0, 0]
    maxstep = [np.infty, np.infty, np.infty, np.infty]
    options = [minstep, maxstep]

    # These are probably not needed
    # t = t(:);
    # x = x(:);
    wt = np.ones(t.shape)

    F = lambda t, *p: p[3] * (t + p[0] * np.sin(p[1] * t + p[2]))
    # dFdp = lambda t : [p[4]*np.sin(p[2]*t+p[3]),p[4]*p[0]*t*np.cos(p[1]*t+p[2]),p[3]*p[0]*np.cos(p[1]*t+p[2]),t+p[0]*np.sin(p[1]*t+p[2])]

    ##This has to be modified!!!!
    p, covp = curve_fit(F, t, x, p0=pin)

    # standard deviation for estimated parameters
    # psd = np.zeros(p.shape)
    psd = np.sqrt(np.diag(covp))

    # omfit=2*omb thus tfit=tb/2
    tb = 2 * (2 * np.pi / p[1])
    dtb = 2 * (2 * np.pi * psd[1] / p[1] ** 2)

    td = a4td / p[3]
    dtd = a4td * psd[3] / p[3] ** 2
    fitlon = F(t, *p)

    print(f"Drfit: td={td:.2f}+-{dtd:.2f}")
    return td, dtd, fitlon
