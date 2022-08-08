import os
import sys
import scipy as sp
from scipy.optimize import curve_fit
import numpy as np


def getbounceperiod(t, x, omb):
    """Calculates the bouncing period using a time series of latitudes.
    omb is the estimate of the bounce frequency."""

    # initial parameters [amplitude, angular_frequency, phase]
    pin = [max(abs(x)), omb, 0]

    # what to fit
    dp = [1, 1, 1]

    stol = 1e-6
    niter = 200
    minstep = [0, 0, 0]
    maxstep = [1e10, 1e10, 1e10]
    options = [minstep, maxstep]

    # weights on ydata
    wt = np.ones(t.shape)

    F = lambda t, *p: p[0] * np.sin(p[1] * t + p[2])
    # dFdp = lambda t,p: [np.sin(p[1]*t+p[2]),p[0]*t*np.cos(p[1]*t+p[2]),p[0]*np.cos(p[1]*t+p[2])]

    p, covp = curve_fit(F, t, x, p0=pin, sigma=wt)

    # standard deviation for estimated parameters
    psd = np.sqrt(np.diag(covp))

    tb = 2 * np.pi / p[1]
    dtb = 2 * np.pi * psd[1] / p[1] ** 2
    mplat = max(abs(x))
    fitlat = F(t, *p)

    print(
        f"Bounce: tb={tb:.2f} +- {dtb:.2f} s lm={mplat:.2f} deg A={p[0]:.2f} +- {psd[0]:.2f}"
    )
    return tb, dtb, mplat, fitlat
