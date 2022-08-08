
from pymagdisc.bfield.magneticField3D import mdiscMagneticField3D
from pymagdisc.vis.plot_compare_mdisc_dipole import plot_compare_mdisc_dipole
from pymagdisc.vis.plot_compare_contours_alpha import plot_compare_contours_alpha
from pymagdisc.data.load_data import load_model
from pymagdisc import config
from pymagdisc.tracer.mdbtracer import MDBTracer
from scipy.interpolate import RectBivariateSpline as interp2
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    mdfile = 'jup_mdisc_kh3e7_rmp90fix.mat' # magnetodisc mat-file
    partype = 'p';                          # proton
    Ep = 100;                               # energy 100 MeV
    Ri = 4;                                 # initial equatorial distance 4*Rp
    ai = 30;                                # initial pitch angle 30 degrees
    timespec = [0,0,2,0];                   # run for 2 dipole bounce periods
    savefile = 'my_jup_mdisc_sim';          # name for result mat-file
    npertc = 5                              # number of Boris iterations per gyroperiod
    MD = load_model(config.PATH_TO_DATA + mdfile)
    tracer_jup = MDBTracer(mdfile=mdfile,partype=partype,Ep=Ep,Ri=Ri,ai=ai,timespec=timespec,savefile=savefile)
    results = tracer_jup.run_simul()
