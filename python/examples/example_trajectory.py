import matplotlib.pyplot as plt
import sys
import os

sys.path.append("../")
import pymagdisc
from pymagdisc.bfield.magneticField3D import mdiscMagneticField3D
from pymagdisc.vis.plot_compare_mdisc_dipole import plot_compare_mdisc_dipole
from pymagdisc.vis.plot_compare_contours_alpha import plot_compare_contours_alpha
from pymagdisc.data.load_data import load_model
from pymagdisc import config
from pymagdisc.tracer.mdbtracer import MDBTracer
from scipy.interpolate import RectBivariateSpline as interp2
from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    ### Change the following set of parameters
    mdfile = "jup_mdisc_kh3e7_rmp90.mat"  # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
    partype = "S+"
    # particle type (select: p, e, O+, O++, S+, S++, S+++)
    Ep = 25
    # initial energy (in MeV)
    Ri = 70
    # initial equatorial distance (in Jovian radii)
    ai = 10
    # initial pitch angle (in degrees)
    timespec = [0, 0, 20, 0]
    # run for x dipole bounce periods (modify the third element)
    npertc = 20  # number of Boris iterations per gyroperiod (~ points used per bounce)
    mdisctype = (
        "exp"  # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
    )

    # --- Do not modify the following unless necessary! ---
    # The code block below creates a text ("string" in python terminology), based on the parameters above
    runname = (
        partype + "Ri" + str(Ri) + "_Ep" + str(Ep) + "_ai" + str(ai) + "" + mdisctype
    )
    # This creates the name of the file
    # Choose to save the results in a "mat" file
    savefile = runname + ".mat"

    # Load the magnetic field model, run the code and save the results in the file specified above
    MD = load_model("../../" + mdfile)
    tracer_jup = MDBTracer(
        mdfile=mdfile,
        partype=partype,
        Ep=Ep,
        Ri=Ri,
        ai=ai,
        timespec=timespec,
        savefile=savefile,
        npertc=npertc,
    )
    results = tracer_jup.run_simul()
    # Save the file in the folder 'results'
    # savemat('./results/'+savefile,results)

    # Plot the trajectory of the particle
    fig = plt.figure(figsize=(4, 4), dpi=120)
    plt.plot(results["Rcylb"] / results["Re"], results["Zb"] / results["Re"], "-")
    plt.xlabel(r"rc ($R_J$)")
    plt.ylabel(r"z ($R_J$)")
    plt.xlim(
        min(results["Zb"] / results["Re"] + 2), max(results["Zb"] / results["Re"] + 2)
    )

    plt.title(r"Trajectory: {}, E={} MeV, Ri={} RJ, ai={}°".format(partype, Ep, Ri, ai))

    circle1 = plt.Circle((0, 0), 1, color="k")
    plt.gca().add_patch(circle1)
    plt.axis("equal")
    # plt.savefig('./plots/trajectory_'+runname+'.png',bbox_inches="tight")
    plt.show()

    for i in range(0, len(results["tb"])):
        if (
            np.sqrt((results["Zb"][i]) ** 2 + (results["Rcylb"][i]) ** 2) / results["Re"]
            <= 1.0
        ):
            print("Particle lost !")

    fig = plt.figure(figsize=(4, 4), dpi=120)
    plt.plot(results["tb"] / 60, results["pitchan"])
    plt.xlabel("Time (min)")
    plt.ylabel(r"$\alpha$")
    plt.title("Pitch angle: {}, E={} MeV, Ri={} RJ, ai={}°".format(partype, Ep, Ri, ai))
    # plt.savefig('./plots/pitch_angle'+runname+'.png',bbox_inches="tight")
    plt.show()
