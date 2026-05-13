
from pathlib import Path

import numpy as np
from scipy.io import savemat
import matplotlib.pyplot as plt

from pymagdisc import config
from pymagdisc.data.load_data import load_model
from pymagdisc.tracer.mdbtracer import MDBTracer


if __name__ == "__main__":

    ### Change the following set of parameters

    # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
    mdfile = "jup_mdisc_kh3e7_rmp90.mat"  
    # particle type (select: p, e, O+, O++, S+, S++, S+++)
    partype = "S+"
    # initial energy (in MeV)
    Ep = 25
    # initial equatorial distance (in Jovian radii)
    Ri = 70
    # initial pitch angle (in degrees)
    ai = 10
    # run for x dipole bounce periods (modify the third element)
    timespec = [0, 0, 20, 0]
    # number of Boris iterations per gyroperiod (~ points used per bounce)
    npertc = 20  
    # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
    mdisctype = "exp"

    # --- Do not modify the following unless necessary! ---

    # The code block below creates a text ("string" in python terminology), based on the parameters above
    runname = partype + "Ri" + str(Ri) + "_Ep" + str(Ep) + "_ai" + str(ai) + "" + mdisctype

    # This creates the name of the file
    # Choose to save the results in a "mat" file
    savefile = runname + ".mat"

    # Load the magnetic field model, run the code and save the results in the file specified above
    md_path = Path(config.PATH_TO_DATA) / mdfile
    MD = load_model(md_path)

    tracer = MDBTracer(
        mdfile=mdfile,
        partype=partype,
        Ep=Ep,
        Ri=Ri,
        ai=ai,
        timespec=timespec,
        savefile=savefile,
        npertc=npertc)

    res = tracer.run_simul()

    # Create output directory if needed
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)

    # Save output data
    output_path = results_dir / f"{savefile}.mat"
    savemat(output_path,res)

    # Plot the trajectory of the particle
    fig = plt.figure(figsize=(4, 4), dpi=120)
    plt.plot(res["Rcylb"] / res["Re"], res["Zb"] / res["Re"], "-")
    plt.xlabel(r"rc ($R_J$)")
    plt.ylabel(r"z ($R_J$)")
    plt.xlim(min(res["Zb"] / res["Re"] + 2), max(res["Zb"] / res["Re"] + 2))

    plt.title(r"Trajectory: {}, E={} MeV, Ri={} RJ, ai={}°".format(partype, Ep, Ri, ai))

    circle1 = plt.Circle((0, 0), 1, color="k")
    plt.gca().add_patch(circle1)
    plt.axis("equal")
    # plt.savefig('./plots/trajectory_'+runname+'.png',bbox_inches="tight")
    plt.show()

    for i in range(0, len(res["tb"])):
        if (
            np.sqrt((res["Zb"][i])** 2 + (res["Rcylb"][i])** 2) / res["Re"]
            <= 1.0
        ):
            print("Particle lost !")

    fig = plt.figure(figsize=(4, 4), dpi=120)
    plt.plot(res["tb"] / 60, res["pitchan"])
    plt.xlabel("Time (min)")
    plt.ylabel(r"$\alpha$")
    plt.title("Pitch angle: {}, E={} MeV, Ri={} RJ, ai={}°".format(partype, Ep, Ri, ai))
    # plt.savefig('./plots/pitch_angle'+runname+'.png',bbox_inches="tight")
    plt.show()

