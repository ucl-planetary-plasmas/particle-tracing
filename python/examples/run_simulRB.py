################################
### Save the simulation data ###
################################
import os
os.sys.path.append('../')

import pickle as pkl
import numpy as np
from pymagdisc.tracer.mdbtracer import MDBTracer

## Folder where to save the results
dossier = "res_simul"
if not os.path.exists(dossier):
    os.makedirs(dossier)

if __name__ == "__main__":
    ### Change the following set of parameters

    # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
    MDfile = ["sat_mdisc_kh2e6_rmp25.mat","jup_mdisc_kh3e7_rmp60.mat","jup_mdisc_kh3e7_rmp90.mat"]
    # particle type (select: p, e, O+, O++, S+, S++, S+++)
    partype = "p"
    # initial energy (in MeV)
    E = np.logspace(-3,.5,36)
    # initial equatorial distance (in Jovian radii)
    R = np.linspace(5,70,66)
    # initial pitch angle (in degrees)
    a = np.linspace(10,80,15)
    # run for x dipole bounce periods (modify the third element)
    timespec = [0, 0, 10, 0]
    # number of Boris iterations per gyroperiod (~ points used per bounce)
    npertc = 15
    # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
    mdisctype_list = ["exp","comp","exp"]

    for i in range(len(MDfile)):
        mdfile = MDfile[i]
        mdisctype = mdisctype_list[i]
        directory = dossier + "/" + mdfile[:-4]
        if not os.path.exists(directory):
            os.makedirs(directory)
        for Ep in E:
            for Ri in R:
                for ai in a:
                    # The code block below creates a text ("string" in python terminology), based on the parameters above
                    runname = (
                        partype + "_Ri" + str(Ri) + "_Ep" + str(np.round(Ep,decimals=2)) + "_ai" + str(ai)
                    )
                
                    # This creates the name of the file
                    # Choose to save the results in a "mat" file
                    savefile = runname + '.pkl'
                    filepath = directory+"/"+savefile
                    
                    if os.path.exists(filepath):
                        print("File already exists")
                    else:
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
                        res = tracer_jup.run_simul()
                        ## Save the file
                        with open(filepath, 'wb') as f:
                            pkl.dump(res, f)

