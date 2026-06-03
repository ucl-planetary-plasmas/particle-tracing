################################
### Save the simulation data ###
################################

import os
os.sys.path.append('../')

import pickle as pkl
import numpy as np
from pymagdisc.tracer.mdbtracer import MDBTracer

## Folder where to save the results
directory = "res_simul"
if not os.path.exists(directory):
    os.makedirs(directory)

if __name__ == "__main__":
    ### Change the following set of parameters

    # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
    mdfile = "sat_mdisc_kh2e6_rmp25.mat"
    # particle type (select: p, e, O+, O++, S+, S++, S+++)
    partype = "p"
    # initial energy (in MeV)
    E = np.logspace(0,1,11)
    # initial equatorial distance (in Jovian radii)
    Ri = 30
    # initial pitch angle (in degrees)
    ai = 30
    # run for x dipole bounce periods (modify the third element)
    timespec = [0, 0, 10, 0]
    # number of Boris iterations per gyroperiod (~ points used per bounce)
    npertc = 15
    # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
    mdisctype = "exp"

    # --- Do not modify the following unless necessary! ---
    
    for Ep in E:
        print(Ep)
        # The code block below creates a text ("string" in python terminology), based on the parameters above
        runname = (
            partype + "_Ri" + str(Ri) + "_Ep" + str(np.round(Ep,decimals=2)) + "_ai" + str(ai) + "_" + mdfile[0:3] + mdfile[12:15] + mdfile[19:21] + mdisctype
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

#%%
with open(filepath, 'rb') as f:
    data = pkl.load(f)
