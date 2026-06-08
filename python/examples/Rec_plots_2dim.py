########################
### Recurrence plots ###
########################   
import os
os.sys.path.append('../')

import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

from pymagdisc.tracer.mdbtracer import MDBTracer
from pyunicorn.timeseries import RecurrencePlot as RecPlot

##############################################
### Change the following set of parameters ###
##############################################
mdfile = "jup_mdisc_kh3e7_rmp60.mat"  # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
partype = "p"
# particle type (select: p, e, O+, O++, S+, S++, S+++)
E = np.logspace(-2,1,31)
# initial energy (in MeV)
R = np.linspace(15,70,51)
# initial equatorial distance (in Jovian radii or Kronian radii, depending on mdfile)
ai = 30
# initial pitch angle (in degrees)
timespec = [0, 0, 10, 0]
# run for x dipole bounce periods (modify the third element)
npertc = 15  # number of Boris iterations per gyroperiod (~ points used per bounce)
mdisctype = (
    "exp"  # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
)

## The different measures we make with recurrence plots
RR_mu_tot,DET_mu_tot,LAM_mu_tot,ENT_mu_tot,DIV_mu_tot,ADL_mu_tot,TT_mu_tot = [],[],[],[],[],[],[]
RR_lat_tot,DET_lat_tot,LAM_lat_tot,ENT_lat_tot,DIV_lat_tot,ADL_lat_tot,TT_lat_tot = [],[],[],[],[],[],[]
RR_a_tot,DET_a_tot,LAM_a_tot,ENT_a_tot,DIV_a_tot,ADL_a_tot,TT_a_tot = [],[],[],[],[],[],[]


i,j=0,0
for Ep in E:
    i+=1
    j=0
    RR_mu,DET_mu,LAM_mu,ENT_mu,DIV_mu,ADL_mu,TT_mu = [],[],[],[],[],[],[]
    RR_lat,DET_lat,LAM_lat,ENT_lat,DIV_lat,ADL_lat,TT_lat = [],[],[],[],[],[],[]
    RR_a,DET_a,LAM_a,ENT_a,DIV_a,ADL_a,TT_a = [],[],[],[],[],[],[]
    for Ri in R:
        j+=1
        print(i,j)
        # --- Do not modify the following unless necessary! ---
        # The code block below creates a text ("string" in python terminology), based on the parameters above
        runname = (
            partype + "_Ri" + str(int(Ri)) + "_Ep" + str(np.round(Ep,decimals=5)) + "_ai" + str(int(ai))
        )
        # This creates the name of the file
        # Choose to save the results in a "mat" file
        savefile = runname + '.pkl'
        directory = "res_simul/" + mdfile
        filepath = directory+"/"+savefile
        
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f:
                results = pkl.load(f)
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
            results = tracer_jup.run_simul()
            ## Save the file
            with open(filepath, 'wb') as f:
                pkl.dump(results, f)
                 
        ## Factor for resampling data to have approximately the same number of points
        nb_pts = 3000
        divi = int(len(results['muib'])/nb_pts)
        divi = divi + 1*(divi==0) ## In case the data doesn't have nb_pts points
        ## Resampling
        Mu = results['muib'][::divi]
        Lat = results['latb'][::divi]
        Alpha = results['aib'][::divi]
       
        ## Parameters for the recurrence plots
        m,t=4,3 ## Embedded dimensions and time delay
        # Precision for each plot
        trh_mu = 0.03*(max(abs(results['muib']))-np.mean(results['muib'])) 
        trh_lat = 0.03*max(abs(results['latb']))
        trh_a = 0.03*max(abs(results['aib']))
        # Recurrence plots
        rp_mu = RecPlot(Mu,dim=m,tau=t,threshold=trh_mu)
        rp_lat = RecPlot(Lat,dim=m,tau=t,threshold=trh_lat)
        rp_a = RecPlot(Alpha,dim=m,tau=t,threshold=trh_a)
        
        ## RQA
        # Recurrence Rate
        RR_mu.append(rp_mu.recurrence_rate())
        RR_lat.append(rp_lat.recurrence_rate())
        RR_a.append(rp_a.recurrence_rate())
    
        # Determinism
        DET_mu.append(rp_mu.determinism())
        DET_lat.append(rp_lat.determinism())
        DET_a.append(rp_a.determinism())
        
        # Laminarity
        LAM_mu.append(rp_mu.laminarity())
        LAM_lat.append(rp_lat.laminarity())
        LAM_a.append(rp_a.laminarity())
    
        # Entropy
        ENT_lat.append(rp_lat.white_vert_entropy())
        ENT_mu.append(rp_mu.white_vert_entropy())
        ENT_a.append(rp_a.white_vert_entropy())
        
        # Average diagonal length
        ADL_lat.append(rp_lat.average_diaglength())
        ADL_mu.append(rp_mu.average_diaglength())
        ADL_a.append(rp_a.average_diaglength())
    
        # Divergence
        DIV_lat.append(1/rp_lat.max_diaglength())
        DIV_mu.append(1/rp_mu.max_diaglength())
        DIV_a.append(1/rp_a.max_diaglength())
    
        # Traping Time
        TT_mu.append(rp_mu.trapping_time())
        TT_lat.append(rp_lat.trapping_time())
        TT_a.append(rp_a.trapping_time())
    RR_mu_tot.append(RR_mu)
    RR_lat_tot.append(RR_lat)
    RR_a_tot.append(RR_a)
    DET_mu_tot.append(DET_mu)
    DET_lat_tot.append(DET_lat)
    DET_a_tot.append(DET_a)
    LAM_mu_tot.append(LAM_mu)
    LAM_lat_tot.append(LAM_lat)
    LAM_a_tot.append(LAM_a)
    ENT_mu_tot.append(ENT_mu)
    ENT_lat_tot.append(ENT_lat)
    ENT_a_tot.append(ENT_a)
    ADL_mu_tot.append(ADL_mu)
    ADL_lat_tot.append(ADL_lat)
    ADL_a_tot.append(ADL_a)
    DIV_mu_tot.append(DIV_mu)
    DIV_lat_tot.append(DIV_lat)
    DIV_a_tot.append(DIV_a)
    TT_mu_tot.append(TT_mu)
    TT_lat_tot.append(TT_lat)
    TT_a_tot.append(TT_a)
 
## Gathering data in one list
Data={'RR':(RR_mu_tot,RR_lat_tot,RR_a_tot),
      'DET':(DET_mu_tot,DET_lat_tot,DET_a_tot),
      'LAM':(LAM_mu_tot,LAM_lat_tot,LAM_a_tot),
      'ENT':(ENT_mu_tot,ENT_lat_tot,ENT_a_tot),
      'ADL':(ADL_mu_tot,ADL_lat_tot,ADL_a_tot),
      'DIV':(DIV_mu_tot,DIV_lat_tot,DIV_a_tot),
      'TT':(TT_mu_tot,TT_lat_tot,TT_a_tot),
      }

np.savez("../../../../Figures/2D/sat2e625_p_30a_10tb_4m3t_333/Data",Data)  

   
for Param in Data.keys():
    Names = [r'$\mu$',r'$\lambda$',r'$\alpha$']
    plt.rcParams.update({"text.usetex":True,'font.family':'Computer Modern','font.size': 20})
    fig = plt.figure(figsize = (20,15), constrained_layout=False)
    gs = fig.add_gridspec(3,1,hspace=0.35)
    
    for i in range(3):
        f_ax = fig.add_subplot(gs[i,:]) ## subplot
        plt.pcolor(R,E,Data[Param][i])
        if i == 2:
            plt.xlabel(r'R [$\mathrm(R)_j$]')
        plt.ylabel(r'E [MeV]')
        plt.title(Names[i],pad=15)
        plt.xlim(5,30)
        plt.yscale('log')
        plt.colorbar()
        plt.savefig("../../../../Figures/2D/sat2e625_p_30a_10tb_4m3t_333/"+Param,format='pdf')
        plt.show()