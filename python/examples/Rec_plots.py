########################
### Recurrence plots ###
########################
"""The idea is to plot a function:
    R(i,j) = 1 if |x(i)-x(j)| < e and 0 if not
   that represents when their is a recurrence in a function between the times i and j."""
   

import matplotlib.pyplot as plt
import numpy as np
import pickle as pkl
import sys
import os

sys.path.append("../")
from pymagdisc.tracer.mdbtracerRB import MDBTracer
from pyunicorn.timeseries import RecurrencePlot as RecPlot


#%% Function to smooth the data
def lissage(signal_brut,L):
    res = np.copy(signal_brut) # duplication des valeurs
    for i in range (1,len(signal_brut)-1): # toutes les valeurs sauf la première et la dernière
        L_g = min(i,L) # nombre de valeurs disponibles à gauche
        L_d = min(len(signal_brut)-i-1,L) # nombre de valeurs disponibles à droite
        Li=min(L_g,L_d)
        res[i]=np.sum(signal_brut[i-Li:i+Li+1])/(2*Li+1)
    return res


#%%
## The different measures we make with recurrence plots
RR_mu,DET_mu,LAM_mu,ENT_mu,DIV_mu,ADL_mu,TT_mu = [],[],[],[],[],[],[]
RR_lat,DET_lat,LAM_lat,ENT_lat,DIV_lat,ADL_lat,TT_lat = [],[],[],[],[],[],[]
RR_a,DET_a,LAM_a,ENT_a,DIV_a,ADL_a,TT_a = [],[],[],[],[],[],[]

### Change the following set of parameters
mdfile = "sat_mdisc_kh2e6_rmp25.mat"  # magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
partype = "p"
# particle type (select: p, e, O+, O++, S+, S++, S+++)
E = np.logspace(0,1,11)
# initial energy (in MeV)
Ri = 30
# initial equatorial distance (in Jovian radii or Kronian radii, depending on mdfile)
ai = 30
# initial pitch angle (in degrees)
timespec = [0, 0, 10, 0]
# run for x dipole bounce periods (modify the third element)
npertc = 15  # number of Boris iterations per gyroperiod (~ points used per bounce)
mdisctype = (
    "exp"  # Type 'exp' for expanded or 'comp' for compressed magnetosphere models
)

#%%
for Ep in E:
    # --- Do not modify the following unless necessary! ---
    # The code block below creates a text ("string" in python terminology), based on the parameters above
    runname = (
        partype + "_Ri" + str(Ri) + "_Ep" + str(np.round(Ep,decimals=2)) + "_ai" + str(ai) + "_" + mdfile[0:3] + mdfile[12:15] + mdfile[19:21] + mdisctype
    )
    # This creates the name of the file
    # Choose to save the results in a "mat" file
    savefile = runname + '.pkl'
    filepath = "res_simul/"+savefile
    
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
    divi = int(len(results['muib'])/3000)
    divi = divi + 1*(divi==0)
    ## Resampling
    Mu = results['muib'][::divi]
    Lat = results['latb'][::divi]
    Alpha = results['aib'][::divi]
    T = results['tb'][::divi]
   
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
    
## Gathering data in one list
Label=['RR','ADL','DIV','ENT','LAM','DET','TT']
Data = [(RR_mu,RR_lat,RR_a),
        (ADL_mu,ADL_lat,ADL_a),
        (DIV_mu,DIV_lat,DIV_a),
        (ENT_mu,ENT_lat,ENT_a),
        (LAM_mu,LAM_lat,LAM_a),
        (DET_mu,DET_lat,DET_a),
        (TT_mu,TT_lat,TT_a)]

#%%
## If you want to save
## Use the format : planet_kh_rmp_particle_Ri_ai_nbbounceperiod_parametersembeddeding_precisionforeachdata
## Example : 'sat2e625_p_20R_30a_10tb_4m3t_355/'
def save_plots(namefile,e=0,r=0):
    # Create output directory if needed
    dos = "../../../../Figures/"
    os.mkdir(dos+namefile)
    for j in range(7):
        ## Mu
        plt.figure(figsize = (20,15), constrained_layout=False)
        plt.scatter(E,Data[j][0],color='black',marker='+')
        if e:
            plt.semilogx(E,lissage(Data[j][0],5),color='red',label=r'$\mu$')
        if r:
            plt.plot(E,lissage(Data[j][0],5),color='red',label=r'$\mu$')
        plt.legend(loc = 'upper left')
        plt.xlabel(r'E [MeV]'*e+r'R [$R_j$]'*r)
        plt.ylabel(Label[j])
        plt.grid()
        plt.savefig(dos+namefile+Label[j]+'_mu',format='pdf')
        
        ## Lat and Alpha
        plt.figure(figsize = (20,15), constrained_layout=False)
        plt.scatter(E,Data[j][1],color='red',label=r'$\lambda$',marker='+')
        plt.scatter(E,Data[j][2],color='blue',label=r'$\alpha$',marker='+')
        if e: 
            plt.semilogx(E,lissage(Data[j][1],5),color='orange',label=r'$\lambda$')
            plt.semilogx(E,lissage(Data[j][2],5),color='green',label=r'$\alpha$')
        if r: 
            plt.plot(E,lissage(Data[j][1],5),color='orange',label=r'$\lambda$')
            plt.plot(E,lissage(Data[j][2],5),color='green',label=r'$\alpha$')
        plt.legend(loc = 'upper left')
        plt.xlabel(r'E [MeV]'*e+r'R [$R_j$]'*r)
        plt.ylabel(Label[j])
        plt.grid()
        plt.savefig(dos+namefile+Label[j]+'_latalpha',format='pdf')
        
        ## Normalized and lissed
        plt.figure(figsize = (20,15), constrained_layout=False)
        if e:
            plt.semilogx(E,lissage(Data[j][0],5)/max(Data[j][0]),color='black',ls='-.',label=r'$\mu$')
            plt.semilogx(E,lissage(Data[j][1],5)/max(Data[j][1]),color='blue',ls='-.',label=r'$\lambda$')
            plt.semilogx(E,lissage(Data[j][2],5)/max(Data[j][2]),color='red',ls='-.',label=r'$\alpha$')
        if r:
            plt.plot(E,lissage(Data[j][0],5)/max(Data[j][0]),color='black',ls='-.',label=r'$\mu$')
            plt.plot(E,lissage(Data[j][1],5)/max(Data[j][1]),color='blue',ls='-.',label=r'$\lambda$')
            plt.plot(E,lissage(Data[j][2],5)/max(Data[j][2]),color='red',ls='-.',label=r'$\alpha$')
        plt.legend(loc = 'upper left')
        plt.xlabel(r'E [MeV]'*e+r'R [$R_j$]'*r)
        plt.ylabel(Label[j])
        plt.grid()
        plt.savefig(dos+namefile+Label[j],format='pdf')
        plt.close('all')
    np.savez(dos+namefile+'Data',Data,Label)

save_plots('jup3e790_p_10E_30a_10tb_7m7t_333/',r=1)

#%% Tests
j=0
mu = 0
latalpha = 1
Norm = 0
fig4 = plt.figure(figsize = (20,15), constrained_layout=False)
if mu:
    plt.semilogx(E,Data[0][0],ls = '-.',color='black',label='Mu')
    plt.scatter(E,Data[j][0],color='black',marker='+')
    plt.semilogx(E,lissage(Data[j][0],5),color='red',label=r'$\mu$')
elif latalpha:
    #plt.semilogx(E,Data[j][1],ls = '-.',color='red',label=r'$\lambda$')
    plt.scatter(E,Data[j][1],color='red',label=r'$\lambda$',marker='+')
    plt.semilogx(E,lissage(Data[j][1],5),color='orange',label=r'$\lambda$')

    #plt.semilogx(E,Data[j][2],ls = '-.',color='blue',label=r'$\alpha$')
    plt.scatter(E,Data[j][2],color='blue',label=r'$\alpha$',marker='+')
    plt.semilogx(E,lissage(Data[j][2],5),color='green',label=r'$\alpha$')  
elif Norm:
    plt.semilogx(E,lissage(Data[j][0],5)/max(Data[j][0]),color='black',ls='-.',label=r'$\mu$')
    plt.semilogx(E,lissage(Data[j][1],5)/max(Data[j][1]),color='blue',ls='-.',label=r'$\lambda$')
    plt.semilogx(E,lissage(Data[j][2],5)/max(Data[j][2]),color='red',ls='-.',label=r'$\alpha$')

plt.legend(loc = 'upper left')
plt.xlabel('E [MeV]')
plt.ylabel(Label[j])
plt.grid()
plt.show()




   
