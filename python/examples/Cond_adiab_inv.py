######################
## Useful libraries ##
######################
import matplotlib.pyplot as plt  
import numpy as np 
import pickle as pkl 

import sys 
sys.path.append("../")
import os 

## To compute trajectories
from pymagdisc.data.load_data import load_model
from pymagdisc.tracer.mdbtracer import MDBTracer
from pymagdisc import config

#############################################
## Folder to save simulations already done ##
#############################################
dossier = "res_simul"
if not os.path.exists(dossier):
    os.makedirs(dossier)

Db,Df = {},{} ## Dictionnaries to save the data. It is reset to {} here, the creation of a .pkl might be more appropriate.

######################
## Useful functions ##
######################
# cosinus and sinus with an argument in degree
sindeg = lambda x: np.sin(x*np.pi/180)
cosdeg = lambda x: np.cos(x*np.pi/180)
# Function to smooth the data
def lissage(signal_brut,L):
    res = np.copy(signal_brut) # duplication des valeurs
    for i in range (1,len(signal_brut)-1): # toutes les valeurs sauf la première et la dernière
        L_g = min(i,L) # nombre de valeurs disponibles à gauche
        L_d = min(len(signal_brut)-i-1,L) # nombre de valeurs disponibles à droite
        Li=min(L_g,L_d)
        res[i]=np.sum(signal_brut[i-Li:i+Li+1])/(2*Li+1)
    return res
# Function to put the radii in int if necessary
# Used to havethe good names for the files
def round_radius(r):
    if r%1==0.0:
        return int(r)
    return r

#%%###############################
## Simulation of one trajectory ##
##################################

########################
## Initial parameters ##
########################

# magnetodisc mat-file (rmp60 for compressed, rmp90 for expanded models)
mdfile = "sat_mdisc_kh2e6_rmp25.mat"  
# Load the corresponding model
MD = load_model(f"{config.PATH_TO_DATA}{mdfile}")

## Creation of a folder to save data for this mdfile
directory = dossier + "/" + mdfile[:-4]
if not os.path.exists(directory):
    os.makedirs(directory)
    
partype = "p"               # particle type (select: p, e, O+, O++, S+, S++, S+++)
E = np.logspace(-3,.5,36)   # initial energy (in MeV)
R = np.linspace(5,25,41)    # initial equatorial distance (in Jovian radii or Kronian radii, depending on mdfile)
ai = 30                     # initial pitch angle (in degrees)
timespec = [0, 0, 3, 0]     # run for x dipole bounce periods (modify the third element)
npertc = 5                  # number of Boris iterations per gyroperiod (~ points used per bounce)
mdisctype = ("comp")        # Type 'exp' for expanded or 'comp' for compressed magnetosphere models

## New keys in the dictionnaries associated with this mdfile
name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(ai)
for D in [Db,Df]:
    D[name_for_dict] = {}

## Creation of lists of data to compute
Etab1,Etaf1,Etab2,Etaf2,Etab3,Etaf3 = [],[],[],[],[],[]

################################
## Simulation and computation ##
################################
for Ep in E:
    Etab11,Etaf11,Etab22,Etaf22,Etab33,Etaf33 = [],[],[],[],[],[]
    for Ri in R:
        print('Ep = ' + str(np.round(Ep,decimals=5)),'Ri = '+str(int(Ri))) # To check where we are
        
        ##################################
        ## Simulation of the trajectory ##
        ##################################
        # Name of the run, based on the parameters above
        runname = (
                partype + "_Ri" + str(round_radius(Ri)) + "_Ep" + str(np.round(Ep,decimals=5)) + "_ai" + str(int(ai)) + "_tb" + str(timespec[2])
        )
        # This creates the name of the file
        # Choose to save the results in a "pkl" fileù, adapted for python dictionaries
        savefile = runname + '.pkl'
        filepath = directory+"/"+savefile
        
        # Check if the file already exists
        # If yes, we load the data directly
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f:
                res = pkl.load(f)
        # If no, run the simulation of the trajectory
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
            with open(filepath, 'wb') as f: # Save the data to reuse it
                pkl.dump(res, f)
        
        #########################
        ## Useful computations ##
        #########################
        rgbn = np.sqrt(res['rgb'][:,0]**2+res['rgb'][:,1]**2+res['rgb'][:,2]**2)
        rgfn = np.sqrt(res['rgf'][:,0]**2+res['rgf'][:,1]**2+res['rgf'][:,2]**2)
        l = res['latb']
        
        ## Implementation of Young's paper's equations
        # Using rho*gradB/B
        eta = abs(2*np.pi*(rgbn/res['gradBb'][4])*cosdeg(res['aib']))
        etab1 = np.mean(np.sort(eta)[-30:])
        eta = abs(2*np.pi*(rgfn/res['gradBf'][4])*cosdeg(res['aif']))
        etaf1 = np.mean(np.sort(eta)[-30:])
        
        # Using rho/Rc
        eta = abs(2*np.pi*(rgbn/res['curvBb'][0])*cosdeg(res['aib']))
        etab2 = np.mean(np.sort(eta)[-30:]) 
        eta = abs(2*np.pi*(rgfn/res['curvBf'][0])*cosdeg(res['aif']))
        etaf2 = np.mean(np.sort(eta)[-30:]) 
        
        # Using the max between the two first expressions
        maxeps = np.maximum(rgbn/res['gradBb'][4],rgbn/res['curvBb'][0])
        maxeta = abs(2*np.pi*maxeps*cosdeg(res['aib']))
        etab3 = np.mean(np.sort(maxeta)[-30:])
        maxeps = np.maximum(rgfn/res['gradBf'][4],rgfn/res['curvBf'][0])
        maxeta = abs(2*np.pi*maxeps*cosdeg(res['aif']))
        etaf3 = np.mean(np.sort(maxeta)[-30:])
        
        ## Add in different lists
        Etab11.append(etab1)
        Etaf11.append(etaf1)
        Etab22.append(etab2)
        Etaf22.append(etaf2)
        Etab33.append(etab3)
        Etaf33.append(etaf3)
    Etab1.append(Etab11)
    Etaf1.append(Etaf11)
    Etab2.append(Etab22)
    Etaf2.append(Etaf22)
    Etab3.append(Etab33)
    Etaf3.append(Etaf33)
    
## Add data to dictionnaries
Db[name_for_dict]['Eta1'] = Etab1
Df[name_for_dict]['Eta1'] = Etaf1
Db[name_for_dict]['Eta2'] = Etab2
Df[name_for_dict]['Eta2'] = Etaf2
Db[name_for_dict]['Eta3'] = Etab3
Df[name_for_dict]['Eta3'] = Etaf3
    

#%%############
## Show data ##
###############
## It can be done using the data saved in the dictionnaries

## Comparison between a trajectory and RQA
## Studied trajectory
mdfile = "jup_mdisc_kh3e7_rmp90.mat"
ai = 30
name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(ai)

## Which parameters to compare it to
## The RQA measures must have been done previously
whichRQA = 'DIV'    # RR, ENT, DIV, LAM, TT, DET, ADL
whichparam = 2      # 0 for mu, 1 for lambda, 2 for alpha
Paraminit = [r'$\mu$',r'$\lambda$',r'$\alpha$'] # Names for labels
## With or without centrifugal force
CF = False

if CF:
    D = Df[name_for_dict]
    with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(ai)+'_CF/Data.pkl','rb') as f: ## Load the RQA data
        Data = pkl.load(f)
else:
    D = Db[name_for_dict]
    with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(ai)+'_tb10_m4t3_333/Data.pkl','rb') as f: ## Load the RQA data
        Data = pkl.load(f)
 

plt.rcParams.update({"text.usetex":True,'font.family':'Computer Modern','font.size': 20})
fig = plt.figure(figsize = (20,15), constrained_layout=False)

plt.pcolor(R,E,np.array(Data[whichRQA][whichparam])) ## RQA associated with the parameter chosen

#CS = plt.contour(R,E,D['Eta1'],levels=[.01,.02,.03,.05,.07,.1,.2,.5,1,2],colors='blue')
#plt.clabel(CS)
#CS = plt.contour(R,E,D['Eta2'],levels=[.01,.02,.03,.05,.07,.1,.2,.5,1,2],colors='red')
#plt.clabel(CS)
CS = plt.contour(R,E,D['Eta3'],levels=[.01,.02,.03,.05,.07,.1,.2,.5,1,2],colors='white')
plt.clabel(CS)
plt.vlines(100,1,2,label=r'isolines of $\eta$',color='white')

plt.xlabel(r'R [$\mathrm{R}_j$]'*(mdfile[:3]=='jup')+r'R [$\mathrm{R}_s$]'*(mdfile[:3]=='sat'))
plt.ylabel(r'E [MeV]')
plt.title(r'Exp. model of Jovian magn., '+whichRQA + r' of '+Paraminit[whichparam]+r', initial pitch angle $\alpha_i =$'+str(ai)
          +r', centrifugal force added'*CF,pad=20)
plt.xlim(min(R),max(R))
plt.yscale('log')
plt.legend(loc='upper right')
plt.show()

#%% To plot different RQA at the same time

fig = plt.figure(figsize = (20,15), constrained_layout=False)
gs = fig.add_gridspec(2,2,wspace=0.1,hspace=0.2)

mdfile = "jup_mdisc_kh3e7_rmp90.mat"
a = 30
name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(a)
CF = False
Paraminit = [r'$\mu$',r'$\lambda$',r'$\alpha$']
W = ['RR','ENT','DIV','ENT'] # Quantificators of the RQA
P = [0,0,2,2]                # Parameter of the motion associated (0 for mu, 1 for lambda, 2 for alpha)

for i in range(len(W)):
    w = W[i]
    p = P[i]
    if CF: 
        D = Df[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(ai)+'_CF/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    else:
        D = Db[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_tb10_m4t3_333/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
            
    f = fig.add_subplot(gs[i//2,i%2])
    plt.pcolor(R,E,np.array(Data[w][p]))
    plt.colorbar()
    #CS = plt.contour(R,E,D['Eta3'],levels=[.02,.05,.1,.2,.5,1],colors='black')
    #plt.clabel(CS)
    if i//2==1:
        plt.xlabel(r'$R_{\mathrm{eq}}$ [$\mathrm{R}_j$]'*(mdfile[:3]=='jup')+r'$R_{eq}}$ [$\mathrm{R}_s$]'*(mdfile[:3]=='sat'))
    if i%2==0:
        plt.ylabel(r'E [MeV]')
    plt.yscale('log')
    plt.title(w + r' for '+Paraminit[p]+r' CF'*CF,pad=10)
plt.show()

#%% To plot the four pitch angle together
mdfile = "jup_mdisc_kh3e7_rmp90.mat"

## Which parameters to compare it to
whichRQA = 'DIV' # RR, ENT, DIV, LAM, TT, DET, ADL
whichparam = 2  # 0 for mu, 1 for lambda, 2 for alpha
Paraminit = [r'$\mu$',r'$\lambda$',r'$\alpha$']

fig = plt.figure(figsize = (20,15), constrained_layout=False)
gs = fig.add_gridspec(2,2,wspace=0.15,hspace=0.25)
A = [10,30,50,70]
CF = False

for i in range(len(A)):
    a = A[i]
    name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(a)
    if CF:
        D = Df[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_CF/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    else:
        D = Db[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_tb10_m4t3_333/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    f = fig.add_subplot(gs[i//2,i%2])
    plt.pcolor(R,E,np.array(Data[whichRQA][whichparam]))
    #plt.colorbar()
    #CS = plt.contour(R,E,D['Eta3'],levels=[.01,.02,.05,.1,.2,.5,1],colors='white')
    #plt.clabel(CS)
    if i//2==1:
        plt.xlabel(r'$R_{\mathrm{eq}}$ [$\mathrm{R}_j$]'*(mdfile[:3]=='jup')+r'R [$\mathrm{R}_s$]'*(mdfile[:3]=='sat'))
    if i%2==0:
        plt.ylabel(r'E [MeV]')
    plt.yscale('log')
    plt.title(r'$\alpha_i =$'+str(a)+r' CF'*CF,pad=10)
plt.show()

#%% To plot different mdfiles
fig = plt.figure(figsize = (20,15), constrained_layout=False)
gs = fig.add_gridspec(1,2,wspace=0.3)

## Which parameters to compare it to
whichRQA = 'RR' # RR, ENT, DIV, LAM, TT, DET, ADL
whichparam = 0  # 0 for mu, 1 for lambda, 2 for alpha
Paraminit = [r'$\mu$',r'$\lambda$',r'$\alpha$']

Mdfiles = ["jup_mdisc_kh3e7_rmp90.mat","jup_mdisc_kh3e7_rmp60.mat","sat_mdisc_kh2e6_rmp25.mat"]
Eng = [np.logspace(-2,.9,30),np.logspace(-2,1.8,39),np.logspace(-3,.5,36)]
Ray = [np.linspace(10,70,61),np.linspace(10,60,51),np.linspace(5,25,41)]
a = 30
CF = False

for i in range(1,len(Mdfiles)):
    mdfile = Mdfiles[i]
    name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(a)
    E = Eng[i]
    R = Ray[i]
    if CF:
        D = Df[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_CF/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    else:
        D = Db[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_tb10_m4t3_333/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    f = fig.add_subplot(gs[0,i-1])
    #plt.gca().set_aspect(5)
    plt.pcolor(R,E,np.array(Data[whichRQA][whichparam]))
    #plt.colorbar(shrink=0.6)
    CS = plt.contour(R,E,D['Eta3'],levels=[.01,.02,.05,.1,.2,.5,1],colors='white')
    plt.clabel(CS)
    plt.vlines(100,1,2,label=r'isolines of $\eta$',color='white')
    plt.xlim(min(R),max(R))
    plt.xlabel(r'$R_{\mathrm{eq}}$ [$\mathrm{R}_j$]'*(mdfile[:3]=='jup')+r'$R_{\mathrm{eq}}$ [$\mathrm{R}_s$]'*(mdfile[:3]=='sat'))
    if i==1:
        plt.ylabel(r'E [MeV]')
    plt.yscale('log')
    plt.title(mdfile[:3]+mdfile[19:21], pad=10)#+r', $\alpha_i =$'+str(a)+r' CF'*CF,pad=10)
    plt.legend(loc = 'lower right')
    
#%% To plot different mdfiles
fig = plt.figure(figsize = (20,15), constrained_layout=False)
gs = fig.add_gridspec(1,2,wspace=0.3)

## Which parameters to compare it to
whichRQA = 'RR' # RR, ENT, DIV, LAM, TT, DET, ADL
whichparam = 0  # 0 for mu, 1 for lambda, 2 for alpha
Paraminit = [r'$\mu$',r'$\lambda$',r'$\alpha$']

mdfile = "jup_mdisc_kh3e7_rmp90.mat"
E = np.logspace(-2,.9,30)
R = np.linspace(10,70,61)
a = 30
name_for_dict = mdfile[:3]+mdfile[19:21]+'_a'+str(a)

for CF in [0,1]:
    if CF:
        D = Df[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_CF/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    else:
        D = Db[name_for_dict]
        with open('Data_RecPlot2D/'+mdfile[:3]+mdfile[12:15]+mdfile[19:21]+'_p_a'+str(a)+'_tb10_m4t3_333/Data.pkl','rb') as f: ## Load the RQA data
            Data = pkl.load(f)
    f = fig.add_subplot(gs[0,CF])
    #plt.gca().set_aspect(5)
    plt.pcolor(R,E,np.array(Data[whichRQA][whichparam]))
    #plt.colorbar(shrink=0.6)
    CS = plt.contour(R,E,D['Eta3'],levels=[.01,.02,.05,.1,.2,.5,1],colors='white')
    plt.clabel(CS)
    plt.vlines(100,1,2,label=r'isolines of $\eta$',color='white')
    plt.xlim(min(R),max(R))
    plt.xlabel(r'$R_{\mathrm{eq}}$ [$\mathrm{R}_j$]'*(mdfile[:3]=='jup')+r'$R_{\mathrm{eq}}$ [$\mathrm{R}_s$]'*(mdfile[:3]=='sat'))
    if i==1:
        plt.ylabel(r'E [MeV]')
    plt.yscale('log')
    plt.title(r'With '*CF+r'Without '*(1-CF)+r'centrifugal force', pad=10)#+r', $\alpha_i =$'+str(a)+r' CF'*CF,pad=10)
    plt.legend(loc = 'lower right')

