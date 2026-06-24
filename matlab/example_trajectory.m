function example_trajectory

mdfile = '../data/jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
partype = 'S+';                           % proton
Ep = 25;                                  % energy [MeV]
Ri = 70;                                  % initial equatorial distance [Rp]
ai = 10;                                  % initial equatorial pitch angle [deg]
timespec = [0,0,20,0];                    % run for 20 dipole bounce periods
savefile = [];                            % name for result mat-file
pauseOn = true;                           % flag to pause on/off.
plotOn = true;                            % flag to plot on/off.
npertc = 20;                              % number of iterations/equatorial gyroperiod


mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc);

