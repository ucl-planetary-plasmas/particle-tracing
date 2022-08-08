function example_trajectory

mdfile = '../jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
partype = 'S+';                    % proton
Ep = 25;                         % energy 100 MeV
Ri = 70;                           % initial equatorial distance 4*Rp
ai = 10;                          % initial pitch angle 30 degrees
timespec = [0,0,20,0];             % run for 2 dipole bounce periods
savefile = [];                    % name for result mat-file
pauseOn = false;                  % flag to pause on/ pause off.
npertc = 20;                      % number of Boris iterations per gyroperiod


mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,npertc);

