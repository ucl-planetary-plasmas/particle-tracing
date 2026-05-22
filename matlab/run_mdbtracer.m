function run_mdbtracer

mdfile = '../data/jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
partype = 'p';                    % proton
Ep = 100;                         % energy 100 MeV
Ri = 4;                           % initial equatorial distance 4*Rp
ai = 30;                          % initial pitch angle 30 degrees
timespec = [0,0,2,0];             % run for 2 dipole bounce periods
savefile = 'my_jup_mdisc_sim';    % name for result mat-file
pauseOn = false;                  % flag to pause on/ pause off.
plotOn = false;                   % flag to plot on/ plot off.
npertc = 5;                       % number of Boris iterations per gyroperiod


mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc);

