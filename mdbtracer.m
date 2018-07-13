function mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,npertc)
% function mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,npertc)
%
% Magnetodisc Boris particle tracer
%
% mdfile   : magnetodisc mat-file filename.
% partype  : particle type 'p' -> proton, 'e' -> electron
% Ep       : particle energy in MeV.
% Ri       : initial particle equatorial position (in planet's radius Rp).
% ai       : initial particle pitch angle (0..180 degrees).
% timespec : defined as tmax = sum(timespec.*[1,tc,tb,td])
%            where 1 is in units of seconds, tc in units of gyroperiod at Ri,
%            tb in units of bounce period and td in units of drift period.
% savefile : filename to save the simulation data,
%            if not given, do not save simulation data.
% pauseOn  : flag to pause on/ pause off.
%            if not given, pause on.
% npertc   : optional number of Boris iterations per gyroperiod (default 5)
%
% For instance 
%
%   mdfile = 'jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
%   partype = 'p';                    % proton
%   Ep = 100;                         % energy 100 MeV
%   Ri = 4;                           % initial equatorial distance 4*Rp
%   ai = 30;                          % initial pitch angle 30 degrees
%   timespec = [0,0,2,0];             % run for 2 dipole bounce periods
%   savefile = 'my_jup_mdisc_sim';    % name for result mat-file
%
%   mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile);
%
% runs the Boris algorithm, and guiding centre approximation in a Jupiter
% magnetodisc for a proton of 100 MeV at initial position 4*Rp and
% pitch angle 30 deg, for 2*tb, i.e. two dipole bounce periods, and save the
% simulation data in file 'my_jup_mdisc_sim.mat'.
%
% Similarly for Saturn
%
%   mdfile = 'sat_mdisc_kh2e6_rmp25'; % magnetodisc mat-file
%   partype = 'p';                    % proton
%   Ep = 100;                         % energy 100 MeV
%   Ri = 4;                           % initial equatorial distance 4*Rp
%   ai = 30;                          % initial pitch angle 30 degrees
%   timespec = [0,0,2,0];             % run for 2 dipole bounce periods
%
%   mdbtracer(mdfile,partype,Ep,Ri,ai,timespec);
%
% but don't save the simulation data.

%
% $Id: mdbtracer.m,v 1.1 2018/07/13 16:50:52 patrick Exp $
%
% Copyright (c) 2018 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% close all figures
close all

% planet magnetodisc global filename
global MDFile
MDFile = mdfile;
% clear function with persistent variables
clear MDiscField
load(MDFile,'MD')

% Planet's equatorial radius in m
Re = MD.planet.r0;
% magnetic equator field strength in T
B0 = MD.planet.B0;
% corrected to match Bthdip(mu=0,r=1)
B0 = B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
Md = [0,0,B0*Re^3];  % Magnetic moment \mu_0 M/(4\pi) in T m^3
Rm = [0,0,0];        % Centered moment

% SI constants
e  = 1.60217733e-19; % elementary charge C
eV = 1.60217733e-19; % conversion factor from eV to J
c  = 299792458.0;    % speed of light m s^-1

switch partype,
  case 'p', % proton
    mp = 1.6726231e-27;  % proton mass kg
    qp = e;
  case 'e', % electron
    mp = 9.1093897e-31;  % electron mass kg
    qp = -e;
end

EMeV = Ep*1e6*eV;      % particle energy in MeV
vp = sqrt(2/mp*EMeV);  % classic velocity 
fprintf('classic v=%.4g m/s E=%.4g MeV\n',vp,.5*mp*vp^2/(1e6*eV))
vp = c*sqrt(EMeV*(EMeV+2*mp*c^2))/(EMeV+mp*c^2); % relativistic velocity
fprintf('relativ v=%.4g m/s E=%.4g MeV\n',vp,(1/sqrt(1-vp^2/c^2)-1)*mp*c^2/(1e6*eV));
% velocity squared
V2 = vp^2;

% relativistic factor
gamma = 1/sqrt(1-V2/c^2);

% particle charge-to-mass ratio
qOverM = qp/mp;

close all
% initial conditions
[X,mu,Tc,Tb,Td,Lm] = init(Ri*Re,vp,ai);

% factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
facdv = mu/gamma^2/mp; 
tspan = [0,sum(timespec.*[1,Tc,Tb,Td])];
fprintf(1,'tmax=%.2f (%.1f tc, %.2f tb, %.2f td)\n',...
        tspan(end)./[1,Tc,Tb,Td]);

if ~exist('savefile','var'), % define savefile if necessary
  savefile = [];
end

if ~exist('pauseOn','var'), % define pauseOn to true by default
  pauseOn = true;
end

if ~exist('npertc','var'), % define npertc=5 by default
  npertc = 5;
end

trace(X,mdfile,Ep,Ri,ai,timespec,Tc,Tb,Td,Lm,savefile,pauseOn,npertc);

function [Xo,mu,tc,tb,td,lm]=init(R,v,alpha)

% initial condition x0,y0,z0,vx0,vy0,vz0
Xo = [R;            0;            0;...
      0;v*sind(alpha);v*cosd(alpha)];

% magnetic field at initial position
B = mdiscMagneticField3D({Xo(1)/Re,Xo(2)/Re,Xo(3)/Re});
b = sqrt(B{4});

fprintf(1,'R=%.2f pitch angle=%.0f B=%.5g nT\n', R/Re, alpha, 1e9*b)

% v parallel and perpendicular
vpar = sum(Xo(4:6).*[B{1};B{2};B{3}]/b);
vper = sqrt(V2-vpar^2);

% gyro radius
rho = abs(gamma/qOverM*vper/b);

% first invariant mu
mu = gamma^2*mp*vper^2/(2*b);

[tc,tb,td] = periods(R,alpha);

% dipole mirror point latitude
lm = mirrorlat(alpha);

L = R/Re;
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
alphaloss = asind(sqrt(1/sqrt(4*L^6-3*L^5)));
fprintf(1,'|alpha| %.2f |pi-alpha| %.2f loss cone %.2f deg\n', ...
        abs(alpha),abs(180-alpha),alphaloss);
if abs(alpha)<alphaloss || abs(180-alpha)<alphaloss,
  fprintf(1,'warning: particle in loss cone!\n');
end
	
end

function trace(Xo,mdfile,Ep,Ri,ai,timespec,Tc,Tb,Td,Lm,savefile,pauseOn,npertc)

% Planet's surface
ts = linspace(0,2*pi,100);
Xp = cos(ts);
Yp = sin(ts);

% *dipole* mirror point 
rm = {Xo(1)*cosd(Lm)^3,Xo(1)*cosd(Lm)^2*sind(Lm),0};
if sqrt(rm{1}^2+rm{2}^2)/Re < 1,
  error(sprintf('%s R=%.2f < 1\n%s',...
        'Estimated mirror point within planet',...
	      sqrt(rm{1}^2+rm{2}^2)/Re,...
        'MDisc interpolation problem!'));
end
% B from magnetodisc
%B = mdiscMagneticField3D({rm{1}/Re,rm{2}/Re,rm{3}/Re});
% B from dipole
B = dipoleMagneticField3D(Md,Rm,rm);
% Gyro period at dipole mirror point
tm = 2*pi/(qOverM/gamma*sqrt(B{4}));
% number of iteration such that dt \sim max([tm/3,Tc/5])
tb = linspace(tspan(1),tspan(end),ceil(diff(tspan)./[Tc/npertc]))';
dt = diff(tb(1:2));
% tc/tm ratio of gyro period at equator and dipole mirror point 
fprintf(1,'tm=%.2g, tm/dt=%.2f tc/tm=%.2f, tc/dt=%.2f\n',...
        [tm,tm/dt,Tc/tm,Tc/dt]);

if is_octave(),
  fflush(1);
end

tic
Xb = zeros(length(tb),length(Xo));
Xb(1,:) = BorisInit(Xo,dt);
for i=1:length(tb)-1
  Xb(i+1,:) = BorisIter(Xb(i,:),dt);
end
toc

Zb    = Xb(:,3);
Rcylb = sqrt(Xb(:,1).^2+Xb(:,2).^2);
Rtotb = sqrt(Xb(:,1).^2+Xb(:,2).^2+Xb(:,3).^2);

figure
plot(Rcylb/Re,Zb/Re,Xp,Yp),
axis equal,
xlabel('\rho'), ylabel('z')
legend({'FD'})
if pauseOn,
  fprintf(1,'Ok, press return\n'),
  pause
  fprintf(1,'\n');
else
  drawnow
end

% Kinetic energy in eV for Full Dynamic Boris
v2b = sum(Xb(:,[4:6])'.^2);
%Eb = 1/2*mp*v2b/eV;  % classic
Eb = (1./sqrt(1-v2b/c^2)-1)*mp*c^2/eV; % relativistic
% mean and variance
meanEb = mean(Eb);
stdEb = std(Eb);
fprintf(1,'<Eb> = %10g MeV, std(E) = %10g eV\n', meanEb/1e6, stdEb);

plot(tb,Eb-meanEb);
xlabel('time'), ylabel('E-<E>')
legend({'FD'})
if pauseOn,
  fprintf(1,'Ok, press return\n'),
  pause
  fprintf(1,'\n');
else
  drawnow
end

% Compute b dot gradB and mu for Full Dynamic Boris
[B,gradB] = mdiscMagneticField3D({Xb(:,1)/Re,Xb(:,2)/Re,Xb(:,3)/Re});
b = sqrt(B{4});
BgBb = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./b;
%Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./b;
V2b = Xb(:,4).^2+Xb(:,5).^2+Xb(:,6).^2;
Vparb = (B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6))./b;
Vperb = sqrt(V2b-Vparb.^2);
Vper2b = V2-Vparb.^2;
% Instantaneous first invariant Eq.14 
muib = gamma^2*mp*Vper2b./(2*b);

subplot(211), 
plot(tb,muib,tb,facdv*gamma^2*mp*ones(size(tb))),
ylabel('\mu 1st invariant')
legend({'FD','Initial'})

% Compute latitude and longitude
latb = atan2d(Zb,Rcylb);        % atan(Rcyl/Z)
lonb = atan2d(Xb(:,2),Xb(:,1)); % atan(X/Y);

subplot(212), 
%plot(tb,latb), xlabel('time'); ylabel('Latitude')
aip = abs(ai);
x=[0;sind(aip);0;sind(aip)]; y=[0;cosd(aip);0;-cosd(aip)];
plot(Vperb./sqrt(V2b),Vparb./sqrt(V2b),x,y),
xlabel('v_\perp'); ylabel('v_{||}'); 
legend({'FD'})
if pauseOn,
  fprintf(1,'Ok, press return\n'), 
  pause
  fprintf(1,'\n');
else
  drawnow
end

% L-Shell with L
t = linspace(-pi/2,pi/2,100);
% estimate L from instantaneous L
%coslat = Rcyl./R; Ls = R./coslat.^2/Re;
Lsb = Rtotb.^3./Rcylb.^2/Re;
Le = mean(Lsb);
fprintf(1,'L estimated = %.2f\n', Le);
if 0,
	f=figure; 
	plot(tb,Lsb), 
	pause, close(f);
end
r = Le *cos(t).^2;
ir = find(r>=1);
xe = r(ir).*cos(t(ir));
ye = r(ir).*sin(t(ir));

% plot trajectory Rcyl vs Z
subplot(211), 
plot(Rcylb/Re,Zb/Re,xe,ye,'-.',Xp,Yp), 
xlabel('\rho'),
ylabel('z'),
title(sprintf('X_0=%.2f L=%.2f', Xo(1)/Re, Le))
legend({'FD','L-shell'});
%set(gca,'xlim',[0 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal,

subplot(212),
if 1,
% plot trajectory X vs Y
plot(Xb(:,1)/Re,Xb(:,2)/Re,Xp,Yp), 
xlabel('x'),
ylabel('y'),
else
% 3d plot trajectory
plot3(Xb(:,1)/Re,Xb(:,2)/Re,Xb(:,3)/Re), 
xlabel('x'),
ylabel('y'),
zlabel('z'),
end
title(sprintf('<E>=%2g MeV std(E)=%2g eV',meanEb/1e6,stdEb));
%set(gca,'xlim',[-4.5 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal
if pauseOn,
  fprintf(1,'Ok, press return\n'),
  pause
  fprintf(1,'\n');
else
  drawnow
end


% identify zero crossings
izc = find(latb(1:end-1).*latb(2:end)<0);
Tbe = 2*mean(diff(tb(izc)));
% Taking average of Tb and Tbe with weight 1 and 2
Tbi= (Tb+2*Tbe)/3;
Tbi= Tbe;
fprintf(1,'**** Tbd=%.2f Tbe=%.2f Tbi=%.2f\n',Tb,Tbe,Tbi);
[kvgb,tbb,dtbb,lmb,fitlat] = getbounceperiod(tb,latb,2*pi/Tbi);

subplot(211),
plot(tb,latb,tb,fitlat), xlabel('time'); ylabel('Latitude')
legend({'FD','FIT'})

[kvgd,tdb,dtdb,fitlon] = getdriftperiod(tb,lonb,2*pi/tbb,360);

subplot(212),
plot(tb,lonb,tb,fitlon), xlabel('time'); ylabel('Longitude')
legend({'FD','FIT'})

drawnow

if ~isempty(savefile), % save all trajectories
  Be = mdiscMagneticField3D({Re/Re,0,0});
  Be = sqrt(Be{4});
  save(savefile,'mdfile','Re','Be','Ep','Ri','ai','timespec',...
       'Tc','Tb','Td','Lm','tb','Xb',...
       'Zb','Rcylb','Rtotb','Eb','muib','latb','lonb',...
       'kvgb','tbb','dtbb','lmb','fitlat',...
       'kvgd','tdb','dtdb','fitlon');
end

end

function rn = BorisIter(r,h)

rn = zeros(size(r));

[x,y,z,vx,vy,vz] = deal(r(1),r(2),r(3),r(4),r(5),r(6));

B = mdiscMagneticField3D({x/Re,y/Re,z/Re});

[Bx,By,Bz] = deal(B{1:3});

fac = qOverM/gamma;

%T = -[Bx,By,Bz]'/sqrt(B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
T = fac*[Bx,By,Bz]'*h/2;
S = 2.0 * T /(1+T'*T);

V = [vx,vy,vz]';
R = [x,y,z]';

v = V + cross(V,T);
V = V + cross(v,S);

R = R + V*h;

rn = [R;V];

end

function r0 = BorisInit(r,h)

r0 = zeros(size(r));

[x,y,z,vx,vy,vz] = deal(r(1),r(2),r(3),r(4),r(5),r(6));

B = mdiscMagneticField3D({x/Re,y/Re,z/Re});

[Bx,By,Bz] = deal(B{1:3});

fac = qOverM/gamma;

% dt is now -dt/2 (to go backward in time for dt/2)
h = -h/2;
T = fac*[Bx,By,Bz]'*h/2;
S = 2.0 * T /(1+T'*T);

V = [vx,vy,vz]';
R = [x,y,z]';

v = V + cross(V,T);
V = V + cross(v,S);

r0 = [R;V];

end

function [tc,tb,td] = periods(R,alpha)

B = mdiscMagneticField3D({R/Re,0,0});
Be = mdiscMagneticField3D({Re/Re,0,0});

B  = sqrt(B{4});
Be = sqrt(Be{4});

% Angular gyro frequency
Omegac = abs(qOverM)/gamma*B;

% Gyro period
tc = 2*pi/Omegac;

% Bouncing period
% factor Jupiter 1.2478 = Re/c*sqrt(2)*3.7
% factor Saturn  1.0521 = Re/c*sqrt(2)*3.7
fac = Re/c*sqrt(2)*3.7;
tb = fac*(R/Re)*c/vp*(1-0.4635*sind(alpha)^0.75);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
%tb = R*sqrt(2)/vp*(3.7-1.6*sind(alpha));

% Drift period 
td = 2*pi*abs(qOverM)*Be*Re^3/vp^2/R*(1-1/3*(sind(alpha))^0.62);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
%td = pi*qp*Be*Re^3/(3*EMeV*R)/(0.35+0.15*sind(alpha));
%td = 1.05/(EMeV/(1e6*eV))/(R/Re)/(1+0.43*sind(alpha))*3600;

fprintf(1,'tc=%6.2f s, tb=%6.2f s, td=%6.2f s\n', tc,tb,td);

end

function lm = mirrorlat(alpha)

% numerically 
opts = optimset('fzero');
myfun = @(x,a) cosd(x).^6-sind(a)^2.*sqrt(1+3*sind(x).^2);
lm = fzero(@(x) myfun(x,alpha),[0,90],opts);

fprintf(1,'pitcheq=%.2f lm=%.2f deg\n',alpha,lm);

end

end
