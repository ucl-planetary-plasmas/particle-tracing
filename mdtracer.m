function mdtracer(mdfile,Ep,Ri,ai,timespec,savefile)
% function mdtracer(mdfile,Ep,Ri,ai,timespec,savefile)
%
% Parameterization of Ozturk's 2012 paper 2012
% Trajectories of charged particles trapped in Earth's magnetic field
%
% mdfile matlab magnetodisc filename
% Ep proton energy in MeV
% Ri initial particle equatorial position (in planetary radius)
% ai initial particle pitch angle at equatorial position (in degrees)
% timespec defined as tmax = sum(timespec.*[1,tc,tb,td])
%          where 1 is in units of seconds, tc in units of gyroperiod at Ri, 
%          tb in units of bounce period and td in units of drift period
% savefile is a filename to save the simulation data.
%          if not given, do not save simulation data.
%
% For instance 
%
%   mdtracer('jup_mdisc_kh3e7_rmp90',100,4,30,[0,0,2,0],'my_jup_mdisc_sim');
%
% runs the Boris algorithm, and guiding centre approximation in a Jupiter
% magnetodisc for a proton of 100 MeV at initial position 4 Rp and
% pitch angle 30 deg for 2*tb, i.e. two bounce periods, and save the
% simulation data in file 'my_jup_mdisc_sim.mat'.
%
% Similarly for Saturn
%
%   mdtracer('sat_mdisc_kh2e6_rmp25',100,4,30,[0,0,2,0]);
%
% but don't save the simulation data.

%
% $Id: mdtracer.m,v 1.4 2018/06/18 12:09:36 patrick Exp $
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


% solver to use
solver = 'ode15s';
% ODE solver options
%opts = odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','off','vectorized','on');
opts = odeset('AbsTol',1e-12,'RelTol',1e-11,'Stats','off','vectorized','on');
%opts = odeset('AbsTol',1e-8,'RelTol',1e-8,'Stats','on','vectorized','on');
%opts = odeset('AbsTol',1e-8,'RelTol',1e-6,'Stats','on','vectorized','on');

% close all figures
close all


% planet magnetodisc global filename
global MDFile
MDFile = mdfile;
% clear function with persistent variables
clear MDiscField
load(MDFile,'MD')

% Equatorial radius in m
Re = MD.planet.r0;

% Planet dipole field
% Magnetic equator field strength in T
%B0 = MD.planet.B0;
% corrected to match magnetodisc file
% Bthdip(mu=0,r=1)=8.176230e-01 -> 3.499426e-04, Md=1.2787e20
% Bth   (mu=0,r=1)=8.172328e-01 -> 3.497756e-04, Md=1.2781e20
%B0 = MD.planet.B0 * MD.v2d.Bth(MD.dims.imu0,1);
B0 = MD.planet.B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
Md = [0,0,B0*Re^3];  % Magnetic moment
Rm = [0,0,0];        % Centered moment

% SI constants
mp = 1.6726231e-27;  % particle mass (proton mass)
qp = 1.60217733e-19; % proton charge
eV = 1.60217733e-19; % conversion factor eV to J
c  = 299792458.0;    % speed of light

EMeV = Ep*1e6*eV;    % proton energy in MeV 
vp = sqrt(2/mp*EMeV);  % classic velocity 
fprintf('classic v=%.4g m/s E=%.4g MeV\n',vp,.5*mp*vp^2/(1e6*eV))
vp = c*sqrt(EMeV*(EMeV+2*mp*c^2))/(EMeV+mp*c^2); % relativistic velocity
fprintf('relativ v=%.4g m/s E=%.4g MeV\n',vp,(1/sqrt(1-vp^2/c^2)-1)*mp*c^2/(1e6*eV));
% velocity squared for GC solution
V2 = vp^2;

% relativistic factor
gamma = 1/sqrt(1-V2/c^2);

% proton charge-to-mass ratio
qOverM = qp/mp;

close all
% initial conditions
[X,Xgc,mu,Tc,Tb,Td,Lm] = init(Ri*Re,vp,ai);

% factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
facdv = mu/gamma^2/mp; 
tspan = [0,sum(timespec.*[1,Tc,Tb,Td])];
fprintf(1,'tmax=%.2f (%.1f tc, %.2f tb, %.2f td)\n',...
        tspan(end)./[1,Tc,Tb,Td]);

trace(X,Xgc);

function [Xo,Xogc,mu,tc,tb,td,lm]=init(R,v,alpha)

% initial condition x0,y0,z0,vx0,vy0,vz0
Xo = [R;            0;            0;...
      0;v*sind(alpha);v*cosd(alpha)];

% initial guiding centre x0,y0,z0,vpar0
B = mdiscMagneticField3D(Rm,{Xo(1)/Re,Xo(2)/Re,Xo(3)/Re});
b = sqrt(B{4});

fprintf(1,'R=%.2f pitch angle=%.0f B=%.5g nT\n', R/Re, alpha, 1e9*b)

% smaller vpar gives larger vper thus guiding centre motion
%vpar = 0.9220*sum(Xo(4:6).*[B{1};B{2};B{3}]/b);
vpar = sum(Xo(4:6).*[B{1};B{2};B{3}]/b);
vper = sqrt(V2-vpar^2);
% gyro radius
rho = abs(gamma/qOverM*vper/b);
% add or substract gyroradius depending on sign of q and magnetic moment
Xogc = [R-sign(qOverM*Md(3))*rho;0;0;vpar];
B = mdiscMagneticField3D(Rm,{Xogc(1)/Re,Xogc(2)/Re,Xogc(3)/Re});
b = sqrt(B{4});
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

function trace(Xo,Xogc)

solverh = fcnchk(solver);

% Planet's surface
ts = linspace(0,2*pi,100);
Xp = cos(ts);
Yp = sin(ts);

tic
rm = {Xo(1)*cosd(Lm)^3,Xo(1)*cosd(Lm)^2*sind(Lm),0};
B = mdiscMagneticField3D(Rm,{rm{1}/Re,rm{2}/Re,rm{3}/Re});
% Gyro period at mirror point
tm = 2*pi/(qOverM/gamma*sqrt(B{4}));
%tb = linspace(0,t(end),floor(t(end)/max(diff(t))));
tb = linspace(tspan(1),tspan(end),6*fix(diff(tspan)/tm))';
%tb = t;
hb = diff(tb);
fprintf(1,'tm=%.2f, tc/tm=%.2f, tc/hb(1)=%.2f\n',[tm,fix(Tc/tm),fix(Tc/hb(1))]);
Xb = zeros(length(tb),length(Xo));
Xb(1,:) = BorisInit(Xo,hb(1));
for i=1:length(tb)-1
  Xb(i+1,:) = BorisIter(Xb(i,:),hb(i));
end
toc

Zb    = Xb(:,3);
Rcylb = sqrt(Xb(:,1).^2+Xb(:,2).^2);
Rtotb = sqrt(Xb(:,1).^2+Xb(:,2).^2+Xb(:,3).^2);

tic
[tgc,Xgc] = solverh(@gcdynamic3d,tspan,Xogc,opts);
toc

Zgc    = Xgc(:,3);
Rcylgc = sqrt(Xgc(:,1).^2+Xgc(:,2).^2);
Rtotgc = sqrt(Xgc(:,1).^2+Xgc(:,2).^2+Xgc(:,3).^2);

figure
plot(Rcylb/Re,Zb/Re,Rcylgc/Re,Zgc/Re,Xp,Yp),
axis equal,
xlabel('\rho'), ylabel('z')
legend({'FD','GC'})
fprintf(1,'Ok, press return\n'),
pause
fprintf(1,'\n');

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
fprintf(1,'Ok, press return\n'),
pause
fprintf(1,'\n');

% Compute b dot gradB and mu for Full Dynamic Boris
[B,gradB] = mdiscMagneticField3D(Rm,{Xb(:,1)/Re,Xb(:,2)/Re,Xb(:,3)/Re});
BgBb = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./sqrt(B{4});
Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./B{4};
% Instantaneous first invariant Eq.14 
muib = gamma^2*mp*Vper2b./(2*sqrt(B{4}));

% Compute b dot gradB and mu for Guiding Centre
[B,gradB] = mdiscMagneticField3D(Rm,{Xgc(:,1)/Re,Xgc(:,2)/Re,Xgc(:,3)/Re});
BgBgc = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./sqrt(B{4});
Vper2gc = V2-Xgc(:,4).^2;
% Instantaneous first invariant Eq.14 
muigc = gamma^2*mp*Vper2gc./(2*sqrt(B{4}));

subplot(211), 
plot(tb,muib,tgc,muigc,tgc,facdv*gamma^2*mp*ones(size(tgc))),
ylabel('\mu 1st invariant')
legend({'FD','GC','Initial'})

% Compute latitude and longitude
latb = atan2d(Zb,Rcylb);        % atan(Rcyl/Z)
lonb = atan2d(Xb(:,2),Xb(:,1)); % atan(X/Y);
latgc = atan2d(Zgc,Rcylgc);
longc = atan2d(Xgc(:,2),Xgc(:,1));

subplot(212), 
%plot(tgc,Xgc(:,3)/Re), ylabel('z')
plot(tb,latb,tgc,latgc), xlabel('time'); ylabel('Latitude')
legend({'FD','GC'})

fprintf(1,'Ok, press return\n'), 
pause
fprintf(1,'\n');
%subplot(211), plot(tgc,BgBgc),
%subplot(212), plot(Xgc(:,3)/Re,BgBgc),pause

% L-Shell with L
t = linspace(-pi/2,pi/2,100);
% estimate L from instantaneous L
%coslat = Rcyl./R; Ls = R./coslat.^2/Re;
Lsb = Rtotb.^3./Rcylb.^2/Re;
Lsgc = Rtotgc.^3./Rcylgc.^2/Re;
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
plot(Rcylb/Re,Zb/Re,Rcylgc/Re,Zgc/Re,xe,ye,'-.',Xp,Yp), 
xlabel('\rho'),
ylabel('z'),
title(sprintf('%s x0=%.0f x0gc=%.4f L=%.4f', ...
      solver, Xo(1)/Re, Xogc(1)/Re,Le))
legend({'FD','GC','L-shell'});
%set(gca,'xlim',[0 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal,

subplot(212),
if 1,
% plot trajectory X vs Y
plot(Xb(:,1)/Re,Xb(:,2)/Re,...
     Xgc(:,1)/Re,Xgc(:,2)/Re,Xp,Yp), 
xlabel('x'),
ylabel('y'),
else
% 3d plot trajectory
plot3(Xb(:,1)/Re,Xb(:,2)/Re,Xb(:,3)/Re,...
      Xgc(:,1)/Re,Xgc(:,2)/Re,Xgc(:,3)/Re), 
xlabel('x'),
ylabel('y'),
zlabel('z'),
end
title(sprintf('<E>=$%2g MeV std(E)=%2g eV',meanEb/1e6,stdEb));
%set(gca,'xlim',[-4.5 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal

fprintf(1,'Ok, press return\n'),
pause
fprintf(1,'\n');

subplot(211),
plot(tb,latb,tgc,latgc), xlabel('time'); ylabel('Latitude')
legend({'FD','GC'})

[tbb,dtbb] = getbounceperiod(tb,latb,2*pi/Tb);
[tbgc,dtbgc] = getbounceperiod(tgc,latgc,2*pi/Tb);

subplot(212),
plot(tb,lonb,tgc,longc), xlabel('time'); ylabel('Longitude')
legend({'FD','GC'})


[tdb,dtdb] = getdriftperiod(tb,lonb,2*pi/tbb,360);
[tdgc,dtdgc] = getdriftperiod(tgc,longc,2*pi/tbgc,360);

drawnow

if exist('savefile','var') & ~isempty(savefile), % save all trajectories
  save(savefile,'mdfile','Ep','Ri','ai','timespec',...
       'Tc','Tb','Td','Lm','tb','Xb','tgc','Xgc',...
       'tbb','dtbb','tbgc','dtbgc','tdb','dtdb','tdgc','dtdgc');
end

end

function rn = BorisIter(r,h)

rn = zeros(size(r));

[x,y,z,vx,vy,vz] = deal(r(1),r(2),r(3),r(4),r(5),r(6));

B = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});

[Bx,By,Bz] = deal(B{1:3});

fac = qOverM/gamma;

%T = -tan(B{4}*h/2.0)*[Bx,By,Bz]'/B{4};
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

B = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});

[Bx,By,Bz] = deal(B{1:3});

fac = qOverM/gamma;

T = fac*[Bx,By,Bz]'*h/2;


V = [vx,vy,vz]';
R = [x,y,z]';

V = V - cross(V,T);

% don't forget to keep v constant
V = V * vp/sqrt(sum(V.^2));

r0 = [R;V];

end

function dr_dt = gcdynamic3d(t,r)

[x,y,z,vpar] = deal(r(1,:),r(2,:),r(3,:),r(4,:));

[B,gradB] = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});
[Bx,By,Bz,B2] = deal(B{:});
[gradBx,gradBy,gradBz] = deal(gradB{:});

B = B2.^.5;
iB = 1./B;
bx = Bx.*iB; by = By.*iB; bz = Bz.*iB;

fac = gamma./(2*qOverM*B2).*(V2 + vpar.^2);
dx = fac.*(by.*gradBz - bz.*gradBy) + vpar.*bx;
dy = fac.*(bz.*gradBx - bx.*gradBz) + vpar.*by;
dz = fac.*(bx.*gradBy - by.*gradBx) + vpar.*bz;
%dv = -(V2 - vpar.^2)./(2*B).*(bx.*gradBx+by.*gradBy+bz.*gradBz);
dv = -facdv*(bx.*gradBx+by.*gradBy+bz.*gradBz);

dr_dt = [dx;dy;dz;dv];

end

function [tc,tb,td] = periods(R,alpha)

B = mdiscMagneticField3D(Rm, {R/Re,0,0});
Be = mdiscMagneticField3D(Rm, {Re/Re,0,0});

B  = sqrt(B{4});
Be = sqrt(Be{4});

% Angular gyro frequency
Omegac = qOverM/gamma*B;

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
td = 2*pi*qOverM*Be*Re^3/vp^2/R*(1-1/3*(sind(alpha))^0.62);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
%td = pi*qp*Be*Re^3/(3*EMeV*R)/(0.35+0.15*sind(alpha));
%td = 1.05/(EMeV/(1e6*eV))/(R/Re)/(1+0.43*sind(alpha))*3600;

fprintf(1,'tc=%6.2f s, tb=%6.2f s, td=%6.2f s\n', tc,tb,td);

end

function lm = mirrorlat(alpha)

if 0,
% symbolically takes longer time and requires symbolic toolbox
M = sym('M','positive');
eqn = M.^6 + 3*M*(sind(alpha)).^4 == 4*(sind(alpha)).^4;
COS2M = vpasolve(eqn,M,0.5);
COSM = double(sqrt(COS2M(2)));
lm = acosd(COSM);
fprintf(1,'pitcheq=%.2f lm=%.2f deg\n',alpha,lm);
end

% numerically 
opts = optimset('fzero');
myfun = @(x,a) cosd(x).^6-sind(a)^2.*sqrt(1+3*sind(x).^2);
lm = fzero(@(x) myfun(x,alpha),45,opts);

fprintf(1,'pitcheq=%.2f lm=%.2f deg\n',alpha,lm);

end

end
