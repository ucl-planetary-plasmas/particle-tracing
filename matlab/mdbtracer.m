function mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc)
% function % mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc)
%
% Magnetodisc Boris particle tracer
%
% mdfile   : magnetodisc mat-file filename.
% partype  : particle type 'p' -> proton, 'e' -> electron, 
%            'O+', 'O++', 'S+', 'S++', 'S+++'.
% Ep       : particle energy [MeV].
% Ri       : initial particle equatorial position in planet's radius [Rp].
% ai       : initial particle pitch angle (0...180 deg).
% timespec : max simulation time defined as sum(timespec.*[1,tc,tb,td]) where
%              * 1 is in units of seconds, 
%              * tc in units of gyroperiod at initial Ri,
%              * tb in units of dipole bounce period,
%              * td in units of dipole drift period.
% savefile : filename to save the simulation data,
%            if not given or empty, do not save simulation data.
% pauseOn  : flag to pause on/ pause off.
%            if not given or empty, pause on.
% plotOn   : flag to plot on/ plot off.
%            if not given or empty, plot on.
% npertc   : optional number of Boris iterations per gyroperiod at Ri,
%            default value 5.
%
% For instance 
%
%   mdfile = 'jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
%   partype = 'p';                    % proton
%   Ep = 100;                         % energy [MeV]
%   Ri = 4;                           % initial equatorial distance [Rp]
%   ai = 30;                          % initial equatorial pitch angle [deg]
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
%   Ep = 100;                         % energy [MeV]
%   Ri = 4;                           % initial equatorial distance [Rp]
%   ai = 30;                          % initial equatorial pitch angle [deg]
%   timespec = [0,0,2,0];             % run for 2 dipole bounce periods
%
%   mdbtracer(mdfile,partype,Ep,Ri,ai,timespec);
%
% but don't save the simulation data.

%
% $Id: mdbtracer.m,v 1.16 2026/07/08 17:33:36 patrick Exp $
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

if ~exist('savefile','var'), % define savefile if not given
  savefile = [];
end

if ~exist('pauseOn','var') | isempty(pauseOn), % define pauseOn default
  pauseOn = true;
end

if ~exist('plotOn','var') | isempty(plotOn), % define plotOn default
  plotOn = true;
end

if ~exist('npertc','var') | isempty(npertc), % define npertc=5 default
  npertc = 5;
end

% close all figures
if plotOn,
close all
end

% planet magnetodisc global filename
global MDFile
MDFile = mdfile;
% clear function with persistent variables
clear MDiscField
load(MDFile,'MD')

Rmp = max(MD.dims.r);
if Ri>=Rmp,
  error(['Ri=' num2str(Ri) ' is larger than Rmp=' num2str(Rmp)]);
end

% Planet's equatorial radius in m
Re = MD.planet.r0;
% Plasma angular velocity
%Omr = griddedInterpolant(MD.disc.L,MD.planet.Omega*MD.disc.omega); 
%Omega = Omr(Ri);
Omega = MD.planet.Omega*interp1(MD.disc.L,MD.disc.omega,Ri,'linear','extrap');
% Planet's angular velocity
%Omega = MD.planet.Omega;
% No planet's rotation
%Omega = 0;
fprintf(1,'Omega planet = %.4g rad/s, Omega(Ri) = %.4g rad/s\n',...
           MD.planet.Omega, Omega);
% Rotation axis tilt with respect to magnetic dipole 
% for Saturn 0 deg
%tilt = 0*pi/180;
% for Earth and Jupiter 10 deg
tilt = 10*pi/180;
% Rotation vector in the (x,z)-plane
Omp = [Omega*sin(tilt),0,Omega*cos(tilt)]';
fprintf(1,'Dipole tilt = %.1f deg\n',tilt*180/pi);

% magnetic equator field strength in T
B0 = MD.planet.B0;
% corrected to match Bthdip(mu=0,r=1)
B0 = B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
% Dipole magnetic moment (Beq = \mu_0/(4 pi Req^3) Md)  in A m^2
Md = [0,0,B0*Re^3];  
% Centered dipole magnetic moment
Rm = [0,0,0];        

% SI constants
e   = 1.60217733e-19; % elementary charge C
eV  = 1.60217733e-19; % conversion factor from eV to J
c   = 299792458.0;    % speed of light m s^-1
m_e = 9.1093897e-31;  % electron mass kg
m_p = 1.6726231e-27;  % proton mass kg
amu = 1.6605402e-27;  % atomic mass unit kg

switch partype,
  case 'p', % proton
    parname = 'proton';
    mp = m_p;
    qp = e;
  case 'e', % electron
    parname = 'electron';
    mp = m_e;
    qp = -e;
  case 'O+', % oxygen ion
    parname = 'O+';
    mp = 15.999 * amu;
    qp = e;
  case 'O++', % oxygen ion
    parname = 'O++';
    mp = 15.999 * amu;
    qp = 2 * e;
  case 'S+', % sulfur ion
    parname = 'S+';
    mp = 32.07 * amu;
    qp = e;
  case 'S++', % sulfur ion
    parname = 'S++';
    mp = 32.07 * amu;
    qp = 2 * e;
  case 'S+++', % sulfur ion
    parname = 'S+++';
    mp = 32.07 * amu;
    qp = 3 * e;
end

EMeV = Ep*1e6*eV;      % particle energy in MeV
Vpc = sqrt(2/mp*EMeV);  % classic velocity 
fprintf('particle: %s (m=%.6g amu, Z=%d)\n',parname,mp/amu,qp/e) 
fprintf('classic v=%.4g m/s E=%.4g MeV\n',Vpc,.5*mp*Vpc^2/(1e6*eV))
Vpr = c*sqrt(EMeV*(EMeV+2*mp*c^2))/(EMeV+mp*c^2); % relativistic velocity
fprintf('relativ v=%.4g m/s E=%.4g MeV\n',Vpr,(1/sqrt(1-Vpr^2/c^2)-1)*mp*c^2/(1e6*eV));
% velocity squared
V2 = Vpr^2;

% relativistic factor
gamma = 1/sqrt(1-V2/c^2);

% particle charge-to-mass ratio
qOverM = qp/mp;

% initial conditions
[X,mu,Tc,Tb,Td,Lm,rg] = init(Ri*Re,Vpr,ai);

% factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
facdv = mu/gamma^2/mp; 
tspan = [0,sum(timespec.*[1,Tc,Tb,Td])];
fprintf(1,'tmax=%.10g (%.10g tc, %.10g tb, %.10g td)\n',...
        tspan(end)./[1,Tc,Tb,Td]);

trace(X,mdfile,Ep,Ri,ai,timespec,Tc,Tb,Td,Lm,rg,savefile,pauseOn,plotOn,npertc);

function [Xo,mu,tc,tb,td,lm,rg]=init(R,v,alpha)

% initial condition x0,y0,z0,vx0,vy0,vz0
Xo = [R,            0,            0,...
      0,v*sind(alpha),v*cosd(alpha)];

% magnetic field at initial position
B = mdiscMagneticField3D(Rm,{Xo(1)/Re,Xo(2)/Re,Xo(3)/Re});
Bm = B{4};

fprintf(1,'Ri=%.2f pitch angle=%.0f B=%.5g nT\n', R/Re, alpha, 1e9*Bm)

% v parallel and perpendicular amplitudes
vpar = dot(Xo(1,4:6),[B{1};B{2};B{3}])/Bm;
vper = sqrt(V2-vpar^2);

% gyro radius at equator
rg = abs(gamma/qOverM*vper/Bm);
fprintf(1,'rg=%.4gRe\n', rg/Re);

% first invariant mu
mu = gamma^2*mp*vper^2/(2*Bm);

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

function trace(Xo,mdfile,Ep,Ri,ai,timespec,Tc,Tb,Td,Lm,rg,savefile,pauseOn,plotOn,npertc)

% Planet's surface
tp = linspace(0,2*pi,100);
Xp = cos(tp);
Yp = sin(tp);

% *dipole* mirror point 
rm = {Xo(1)*cosd(Lm)^3,Xo(1)*cosd(Lm)^2*sind(Lm),0};
fprintf(1,'Mirror point = [%12g, %12g, %12g] Re\n', cell2mat(rm)/Re);
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
tm = 2*pi/(qOverM/gamma*B{4});
% number of iteration such that dt \sim max([tm/3,Tc/5])
nsteps = ceil(diff(tspan)/(Tc/npertc));
fprintf(1,'%.10g %.10g %.10g %.10g\n',tspan(1),tspan(end),Tc,nsteps)
ts = linspace(tspan(1),tspan(end),nsteps)';
dt = diff(ts(1:2));
% tc/tm ratio of gyro period at equator and dipole mirror point 
fprintf(1,'tm=%.6g, tm/dt=%.6f tc/tm=%.6f, tc/dt=%.6f\n',...
        [tm,tm/dt,Tc/tm,Tc/dt]);

if exist('is_octave')==2 && is_octave(),
  fflush(1);
end

s1.on = true; 
s1.name = 'Boris';
s1 = initVar(s1,ts,dt);

s2.on = true;
s2.name = 'BorisCF';
s2 = initVar(s2,ts,dt);

s1.X(1,:) = BorisInit(Xo,dt);
%s1.rg(1,:) = [rg,0,0];
%s1.Xg(1,:) =  s1.X(1,1:3)-s1.rg(1,:);
%s1.rg(1,:)/Re,s1.Xg(1,:)/Re

Rp = s1.X(1,1:3); Vp = .5*(s1.X(1,4:6)+s1.X(1,4:6));
[s1.Xg(1,:),s1.rg(1,:)] = getGyroRadius(Rp,Vp);
s1.rg(1,2:3) = 0; % force initially along x-axis
%s1.rg(1,:)/Re,s1.Xg(1,:)/Re

s1 = runSimul(s1);

% Same init for s2
s2.X(1,:) = s1.X(1,:);
s2.rg(1,:) = s1.rg(1,:);
s2.Xg(1,:) = s1.Xg(1,:);

s2 = runSimulCF(s2);

%s1.Xg(1:5,:)/Re,s1.rg(1:5,:)/Re
%s2.Xg(1:5,:)/Re,s2.rg(1:5,:)/Re

s1 = computeDiagnostics(s1);

s2 = computeDiagnostics(s2);

if plotOn,
figure
subplot(211)
plot(s1.Rcyl/Re,s1.Z/Re,s2.Rcyl/Re,s2.Z/Re,Xp,Yp),
axis equal,
xlabel('\rho'), ylabel('z')
legend({'Boris','BorisCF'},'Location','Best')
subplot(212)
plot(s1.Rgcyl/Re,s1.Zg/Re,s2.Rgcyl/Re,s2.Zg/Re,Xp,Yp),
axis equal,
xlabel('\rho'), ylabel('z')
legend({'Boris-GC','BorisCF-GC'},'Location','Best')
if pauseOn,
  fprintf(1,'Ok, press return\n'),
  pause
  fprintf(1,'\n');
else
  drawnow
end
end

% Kinetic energy in eV for CF Dynamic Boris
%save('v2b_mat','s1.v2','s1.v2c')
%plot(s1.t,sqrt(s1.v2),s2.t,sqrt(s2.v2)); title('|v|^2'), pause

if plotOn,
subplot(311)
plot(s1.t,s1.E,s2.t,s2.E,s2.t,s2.Ek); 
title('Energy'); 
legend({'Boris','BorisCF','BorisCF-Ek'},'Location','Best')
subplot(312)
plot(s2.t,s2.Ep),
legend({'BorisCF-Ep'},'Location','Best')
subplot(313)
plot(s1.t,(s1.E-s1.meanE)/s1.meanE,s2.t,(s2.E-s2.meanE)/s2.meanE);
xlabel('time'), ylabel('(E-<E>)/<E>')
legend({'Boris','BorisCF'},'Location','Best')
if pauseOn,
  fprintf(1,'Ok, press return\n'),
  pause
  fprintf(1,'\n');
else
  drawnow
end
end

if plotOn,
subplot(311), 
plot(s1.t,s1.mui,s2.t,s2.mui,ts,facdv*gamma^2*mp*ones(size(ts))),
xlabel('time')
ylabel('\mu 1st invariant')
legend({'Boris','BorisCF','Initial'},'Location','Best')

subplot(312),
plot(s1.t,s2.ai,s2.t,s2.ai),
xlabel('time')
set(gca,'ylim',[0,180],'ytick',[0,90,180])
ylabel('\alpha pitch angle')
legend({'Boris','BorisCF'},'Location','Best')

subplot(313), 
aip = abs(ai);
xd=[0;sind(aip);0;sind(aip)]; 
yd=[0;cosd(aip);0;-cosd(aip)];
[~,im]=min(s1.Vper);
x1=[0;s1.Vper(im);0;s1.Vper(im)]/sqrt(s1.V2(im)); 
y1=[0;s1.Vpar(im);0;-s1.Vpar(im)]/sqrt(s1.V2(im));
[~,im]=min(s2.Vper);
x2=[0;s2.Vper(im);0;s2.Vper(im)]/sqrt(s2.V2(im)); 
y2=[0;s2.Vpar(im);0;-s2.Vpar(im)]/sqrt(s2.V2(im));
plot(s1.Vper./sqrt(s1.V2),s1.Vpar./sqrt(s1.V2),...
     s2.Vper./sqrt(s2.V2),s2.Vpar./sqrt(s2.V2),...
     x1,y2,x2,y2,xd,yd),
xlabel('v_\perp'); ylabel('v_{||}'); 
legend({'Boris','BorisCF','BMP','FMP','DMP'},'Location','Best')
if pauseOn,
  fprintf(1,'Ok, press return\n'), 
  pause
  fprintf(1,'\n');
else
  drawnow
end
end

% L-Shell with L
tl = linspace(-pi/2,pi/2,100);
if 0,
	f=figure; 
	plot(s1.t,s1.L,s2.t,s2.L), 
	pause, close(f);
end
r1 = s1.meanL *cos(tl).^2;
ir = find(r1>=1);
x1 = r1(ir).*cos(tl(ir));
y1 = r1(ir).*sin(tl(ir));

r2 = s2.meanL *cos(tl).^2;
ir = find(r2>=1);
x2 = r2(ir).*cos(tl(ir));
y2 = r2(ir).*sin(tl(ir));

if plotOn,
% plot trajectory Rcyl vs Z
subplot(211), 
plot(s1.Rcyl/Re,s1.Z/Re,s2.Rcyl/Re,s2.Z/Re,x1,y1,'-.',x2,y2,'-.',Xp,Yp), 
xlabel('\rho'),
ylabel('z'),
title(sprintf('X_0=%.2f L=%.2f', Xo(1)/Re, s1.meanL))
legend({'Boris','BorisCF','L-shell'},'Location','Best');
%set(gca,'xlim',[0 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal,

subplot(212),
if 1,
% plot trajectory X vs Y
plot(s1.X(:,1)/Re,s1.X(:,2)/Re,s2.X(:,1)/Re,s2.X(:,2)/Re,...
     s1.Xg(:,1)/Re,s1.Xg(:,2)/Re,s2.Xg(:,1)/Re,s2.Xg(:,2)/Re,Xp,Yp), 
xlabel('x'),
ylabel('y'),
legend({'Boris','BorisCF','Boris-GC','BorisCF-GC'},'Location','Best');
else
% 3d plot trajectory
plot3(s1.X(:,1)/Re,s1.X(:,2)/Re,s1.X(:,3)/Re), 
xlabel('x'),
ylabel('y'),
zlabel('z'),
end
title(sprintf('<E>=%2g MeV std(E)=%2g eV',s1.meanE/1e6,s1.stdE));
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
end

if ~isempty(s1.latfit) && ~isempty(s2.latfit) && plotOn,
subplot(211),
plot(s1.t,s1.lat,s1.t,s1.latfit.f,s2.t,s2.lat,s2.t,s2.latfit.f),
xlabel('time'); ylabel('Latitude')
legend({'Boris','Boris-FIT','BorisCF','BorisCF-FIT'})

subplot(212),
plot(s1.t,s1.lon,s1.t,s1.lonfit.f,s2.t,s2.lon,s2.t,s2.lonfit.f), 
xlabel('time'); ylabel('Longitude')
legend({'Boris','Boris-FIT','BorisCFF','BorisCF-FIT'})

drawnow
end

if ~isempty(savefile), % save all trajectories
  Be = mdiscMagneticField3D(Rm,{Re/Re,0,0});
  Be = sqrt(Be{4});
	fprintf(1,'Results saved to: %s\n',savefile);
  save(savefile,'mdfile','Re','Be','partype','parname','Ep','Ri','ai',...
       'gamma','qOverM','Vpc','Vpr',...
       'Omega','Omp','timespec','npertc',...
       'Tc','Tb','Td','Lm','rg','ts','s1','s2');
end

end

function s = runSimulCF(s)

backspaces = '';
niter = length(s.t)-1;
frq = fix(niter/20); % 20 -> every 5% (100/20)

tic
for i=1:niter,
  s.X(i+1,:) = BorisIterCF(s.X(i,:), s.dt);
	Rp = s.X(i,1:3); Vp = .5*(s.X(i,4:6)+s.X(i+1,4:6));
  [s.Xg(i+1,:),s.rg(i+1,:)] = getGyroRadius(Rp,Vp);
	if mod(i,frq)==0,
	  progStr = sprintf('Progress %s: %d/%d [%3d%%%%]',...
		                  s.name,i,niter,round(i/niter*100));
		fprintf([backspaces,progStr]);
		backspaces = repmat('\b',1,length(progStr));
	end
end
elapsed = toc;
progStr = sprintf('Progress %s: %d/%d [%3d%%%%] (%.4g s, %.2f it/s)\n',...
                  s.name,i,niter,round(i/niter*100), elapsed,niter/elapsed);
fprintf([backspaces,progStr]);

end 

function s = runSimul(s)

backspaces = '';
niter = length(s.t)-1;
frq = fix(niter/20); % 20 -> every 5% (100/20)

tic
for i=1:niter,
  s.X(i+1,:) = BorisIter(s.X(i,:), s.dt);
	Rp = s.X(i,1:3); Vp = .5*(s.X(i,4:6)+s.X(i+1,4:6));
  [s.Xg(i+1,:),s.rg(i+1,:)] = getGyroRadius(Rp,Vp);
	if mod(i,frq)==0,
	  progStr = sprintf('Progress %s: %d/%d [%3d%%%%]',...
		                  s.name,i,niter,round(i/niter*100));
		fprintf([backspaces,progStr]);
		backspaces = repmat('\b',1,length(progStr));
	end
end
elapsed = toc;
progStr = sprintf('Progress %s: %d/%d [%3d%%%%] (%.4g s, %.2f it/s)\n',...
                  s.name,i,niter,round(i/niter*100), elapsed,niter/elapsed);
fprintf([backspaces,progStr]);

end

function s = computeDiagnostics(s)

s.Z    = s.X(:,3);
s.Rcyl = sqrt(sum(s.X(:,1:2).^2,2));
s.Rtot = sqrt(sum(s.X(:,1:3).^2,2));

s.Zg = s.Xg(:,3);
s.Rgcyl = sqrt(sum(s.Xg(:,1:2).^2,2));
s.Rgtot = sqrt(sum(s.Xg(:,1:3).^2,2));

% gyro radius
s.rhog = sqrt(sum(s.rg(:,1:3).^2,2));

%s.Vx = 0.5*[s.X(1:end-1,4)+s.X(2:end,4);s.X(end,4)+Xend(4)];
%s.Vy = 0.5*[s.X(1:end-1,5)+s.X(2:end,5);s.X(end,5)+Xend(5)];
%s.Vz = 0.5*[s.X(1:end-1,6)+s.X(2:end,6);s.X(end,6)+Xend(6)];

s.v2 = sum(s.X(:,4:6).^2,2);
%s.v2c = s.v2-mean(s.v2);
%fprintf(1,'%s: v2 %10g %10g\n',s.name, mean(s.v2), var(s.v2))
%fprintf(1,'%s: v2c %10g %10g\n',s.name, mean(s.v2c), var(s.v2c))

%s.v21 = s.Vx.^2+s.Vy.^2+s.Vz.^2;
%size(s.v2),size(s.v21)
%plot(s.t,sqrt(s.v2b),s.t,sqrt(s.v21)); title('|v|^2'), pause

%s.Ek = 1/2*mp*s.v2b/eV;  % classic
s.Ek = (1./sqrt(1-s.v2/c^2)-1)*mp*c^2/eV; % relativistic

if contains(s.name,'CF'),
  % Only valid for Omega along z-axis
  %s.Ep = -.5*mp*Omega^2*Rcylf.^2/eV; % relativistic
  s.Ep = -.5*mp*((Omp(:)'*Omp(:))*s.Rtot.^2-...
              (Omp(:)'*s.X(:,1:3)')'.^2)/eV; % relativistic
  % Adjust potential energy of rotation so that minimum is zero
  s.Ep = s.Ep-min(s.Ep(:));
  s.E = s.Ek+s.Ep;
else
  s.E = s.Ek;
end
% mean and variance
s.meanE = mean(s.E);
s.stdE = std(s.E);
fprintf(1,'%s <E> = %10g MeV, std(E) = %10g eV std(E)/<Eb> = %10g\n',...
        s.name, s.meanE/1e6, s.stdE, s.stdE/s.meanE);

% grad B and curvature parameters 
[s.B,s.gradB,s.curvB] = mdiscMagneticField3D(Rm,...
                     {s.X(:,1)/Re,s.X(:,2)/Re,s.X(:,3)/Re});

% magnetic-field scale height/scale length
s.Lb = s.gradB{5};
% magnetic-field inhomogeneity parameter or magnetization parameter
% or adiabaticity parameter or gyrokinetic ordering parameter
s.Eb = s.rhog./s.Lb;

% local radius of curvature of the magnetic field line
s.Rc = s.curvB{5};
% curvature parameter or curvature magnetization parameter
% or magnetic-field inhomogeneity parameter
% or guiding-center ordering parameter due to curvature
s.Ec = s.rhog./s.Rc;

Bm = s.B{4};
% apparently slower
%B = [s.B{1},s.B{2},s.B{3}];
%gradB = [s.gradB{1},s.gradB{2},s.gradB{3}];
%s.bgradB = sum(B.*gradB,2)./Bm;
%s.Vpar = sum(B.*s.X(:,4:6),2)./Bm;
s.bgradB = (s.B{1}.*s.gradB{1}+s.B{2}.*s.gradB{2}+s.B{3}.*s.gradB{3})./Bm;

%Vper2 = V2-(B{1}.*X(:,4)+B{2}.*X(:,5)+B{3}.*X(:,6)).^2./b;
s.V2 = sum(s.X(:,4:6).^2,2);
s.Vpar = (s.B{1}.*s.X(:,4)+s.B{2}.*s.X(:,5)+s.B{3}.*s.X(:,6))./Bm;
%Vparb1 = (B{1}.*Vx+B{2}.*Vy+B{3}.*Vz)./Bm;
%plot(tb,(Vparb-Vparb1)./Vparb), title('vpar'); pause
s.Vper = sqrt(s.V2-s.Vpar.^2);
%Vperb1 = sqrt(V2b-Vparb1.^2);
s.Vper2 = s.V2-s.Vpar.^2;
% Instantaneous first invariant Eq.14 
s.mui = gamma^2*mp*s.Vper2./(2*Bm);
% pitch angle
s.ai = atan2d(s.Vper,s.Vpar);

% Compute latitude and longitude
s.lat = atan2d(s.Z,s.Rcyl);        % atan(Rcyl/Z)
s.lon = 180/pi*unwrap(atan2(s.X(:,2),s.X(:,1))); % atan(X/Y);

% estimate L from instantaneous L
%coslat = Rcyl./R; Ls = R./coslat.^2/Re;
s.L = s.Rtot.^3./s.Rcyl.^2/Re;
s.meanL = mean(s.L);
fprintf(1,'%s: L estimated = %.2f\n', s.name, s.meanL);

% bounce period fit
fprintf(1,'**** Bounce period %s\n',s.name);
% identify zero crossings
izc = find(s.lat(1:end-1).*s.lat(2:end)<0);
% average time between zero-crossings times 2 = exact period
s.Tbe = 2*mean(diff(s.t(izc)));
% initial period value for fit
% average of dipole Tb and Tbe with weight 1 and 2
s.Tbi= (Tb+2*s.Tbe)/3;
% exact period
s.Tbi= s.Tbe;
if isfinite(s.Tbe),
fprintf(1,'Init Tbd=%.2f Tbe=%.2f Tbi=%.2f\n',Tb,s.Tbe,s.Tbi);
s.latfit = getbounceperiod(s.t,s.lat,2*pi/s.Tbi);
s.latfit1 = fitTb(s.t,s.lat,s.Tbi);
else
s.latfit = [];
s.latfit1 = [];
end

% drift period fit
if ~isempty(s.latfit),
fprintf(1,'**** Drift period %s\n',s.name);
s.lonfit = getdriftperiod(s.t,s.lon,2*pi/s.latfit.tb,360);
s.lonfit1 = fitTd(s.t,s.lon,s.latfit.tb,360);
else
s.lonfit = [];
s.lonfit1 = [];
end


end


function s = initVar(s,ts,dt)

s.t = ts;
s.dt = dt;
s.X = zeros(length(ts),6);
s.Xg = zeros(length(ts),3);
s.rg = zeros(length(ts),3);
if exist('gpuArray','builtin'),
  s.t = gpuArray(s.t);
  s.X = gpuArray(s.X);
  s.Xg = gpuArray(s.Xg);
	s.rg = gpuArray(s.rg);
end

end

function rn = BorisIterCF(r,h)

x = r(1); y = r(2); z = r(3);

B = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});

Bxyz = [B{1}, B{2}, B{3}];

fac = qOverM/gamma;

%T = -[Bx,By,Bz]'/B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
T = (fac*h/2) * Bxyz;
S = 2.0 * T /(1+T*T');

R = r(1:3);
V = r(4:6);

% https://en.wikipedia.org/wiki/Centrifugal_force
% Centrifugal force dv/dt = -\omega x (\omega x r) 
% Only valid for Omega along z-axis
%dV = Omega^2*[x,y,0]'*h/2;
dV = -cross(Omp,cross(Omp,R))*h/2;
% https://en.wikipedia.org/wiki/Coriolis_force
% Coriolis force dv/dt = - 2 \omega x dr/dt 
% Only valid for Omega along z-axis
%dV = -2*Omega*[-vy,vx,0]'*h/2;
%dV = -2*cross(Omp,V)*h/2;
% Centrifugal and Coriolis forces
% Only valid for Omega along z-axis
%dV = (Omega*[x,y,0]'*h/2 -2*Omega*[-vy,vx,0]')*h/2;
%dV = -cross(Omp,cross(Omp,R))*h/2-2*cross(Omp,V)*h/2;

vminus = V+dV;
vplus = vminus + cross(vminus,T);
vplus = vminus + cross(vplus,S);

V = vplus+dV;
R = R + V*h;

rn = [R,V];

end

function rn = BorisIter(r,h)

x = r(1); y = r(2); z = r(3);

B = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});

Bxyz = [B{1}, B{2}, B{3}];

fac = qOverM/gamma;

%T = -[Bx,By,Bz]'/B{4}*tan(\theta/2.0) = qB/m\Delta{t}/2;
T = (fac*h/2) * Bxyz;
S = 2.0 * T /(1+T*T');

R = r(1:3);
V = r(4:6);

vprime = V + cross(V,T);
V = V + cross(vprime,S);

R = R + V*h;

rn = [R,V];

end

function r0 = BorisInit(r,h)

% https://physics.stackexchange.com/questions/296863/how-to-initialize-bootstrap-the-boris-algorithm

x = r(1); y = r(2); z = r(3);

B = mdiscMagneticField3D(Rm,{x/Re,y/Re,z/Re});

Bxyz = [B{1}, B{2}, B{3}];

fac = qOverM/gamma;

% dt is now -dt/2 (to go backward in time for dt/2)
h = -h/2;
T = fac*Bxyz*h/2;
S = 2.0 * T /(1+T*T');

V = r(4:6);

vprime = V + cross(V,T);
V = V + cross(vprime,S);

r0 = [r(1:3),V];

end

function [tc,tb,td] = periods(R,alpha)

B = mdiscMagneticField3D(Rm,{R/Re,0,0});
Be = mdiscMagneticField3D(Rm,{Re/Re,0,0});

Bm  = B{4};
Bme = Be{4};

% Angular gyro frequency
Omegac = abs(qOverM)/gamma*Bm;

% Gyro period
tc = 2*pi/Omegac;

% Bouncing period
% factor Jupiter 1.2478 = Re/c*sqrt(2)*3.7
% factor Saturn  1.0521 = Re/c*sqrt(2)*3.7
fac = Re/c*sqrt(2)*3.7;
tb = fac*(R/Re)*c/Vpr*(1-0.4635*sind(alpha)^0.75);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
fprintf(1,'Bounce period Boris\n');%tb = R*sqrt(2)/Vpr*(3.7-1.6*sind(alpha));

% Drift period 
td = 2*pi*abs(qOverM)*Bme*Re^3/Vpr^2/R*(1-1/3*(sind(alpha))^0.62);
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

function [Xg,rg] = getGyroRadius(Rp,Vp)

tol = 1e-10;
maxIter = 100;

Rp = Rp(:);
Vp = Vp(:);

rg = zeros(3,1);
relerr = Inf;

% initial position of guiding centre
Xg = Rp;

fac = -gamma/qOverM;

iter = 0;

while relerr > tol && iter < maxIter

iter = iter + 1;

% magnetic field at guiding centre position
B = mdiscMagneticField3D(Rm,{Xg(1)/Re,Xg(2)/Re,Xg(3)/Re});

Bxyz = [B{1}; B{2}; B{3}];
Bm2 = B{4}.^2; % |B|^2


% parallel and perpendicular projection 
Vppar = (dot(Bxyz,Vp)/Bm2)*Bxyz;
Vpper = Vp-Vppar;

oldrg = rg;

% gyro radius rg = -gamma m/q (Vper x B)/B^2
rg = -fac * cross(Vpper,Bxyz)/Bm2;

delta2 = dot(rg-oldrg, rg-oldrg);
relerr = delta2 / max(dot(rg,rg), eps);

Xg = Rp+rg;

end

if iter == maxIter
  warning('Guiding-centre iteration did not converge.');
end

if 0
plot3(Xg(1)/Re,Xg(2)/Re,Xg(3)/Re,'x')
hold on
quiver3(Xg(1)/Re,Xg(2)/Re,Xg(3)/Re,...
        Bx/sqrt(Bm2),By/sqrt(Bm2),Bz/sqrt(Bm2),.2,'b.')
plot3(Rp(1)/Re,Rp(2)/Re,Rp(3)/Re,'x')
quiver3(Rp(1)/Re,Rp(2)/Re,Rp(3)/Re,Vp(1)/Vpr,Vp(2)/Vpr,Vp(3)/Vpr,.2,'b.')
%quiver3(X(1)/Re,X(2)/Re,X(3)/Re,X(4)/Vpr,X(5)/Vpr,X(6)/Vpr,.2,'r.')
hold off
xlim([.999*min(Rp(1),Xg(1)),1.001*max(Rp(1),Xg(1))]/Re)
%ylim([.99*min(Rp(2),Xg(2)),1.01*max(Rp(2),Xg(2))]/Re)
xlabel('x')
ylabel('y')
zlabel('z')
pause
end


end

end
