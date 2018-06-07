function cmptraj_earth(timespec)
% function cmptraj_earth(timespec)
%
% Parameterization of Ozturk's 2012 paper 2012
% Trajectories of charged particles trapped in Earth's magnetic field
%
% timespec is such that tmax = sum(timespec.*[1,tc,tb,td])
%
% For instance 
%
%   cmptraj_earth([0,0,2,0]);
%
% runs the comparison between full dynamic using ODE solver and the
% Boris algorithm, and guiding centre approximantion for 2xtb, i.e.
% two bounce periods.

%
% $Id: cmptraj_earth.m,v 1.5 2017/10/17 20:20:55 patrick Exp $
%
% Copyright (c) 2009-2016 Patrick Guio <patrick.guio@gmail.com>
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
solver = 'ode45';
solver = 'ode15s';
%solver = 'ode23t';
% ODE solver options
options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','off','vectorized','on');
%options = odeset('AbsTol',1e-10,'RelTol',1e-10,'Stats','off','vectorized','on');
%options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Stats','on','vectorized','on');
%options = odeset('AbsTol',1e-8,'RelTol',1e-6,'Stats','on','vectorized','on');

% close all figures
close all

Re = 6378e3; % Earth equatorial radius in m

% Earth dipole field
B0 = 3.07e-5;        % Magnetic equator field strength in T
Md = [0,0,-B0*Re^3]; % Magnetic moment
%Md = [0,0,-7.9651e15];
Rm = [0,0,0];        % Centered moment

% SI constants
mp = 1.6726231e-27;  % proton mass 
qp = 1.60217733e-19; % proton charge
eV = 1.60217733e-19; % conversion factor eV to J
c  = 299792458.0;    % speed of light

Ep = 10e6*eV;       % 10 MeV proton
vp = sqrt(2/mp*Ep); % classic velocity 
fprintf('classic v=%.4g m/s E=%.4g MeV\n',vp,.5*mp*vp^2/(1e6*eV))
vp = c*sqrt(Ep*(Ep+2*mp*c^2))/(Ep+mp*c^2); % relativistic velocity
fprintf('relativ v=%.4g m/s E=%.4g MeV\n',vp,(1/sqrt(1-vp^2/c^2)-1)*mp*c^2/(1e6*eV));
% velocity squared for GC solution
V2 = vp^2;

% relativistic factor
gamma = 1/sqrt(1-V2/c^2);

% proton charge-to-mass ratio
qOverM = qp/mp;

% pitch angle deg
alpha = 30;
alpha = 15;
%alpha = 4;

close all
% initial conditions
[X,Xgc,mu,tc,tb,td,lm] = init(4*Re,vp,alpha);
% factor mu/gamma^2/m for guiding centre Eq.23 of Ozturk 2012
facdv = mu/gamma^2/mp; 
tspan = [0,sum(timespec.*[1,tc,tb,td])];
fprintf(1,'tmax=%.2f (%.1f tc, %.2f tb, %.2f td)\n',...
        tspan(end)./[1,tc,tb,td]);
trace(X,Xgc);

function [Xo,Xogc,mu,tc,tb,td,lm]=init(R,v,alpha)

% initial condition x0,y0,z0,vx0,vy0,vz0
Xo = [R;            0;            0;...
      0;v*sind(alpha);v*cosd(alpha)];

% initial guiding centre x0,y0,z0,vpar0
B = dipoleMagneticField3D(Md,Rm,{Xo(1),Xo(2),Xo(3)});
b = sqrt(B{4});
% smaller vpar gives larger vper thus guiding centre motion
%vpar = 0.9220*sum(Xo(4:6).*[B{1};B{2};B{3}]/b);
vpar = sum(Xo(4:6).*[B{1};B{2};B{3}]/b);
vper = sqrt(V2-vpar^2);
% gyro radius
rho = abs(gamma/qOverM*vper/b);
% add or substract gyroradius depending on sign of q and magnetic moment
Xogc = [R-sign(qOverM*Md(3))*rho;0;0;vpar];
B = dipoleMagneticField3D(Md,Rm,{Xogc(1),Xogc(2),Xogc(3)});
b = sqrt(B{4});
% first invariant mu
mu = gamma^2*mp*vper^2/(2*b);
[tc,tb,td] = periods(R);
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
[t,X] = solverh(@dynamic3d,tspan,Xo,options); 
toc

Z    = X(:,3);
Rcyl = sqrt(X(:,1).^2+X(:,2).^2);
Rtot = sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);

tic
rm = {Xo(1)*cosd(lm)^3,Xo(1)*cosd(lm)^2*sind(lm),0};
B = dipoleMagneticField3D(Md,Rm,rm);
% Gyro period at mirror point
tm = 2*pi/(qOverM/gamma*sqrt(B{4}));
%tb = linspace(0,t(end),floor(t(end)/max(diff(t))));
tb = linspace(tspan(1),tspan(end),4*fix(diff(tspan)/tm));
%tb = t;
hb = diff(tb);
fprintf(1,'tm=%.2f, tc/tm=%.2f, tc/hb(1)=%.2f\n', [tm,fix(tc/tm),fix(tc/hb(1))]);
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
[tgc,Xgc] = solverh(@gcdynamic3d,tspan,Xogc,options);
toc

Zgc    = Xgc(:,3);
Rcylgc = sqrt(Xgc(:,1).^2+Xgc(:,2).^2);
Rtotgc = sqrt(Xgc(:,1).^2+Xgc(:,2).^2+Xgc(:,3).^2);


figure
plot(Rcyl/Re,Z/Re,Rcylb/Re,Zb/Re,Rcylgc/Re,Zgc/Re,Xp,Yp), 
axis equal,
xlabel('rcyl'), ylabel('z')
legend({'FD','FDB','GC'})
fprintf(1,'Ok, press return\n'), 
pause

% save all trajectories
save cmptraj_earth t X tb Xb tgc Xgc timespec

% Kinetic energy in eV for Full Dynamic
v2 = sum(X(:,[4:6])'.^2);
%E = 1/2*mp*v2/eV;  % classic
E = (1./sqrt(1-v2/c^2)-1)*mp*c^2/eV; % relativistic
% mean and variance
meanE = mean(E);
stdE = std(E);
fprintf(1,'<E>  = %10g MeV, std(E) = %10g eV\n', meanE/1e6, stdE);

% Kinetic energy in eV for Full Dynamic Boris
v2b = sum(Xb(:,[4:6])'.^2);
%E = 1/2*mp*v2b/eV;  % classic
Eb = (1./sqrt(1-v2b/c^2)-1)*mp*c^2/eV; % relativistic
% mean and variance
meanEb = mean(Eb);
stdEb = std(Eb);
fprintf(1,'<Eb> = %10g MeV, std(E) = %10g eV\n', meanEb/1e6, stdEb);

plot(t,E-meanE,tb,Eb-meanEb);
xlabel('time'), ylabel('E-<E>')
legend({'FD','FDB'})
fprintf(1,'Ok, press return\n'), 
pause

% Compute b dot gradB and mu for Full Dynamic
[B,gradB] = dipoleMagneticField3D(Md, Rm, {X(:,1),X(:,2),X(:,3)});
BgB = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./sqrt(B{4});
Vper2 = V2-(B{1}.*X(:,4)+B{2}.*X(:,5)+B{3}.*X(:,6)).^2./B{4};
% Instantaneous first invariant Eq.14 
mui = gamma^2*mp*Vper2./(2*sqrt(B{4}));

% Compute b dot gradB and mu for Full Dynamic Boris
[B,gradB] = dipoleMagneticField3D(Md, Rm, {Xb(:,1),Xb(:,2),Xb(:,3)});
BgBb = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./sqrt(B{4});
Vper2b = V2-(B{1}.*Xb(:,4)+B{2}.*Xb(:,5)+B{3}.*Xb(:,6)).^2./B{4};
% Instantaneous first invariant Eq.14
muib = gamma^2*mp*Vper2b./(2*sqrt(B{4}));

% Compute b dot gradB and mu for Guiding Centre
[B,gradB] = dipoleMagneticField3D(Md, Rm, {Xgc(:,1),Xgc(:,2),Xgc(:,3)});
BgBgc = (B{1}.*gradB{1}+B{2}.*gradB{2}+B{3}.*gradB{3})./sqrt(B{4});
Vper2gc = V2-Xgc(:,4).^2;
% Instantaneous first invariant Eq.14 
muigc = gamma^2*mp*Vper2gc./(2*sqrt(B{4}));

subplot(211), 
plot(t,mui,tb,muib,tgc,muigc,tgc,facdv*gamma^2*mp*ones(size(tgc))),
ylabel('\mu 1st invariant')
legend('FD','FDB','GC','Initial')

% Compute latitude and longitude
lat = atan2d(Z,Rcyl);        % atan(Rcyl/Z)
lon = atan2d(X(:,2),X(:,1)); % atan(X/Y);
latb = atan2d(Zb,Rcylb);        % atan(Rcyl/Z)
lonb = atan2d(Xb(:,2),Xb(:,1)); % atan(X/Y);
latgc = atan2d(Zgc,Rcylgc);
longc = atan2d(Xgc(:,2),Xgc(:,1));

subplot(212), 
%plot(t,X(:,3)/Re,tgc,Xgc(:,3)/Re), ylabel('z')
plot(t,lat,tb,latb,tgc,latgc), xlabel('time'); ylabel('Latitude')
legend('FD','FDB','GC')

fprintf(1,'Ok, press return\n'), 
pause
%subplot(211), plot(t,BgB,tgc,BgBgc),
%subplot(212), plot(X(:,3)/Re,BgB,Xgc(:,3)/Re,BgBgc),pause

% L-Shell with L
t = linspace(-pi/2,pi/2,100);
% estimate L from instantaneous L
%coslat = Rcyl./R; Ls = R./coslat.^2/Re;
Ls = Rtot.^3./Rcyl.^2/Re;
Lsgc = Rtotgc.^3./Rcylgc.^2/Re;
Le = mean(Ls);
fprintf(1,'L estimated = %.2f\n', Le);
if 0,
	f=figure; 
	plot(t,Ls), 
	pause, close(f);
end
r = Le *cos(t).^2;
ir = find(r>=1);
xe = r(ir).*cos(t(ir));
ye = r(ir).*sin(t(ir));

% plot trajectory Rcyl vs Z
subplot(211), 
plot(Rcyl/Re,Z/Re,Rcylb/Re,Zb/Re,Rcylgc/Re,Zgc/Re,xe,ye,'-.',Xp,Yp), 
xlabel('$r_{cyl}=\sqrt{x^2+y^2}$','interpreter','latex'),
ylabel('$z$','interpreter','latex'),
title(sprintf('%s x0=%.0f x0gc=%.4f L=%.4f', ...
      solver, Xo(1)/Re, Xogc(1)/Re,Le),...
      'interpreter','latex')
legend({'FD','FDB','GC','L-shell'},'interpreter','latex');
%set(gca,'xlim',[0 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal,

subplot(212),
if 1,
% plot trajectory X vs Y
plot(X(:,1)/Re,X(:,2)/Re,...
     Xb(:,1)/Re,Xb(:,2)/Re,...
     Xgc(:,1)/Re,Xgc(:,2)/Re,Xp,Yp), 
xlabel('$x$','interpreter','latex'),
ylabel('$y$','interpreter','latex'),
else
% 3d plot trajectory
plot3(X(:,1)/Re,X(:,2)/Re,X(:,3)/Re,...
      Xb(:,1)/Re,Xb(:,2)/Re,Xb(:,3)/Re,...
      Xgc(:,1)/Re,Xgc(:,2)/Re,Xgc(:,3)/Re), 
xlabel('$x$','interpreter','latex'),
ylabel('$y$','interpreter','latex'),
zlabel('$z$','interpreter','latex'),
end
title(sprintf('$\\langle{E}\\rangle=$%2g MeV std$(E)=$%2g eV',meanE/1e6,stdE),...
      'interpreter','latex');
%set(gca,'xlim',[-4.5 4.5])
%set(gca,'ylim',[-4.5 4.5])
axis equal

drawnow

end

function rn = BorisIter(r,h)

rn = zeros(size(r));

[x,y,z,vx,vy,vz] = deal(r(1),r(2),r(3),r(4),r(5),r(6));

B = dipoleMagneticField3D(Md, Rm, {x,y,z});
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

B = dipoleMagneticField3D(Md, Rm, {x,y,z});
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

% 3D ODE function dr/dt=f(r) where r=(x,y,z,vx,vy,vz)
function dr_dt = dynamic3d(t,r)

%fprintf(1,'t=%.2g\n', t)

[x,y,z,vx,vy,vz] = deal(r(1,:),r(2,:),r(3,:),r(4,:),r(5,:),r(6,:));

B = dipoleMagneticField3D(Md, Rm, {x,y,z});
[Bx,By,Bz] = deal(B{1:3});

fac = qOverM/gamma;
dvx = fac*(vy.*Bz - vz.*By);
dvy = fac*(vz.*Bx - vx.*Bz);
dvz = fac*(vx.*By - vy.*Bx);

dr_dt = [vx;vy;vz;dvx;dvy;dvz];

end

function dr_dt = gcdynamic3d(t,r)

[x,y,z,vpar] = deal(r(1,:),r(2,:),r(3,:),r(4,:));

[B,gradB] = dipoleMagneticField3D(Md, Rm, {x,y,z});
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

function [tc,tb,td] = periods(R)

B = dipoleMagneticField3D(Md, Rm, {R,0,0});
Be = dipoleMagneticField3D(Md, Rm, {Re,0,0});

B  = sqrt(B{4});
Be = sqrt(Be{4});

% Angular gyro frequency
Omegac = qOverM/gamma*B;

% Gyro period
tc = 2*pi/Omegac;

% Bouncing period
% factor 0.117 = Re/c*sqrt(2)*3.7
tb = 0.117*(R/Re)*(c/vp)*(1-0.4635*sind(alpha)^0.75);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
%tb = R*sqrt(2)/vp*(3.7-1.6*sind(alpha));

% Drift period 
td = 2*pi*qOverM*Be*Re^3/vp^2/R*(1-1/3*(sind(alpha))^0.62);
% https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
%td = pi*qp*Be*Re^3/(3*Ep*R)/(0.35+0.15*sind(alpha));
%td = 1.05/(Ep/(1e6*eV))/(R/Re)/(1+0.43*sind(alpha))*3600;

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
