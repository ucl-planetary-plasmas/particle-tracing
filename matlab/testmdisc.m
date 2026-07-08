function testmdisc(dipole_on)
% function testmdisc(dipole_on)
%
% dipole_on: optional dipole_on flag
%
% Compare magnetodisc field and equivalent dipolar field Cartesian
% components at the planet's surface.
% Plot both fields and their gradients along L shell contours
% for Jupiter (jup_mdisc_kh3e7_rmp90) 
% and Saturn  (sat_mdisc_kh2e6_rmp25)
%

%
% $Id: testmdisc.m,v 1.4 2026/07/08 17:42:39 patrick Exp $
%
% Copyright (c) 2017 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

clf

global MDFile DipoleOn

% Jupiter magnetodisc
MDFile = 'jup_mdisc_kh3e7_rmp90.mat';

if ~exist('dipole_on','var') | isempty(dipole_on)
  DipoleOn = false;
else
  DipoleOn = dipole_on;
end

% clear function with persistent variables
clear MDiscField
load(MDFile,'MD');

Re = MD.planet.r0;
% Jupiter dipole field
% Magnetic equator field strength in T
% B0 = 4.2800e-04 
% corrected to match magnetodisc file
% Bthdip(mu=0,r=1)=8.176230e-01 -> 3.499426e-04, Md=1.2787e20
% Bth   (mu=0,r=1)=8.172328e-01 -> 3.497756e-04, Md=1.2781e20
%B0 = MD.planet.B0 * MD.v2d.Bth(MD.dims.imu0,1);
B0 = MD.planet.B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
Md = [0,0,B0*Re^3]; % Magnetic moment
Rd = [0,0,0];       % Centered moment

plotmdisc()

pause

MDFile = 'sat_mdisc_kh2e6_rmp25.mat';
% clear function with persistent variables
clear MDiscField
load(MDFile,'MD');

Re = MD.planet.r0;
% Saturn dipole field
% Magnetic equator field strength in T
% B0 = 2.1160e-05;
%B0 = MD.planet.B0;
% corrected to match magnetodisc file
% Bthdip(mu=0,r=1)=8.923767e-01 -> 1.888269e-05, Md=4.1336e18
% Bth   (mu=0,r=1)=8.911721e-01 -> 1.885720e-05, Md=4.1280e18
%B0 = MD.planet.B0 * MD.v2d.Bth(MD.dims.imu0,1);
B0 = MD.planet.B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
Md = [0,0,B0*Re^3];  % Magnetic moment
Rd = [0,0,0];        % Centered moment

plotmdisc()

function plotmdisc()

fprintf(1,'M = %10e A m^2, Re = %10e m\n',Md(3),Re);

p = 0; % azimuth
t = linspace(-pi/2,pi/2,31); % latitude
%t = linspace(0,pi,31); % colatitude

% surface field
[B, gradB, curvB] = mdiscMagneticField3D(Rd,{1,0,0}); 
fprintf(1,'B0 mdisc  (%10e,%10e,%10e)\n', B{1:3})
[B, gradB, curvB] = dipoleMagneticField3D(Md,Rd,{Re,0,0});
fprintf(1,'B0 dipole (%10e,%10e,%10e)\n', B{1:3})
a = 0;
a = [-pi/3, -pi/6, 0, pi/6, pi/3];
cosa = cos(a); sina = sin(a);
r = {Re*cosa,zeros(size(cosa)),Re*sina};
[Bd,gradBd, curvBd] = dipoleMagneticField3D(Md,Rd,r);

%return

r = 1;
x = r*cos(t)*cos(p);
y = r*cos(t)*sin(p);
z = r*sin(t);
B  = mdiscMagneticField3D(Rd,{x,y,z});
Bd = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});

xc = x; yc = y; zc = z;

plot(t,[B{1};Bd{1};B{2};Bd{2};B{3};Bd{3}])
xlabel('\theta');legend({'B_x','B_{xd}','B_y','B_{yd}','B_z','B_{zd}'})

lut = get(gca,'colororder');

pause

clf
hold on
for L=[2,5,10],
  r = L *cos(t).^2;
  ir = find(r>=1);
  x = r(ir).*cos(t(ir))*cos(p);
  y = r(ir).*cos(t(ir))*sin(p);
  z = r(ir).*sin(t(ir));

	% mdisc
  [B,gradB,curvB] = mdiscMagneticField3D(Rd,{x,y,z});
  Bm = B{4};
	gradBm = gradB{4};
	Km= curvB{4};
  quiver(x,z,B{1}./Bm,B{3}./Bm,0.5,'color',lut(1,:));
	quiver(x,z,gradB{1}./gradBm,gradB{3}./gradBm,0.5,'Color',lut(2,:));
	quiver(x,z,curvB{1}./Km,curvB{3}./Km,0.5,'Color',lut(3,:));

  % dipole
  [Bd,gradBd,curvBd] = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});
  Bdm = Bd{4};
	gradBdm = gradBd{4};
	Kdm= curvBd{4};
  quiver(x,z,Bd{1}./Bdm,Bd{3}./Bdm,0.5,'Color',lut(4,:));
	quiver(x,z,gradBd{1}./gradBdm,gradBd{3}./gradBdm,0.5,'Color',lut(5,:));
	quiver(x,z,curvBd{1}./Kdm,curvBd{3}./Kdm,0.5,'Color',lut(6,:));
end
plot(xc,zc,'k-','LineWidth',2), axis equal
hold off
xlabel('x');
ylabel('z'); 
legend({'B_m','\nabla{B_m}','\kappa_m','B_d','\nabla{B_d}','\kappa_d'});

end

end
