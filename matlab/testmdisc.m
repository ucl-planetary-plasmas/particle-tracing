function testmdisc
% funnction testmdisc
%
% Compare magnetodisc field and equivalent dipolar field Cartesian
% components at the planet's surface.
% Plot both fields and their gradients along L shell contours
% for Jupiter (jup_mdisc_kh3e7_rmp90) 
% and Saturn  (sat_mdisc_kh2e6_rmp25)
%

%
% $Id: testmdisc.m,v 1.3 2017/10/18 08:01:55 patrick Exp $
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

global MDFile

% Jupiter magnetodisc
MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
% clear function with persistent variables
clear MDiscField
load(MDFile,'MD');

% Jupiter equatorial radius in m
%Re = 71492.0e3; 
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

% Saturn magnetodisc
global MDFile
MDFile = 'sat_mdisc_kh2e6_rmp25.mat';

plotmdisc()

function plotmdisc()

fprintf(1,'%10e\n',Md(3));

p = 0;
t = linspace(-pi/2,pi/2,31);

[B, gradB]=mdiscMagneticField3D(Rd,{1,0,0}); 
fprintf(1,'B0 mdisc  (%10e,%10e,%10e)\n', B{1:3})
[B, gradB]=dipoleMagneticField3D(Md,Rd,{Re,0,0});
fprintf(1,'B0 dipole (%10e,%10e,%10e)\n', B{1:3})
a = 0;
a = [-pi/3, -pi/6, 0, pi/6, pi/3];
cosa = cos(a); sina = sin(a);
r = {Re*cosa,zeros(size(cosa)),Re*sina};
[Bd,gradBd]=dipoleMagneticField3D(Md,Rd,r);

%return

r = 1;
x = r*cos(t)*cos(p);
y = r*cos(t)*sin(p);
z = r*sin(t);
B  = mdiscMagneticField3D(Rd,{x,y,z});
Bd = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});

xc = x; yc = y; zc = z;

plot(t,[B{1};Bd{1};B{3};Bd{3}])
xlabel('\theta');legend({'B_x','B_{xd}','B_z','B_{zd}'})

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
  [B,gradB] = mdiscMagneticField3D(Rd,{x,y,z});
  b = sqrt(B{1}.^2+B{2}.^2+B{3}.^2);
	db = sqrt(gradB{1}.^2+gradB{2}.^2+gradB{3}.^2);
  quiver(x,z,B{1}./b,B{3}./b,0.5,'color',lut(1,:));
	quiver(x,z,gradB{1}./db,gradB{3}./db,0.5,'Color',lut(2,:));
  [Bd,gradBd] = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});
  bd=sqrt(Bd{1}.^2+Bd{2}.^2+Bd{3}.^2);
	dbd = sqrt(gradBd{1}.^2+gradBd{2}.^2+gradBd{3}.^2);
  quiver(x,z,Bd{1}./bd,Bd{3}./bd,0.5,'Color',lut(3,:));
	quiver(x,z,gradBd{1}./dbd,gradBd{3}./dbd,0.5,'Color',lut(4,:));
end
plot(xc,zc,'k-','LineWidth',2), axis equal
hold off
xlabel('x');ylabel('z'); legend({'B','\nabla B','Bd','\nabla Bd'});

end

end
