function [B,varargout] = mdiscMagneticField3D(rm, r)
% function [B,gradB] = mdiscMagneticField3D(rm, r)
%
%    rm     : position of the magnetic moment
%    r      : cartesion coordinates where to calculate magnetic field in RJ
%    B      : magnetodisc magnetic field (cell array of length(m), each cell is
%             an array of size(r))
%    gradB  : gradient of the magnetodisc magnetic field (cell array of
%             length(m), each cell is an array of size(r))
%
%
% Example to compare Jupiter's magnetodisc to a dipole
%
% global MDFile
% MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
% % clear function with persistent variables
% clear MDiscField
% load(MDFile,'MD')
%
% % equatorial radius in m
% Re = MD.planet.r0;
%
% % magnetic equator field strength in T
% B0 = MD.planet.B0;
% % corrected to match dipole value in the magnetodisc file
% B0 = MD.planet.B0 * MD.v2d.Bthdip(MD.dims.imu0,1);
% Md = [0,0,B0*Re^3];  % Magnetic moment
% Rd = [0,0,0];        % Centered moment
% 
% p = 0;
% t = linspace(-pi/2,pi/2,100);
%
% % unit half circle
% x = cos(t)*cos(p);
% y = cos(t)*sin(p);
% z = sin(t);
% B  = mdiscMagneticField3D(Rd,{x,y,z});
% Bd = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});
% xc = x; yc = y; zc = z;
%
% subplot(211),
% plot(t,[B{1};Bd{1};B{3};Bd{3}])
% xlabel('\theta');legend({'B_x','B_{xd}','B_z','B_{zd}'})
% lut = get(gca,'colororder');
%
% subplot(212)
% hold on
% for L=[2,5,10],
%   r = L *cos(t).^2;
%   ir = find(r>=1);
%   x = r(ir).*cos(t(ir))*cos(p);
%   y = r(ir).*cos(t(ir))*sin(p);
%   z = r(ir).*sin(t(ir));
%   [B,gradB] = mdiscMagneticField3D(Rd,{x,y,z});
%   b = sqrt(B{1}.^2+B{2}.^2+B{3}.^2);
%   db = sqrt(gradB{1}.^2+gradB{2}.^2+gradB{3}.^2);
%   quiver(x,z,B{1}./b,B{3}./b,0.5,'color',lut(1,:));
%   quiver(x,z,gradB{1}./db,gradB{3}./db,0.5,'color',lut(2,:));
%   [Bd,gradBd] = dipoleMagneticField3D(Md,Rd,{Re*x,Re*y,Re*z});
%   bd = sqrt(Bd{1}.^2+Bd{2}.^2+Bd{3}.^2);
%   dbd = sqrt(gradBd{1}.^2+gradBd{2}.^2+gradBd{3}.^2);
%   quiver(x,z,Bd{1}./bd,Bd{3}./bd,0.5,'color',lut(3,:));
%   quiver(x,z,gradBd{1}./dbd,gradBd{3}./dbd,0.5,'color',lut(4,:));
% end
% plot(xc,zc,'k-','LineWidth',2), axis equal
% hold off
% xlabel('x'); ylabel('z'); legend({'B','\nabla B','Bd','\nabla Bd'});

%
% $Id: mdiscMagneticField3D.m,v 1.4 2017/10/18 07:50:02 patrick Exp $
%
% Copyright (c) 2009-2016 Patrick Guio <patrick.guio@gmail.com>
%
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

% \vec{B}(\vec{m},\vec{r} ) = \frac{\mu_0}{4\pi r^3}
%    \left(3(\vec{m}\cdot\vec{\hat{r}}) \vec{\hat{r}} -\vec{m} \right) +
%    \frac{2\mu_0}{3} \vec{m} \delta^3(\vec{r})

% Separate coordinates R=(X,Y,Z)
[X, Y, Z] = deal(r{:});

% Subtract moment position
X = X-rm(1);
Y = Y-rm(2);
Z = Z-rm(3);

% R modulus
R = sqrt(X.^2 + Y.^2 + Z.^2);

% Projection on (x,y) plane
Rcyl = sqrt(X.^2+Y.^2);

% colatitude in radians
theta = pi/2-atan2(Z,Rcyl);

% cos and sin of latitude
coslat = Rcyl./R; 
sinlat = Z./R; 
%sqrt(coslat^2+sinlat^2)

% cos and sin of longitude
coslon = X./Rcyl; 
sinlon = Y./Rcyl; 
%sqrt(coslon^2+sinlon^2)

if nargout<=1,
  [Br, Bt] = MDiscField(R,theta);
else
  [Br, Bt, gradB] = MDiscField(R,theta);
end

B2 =  Br.^2+Bt.^2;

Bx = (Br.*coslat+Bt.*sinlat).*coslon;
By = (Br.*coslat+Bt.*sinlat).*sinlon;
Bz = (Br.*sinlat-Bt.*coslat);

B = {Bx,By,Bz,B2};

if nargout>1,

dBr = gradB{1};
dBt = gradB{2};

gradBx = (dBr.*coslat+dBt.*sinlat).*coslon;
gradBy = (dBr.*coslat+dBt.*sinlat).*sinlon;
gradBz = (dBr.*sinlat-dBt.*coslat);

varargout{1} = {gradBx,gradBy,gradBz};

end
