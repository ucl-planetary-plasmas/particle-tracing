function [B,varargout] = dipoleMagneticField3D(m, rm, r)
% function [B,gradB] = dipoleMagneticField3D(m, rm, r)
%
%    m      : magnetic moment vector
%    rm     : position of the magnetic moment
%    r      : cartesion coordinates where to calculate magnetic field
%    B      : dipole magnetic field (cell array of length(m), each cell is
%             an array of size(r))
%    gradB  : gradient of the dipole magnetic field (cell array of length(m),
%             each cell is an array of size(r))
%
% To plot the L=1 shell of a centred dipole with magnetic moment vector
% [0,0,1]
%
% clf
% Md = [0,0,1]; % magnetic moment aligned to z
% Rd = [0,0,0]; % no offset
% hold on
% for p=linspace(0,2*pi,7),
%   t = linspace(-pi/2,pi/2,100);
%   L = 1;
%   r = L *cos(t).^2;
%   ir = find(r>=0.1);
%   x = r(ir).*cos(t(ir))*cos(p);
%   y = r(ir).*cos(t(ir))*sin(p);
%   z = r(ir).*sin(t(ir));
%   [B,gradB]=dipoleMagneticField3D(Md,Rd,{x,y,z});
%   b=sqrt(B{1}.^2+B{2}.^2+B{3}.^2);
%   quiver3(x,y,z,B{1}./b,B{2}./b,B{3}./b,0.5)
%   db=sqrt(gradB{1}.^2+gradB{2}.^2+gradB{3}.^2);
%   quiver3(x,y,z,gradB{1}./db,gradB{2}./db,gradB{3}./db,0.5)
% end
% hold off
% xlabel('x'),ylabel('y'),zlabel('z')

%
% $Id: dipoleMagneticField3D.m,v 1.2 2017/10/15 20:19:27 patrick Exp $
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

% modulus square
R2 = X.^2 + Y.^2 + Z.^2;

% modulus
R = R2.^0.5;

% modulus cube
R3 = R2.^1.5;

% normalised direction
x = X./R; y = Y./R; z = Z./R;

% scalar product M.r
mDotR = m(1)*x + m(2)*y + m(3)*z;

% Bx, By and Bz
iR3 = 1./R3;
B{1} = iR3.*(3*mDotR.*x - m(1));
B{2} = iR3.*(3*mDotR.*y - m(2));
B{3} = iR3.*(3*mDotR.*z - m(3));

% |B|^2
B{4} = iR3.^2.*(3*mDotR.^2+sum(m.^2));
%B1=B{1}.^2+B{2}.^2+B{3}.^2;
%max(abs(B1-B{4}))
%plot(B1,B{4})

if 0, % test calculation of B in Cartesian from polar coordinates
Rcyl = sqrt(X.^2+Y.^2);
% colatitude
theta = pi/2-atan2(Z,Rcyl);
% cos and sin of colatitude
cost = Z./R;
sint = sqrt(1-cost.^2);
% B radial and tangential
Br = m(3)./R3*2.*cost;
Bt = m(3)./R3.*sint;
% cos and sin of latitude
coslat = Rcyl./R;
sinlat = Z./R;
% cos and sin of longitude
coslon = X./Rcyl;
sinlon = Y./Rcyl;
%b = sqrt(Br.^2+Bt.^2);
% B in polar to cartesian
Bx = (Br.*coslat+Bt.*sinlat).*coslon;
By = (Br.*coslat+Bt.*sinlat).*sinlon;
Bz = (Br.*sinlat-Bt.*coslat);
if 0,
disp('B components')
[B{1};B{3}]
[Bx;Bz]
subplot(211),plot([B{1};B{3}]'), legend({'B_x','B_z'})
subplot(212),plot([Bx;Bz]'), legend({'B_x','B_z'})
pause
end
end

if nargout>1,

% dBx/dx, dBy/dy, dBz/dz
% iR = 1./R;
% iR2 = 1./R2;
% fR6 = 3./R3.^2;
%gradB{1} = fR6.*((m(1).*x./R+mDotR.*(1-x.^2)./R2).*R3-(3*mDotR.*x-m(1)).*X.*R);
%gradB{2} = fR6.*((m(2).*y./R+mDotR.*(1-y.^2)./R2).*R3-(3*mDotR.*y-m(2)).*Y.*R);
%gradB{3} = fR6.*((m(3).*z./R+mDotR.*(1-z.^2)./R2).*R3-(3*mDotR.*z-m(3)).*Z.*R);

% d|B|/dx, d|B|/dy, d|B|/dz
fR6 = 3./R3.^2;
u = sqrt(3*mDotR.^2+sum(m.^2));
fac1 = R2.*mDotR./u;
fac2 = R.*u;
gradB{1} = fR6.*(fac1*m(1).*(1-x.^2) - X.*fac2);
gradB{2} = fR6.*(fac1*m(2).*(1-y.^2) - Y.*fac2);
gradB{3} = fR6.*(fac1*m(3).*(1-z.^2) - Z.*fac2);

varargout{1} = gradB;

%fprintf(1,'gradB(%.2f,%.2f,%.2f)=(%.2f,%.2f,%.2f)\n',X,Y,Z,gradB{:})
%pause

if 0, 
% test calculation of gradient in Cartesian from gradient in polar coordinates
EPS = 1e-6;
rp = R*(1+EPS);
rm = R*(1-EPS);
dr = 2*R*EPS;
% B  = m(3)./R3.*sqrt(1+3*cost.^2)
% gradient with respect to r
dBdr = m(3)*(1./rp.^3-1./rm.^3)./dr.*sqrt(1+3*cost.^2);
%
mup = cost*(1+EPS)+EPS; 
mum = cost*(1-EPS)-EPS;
rdmu = R.*(mup-mum);
dmudt = -sqrt(1-cost.^2);
% gradient with respect to theta
dBdt = m(3)./R3.*(sqrt(1+3*mup.^2)-sqrt(1+3*mum.^2))./rdmu.*dmudt;
gradBx = (dBdr.*coslat+dBdt.*sinlat).*coslon;
gradBy = (dBdr.*coslat+dBdt.*sinlat).*sinlon;
gradBz = (dBdr.*sinlat-dBdt.*coslat);
if 0,
disp('grad B components')
[gradB{1};gradB{2};gradB{3}]
[gradBx;gradBy;gradBz]
subplot(211), plot(Z,[gradB{1};gradB{3}]','-o'),
legend({'B_x','B_z'})
subplot(212), plot(Z,[gradBx;gradBz]','-o'),
legend({'B_x','B_z'})
pause
end
end

end
