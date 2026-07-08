function [B,varargout] = dipoleMagneticField3D(M, rm, r)
% function [B,gradB,curvB] = dipoleMagneticField3D(M, rm, r)
%
%    M      : magnetic moment vector [Mx,My,Mz]
%    rm     : position offset of the magnetic moment [xm,ym,zm]
%    r      : Cartesian coordinates where to calculate magnetic field [Rp]
%    B = {Bx, By, Bz, Bm}
%           : dipole magnetic field 
%             (cell array of length(m), each cell is an array of size(r))
%             {[T], [T], [T], [T]}
%    gradB = {dBm/dx, dBm/dy, dBm/dz, |grad Bm|, LB=Bm/|grad Bm|}
%           : optional magnetic field gradient parameters
%             (cell array of length(m), each cell is an array of size(r))
%            {[T/m], [T/m], [T/m], [T/m], [m]}
%    curvB = {Kx, Ky, Kz, Km, Rc=1/Km}
%           : optional magnetic field curvature parameters 
%             (cell array of length(m), each cell is an array of size(r))
%            {[m-1], [m-1], [m-1], [m-1], [m]}
%
% To plot the L=1 shell of a centred dipole with M = [0,0,1]
%
% clf
% Md = [0,0,1]; % magnetic moment along z-axis
% Rd = [0,0,0]; % no offset
% cols = get(gca,'ColorOrder');
% c1 = cols(1,:);
% c2 = cols(2,:);
% c3 = cols(3,:);
% t = linspace(-pi/2,pi/2,51);
% L = 1;
% r = L *cos(t).^2;
% ir = find(r>=0.1);
% r = r(ir); 
% t = t(ir);
% hold on
% for p=linspace(0,3*pi/2,4),
%   x = r.*cos(t)*cos(p);
%   y = r.*cos(t)*sin(p);
%   z = r.*sin(t);
%   [B,gradB,curvB]=dipoleMagneticField3D(Md,Rd,{x,y,z});
%   Bm=B{4};
%   quiver3(x,y,z,B{1}./Bm,B{2}./Bm,B{3}./Bm,.5,'MaxHeadSize',1,'Color',c1)
%   gradBm=gradB{4};
%   quiver3(x,y,z,gradB{1}./gradBm,gradB{2}./gradBm,gradB{3}./gradBm,1,'MaxHeadSize',.5,'Color',c2)
%   Km=curvB{4};
%   quiver3(x,y,z,curvB{1}./Km,curvB{2}./Km,curvB{3}./Km,.5,'MaxHeadSize',1,'Color',c3)
% end
% hold off
% legend({'B','\nabla{B}','\kappa'})
% xlabel('x'),ylabel('y'),zlabel('z')
% view(30,45)
% axis equal

%
% $Id: dipoleMagneticField3D.m,v 1.4 2026/07/08 17:51:32 patrick Exp $
%
% Copyright (c) 2009 Patrick Guio <patrick.guio@gmail.com>
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

%Rxyz = [r{1}-rm(1), r{2}-rm(2), r{3}-rm(3)];

% R^2 = R.R
R2 = X.^2 + Y.^2 + Z.^2;
% |R|
R = sqrt(R2);
% |R|^3
R3 = R2.*R;
% |R|^4
R4 = R2.*R2;

% 1/R^2
R_4 = 1./R4;
% 1/|R|^5
R_5 = 1./R.^5;
% 1/|R|^8
R_8 = R_4.^2;

% normalised direction r = R/|R|
%x = X./R; y = Y./R; z = Z./R;
% scalar product M.r
%Mdotr = M(1)*x + M(2)*y + M(3)*z;

% scalar product \alpha = M.R
alpha = M(1)*X + M(2)*Y + M(3)*Z;
% \alpha^2
alpha2 = alpha.^2;
% M^2 = M.M
M2 = dot(M, M);

% Bx, By and Bz
%iR3 = 1./R3;
%Bx = iR3.*(3*Mdotr.*x - M(1));
%By = iR3.*(3*Mdotr.*y - M(2));
%Bz = iR3.*(3*Mdotr.*z - M(3));
%[Bx,By,Bz]
Bx = R_5.*(3*alpha.*X - M(1)*R2);
By = R_5.*(3*alpha.*Y - M(2)*R2);
Bz = R_5.*(3*alpha.*Z - M(3)*R2);
%[Bx,By,Bz]

% |B|
%Bm = iR3.*sqrt(3*Mdotr.^2+sum(M.^2));
%Bm
% S = 3\alpha^2+R^2 M^2
S = 3*alpha2 + R2*M2;
S1_2 = sqrt(S);
Bm = R_4.*S1_2;
%Bm

B = {Bx, By, Bz, Bm};

%B1 = sqrt(B{1}.^2+B{2}.^2+B{3}.^2);
%max(abs(B1-B{4}))
%plot(B1,B{4})

if 0, % test B calculation for M along z in polar coordinates
Rcyl = sqrt(X.^2+Y.^2);
% colatitude
theta = pi/2-atan2(Z,Rcyl);
% cos and sin of colatitude
cost = Z./R;
sint = sqrt(1-cost.^2);
% B radial and tangential components
Br = M(3)./R3*2.*cost;
Bt = M(3)./R3.*sint;
% cos and sin of latitude
coslat = Rcyl./R;
sinlat = Z./R;
% cos and sin of longitude
coslon = X./Rcyl;
sinlon = Y./Rcyl;
% B in polar to Cartesian coordinates
Bx = (Br.*coslat+Bt.*sinlat).*coslon;
By = (Br.*coslat+Bt.*sinlat).*sinlon;
Bz = (Br.*sinlat-Bt.*coslat);
if 0,
disp('B components')
[B{1};B{2};B{3}]
[Bx;By;Bz]
[B{1};B{2};B{3}]-[Bx;By;Bz]
subplot(211),plot([B{1};B{3}]'), legend({'B_x','B_z'})
subplot(212),plot([Bx;Bz]'), legend({'B_x','B_z'})
pause
end
end

if nargout>1,

% d|B|/dx, d|B|/dy, d|B|/dz
%fR6 = 3./R3.^2;
%u = sqrt(3*Mdotr.^2+sum(M.^2));
%fac1 = R2.*Mdotr./u;
%fac2 = R.*u;
%gradBx = fR6.*(fac1*M(1).*(1-x.^2) - X.*fac2);
%gradBy = fR6.*(fac1*M(2).*(1-y.^2) - Y.*fac2);
%gradBz = fR6.*(fac1*M(3).*(1-z.^2) - Z.*fac2);

% 1/|R|^6
R_6 = 1./R3.^2;
% 1/sqrt(S)
S_1_2 = 1./S1_2;
% d|B|/dx, d|B|/dy, d|B|/dz
gradBx = R_6.*(R2.*(M2*X+3*alpha.*M(1))-4*S.*X).*S_1_2;
gradBy = R_6.*(R2.*(M2*Y+3*alpha.*M(2))-4*S.*Y).*S_1_2;
gradBz = R_6.*(R2.*(M2*Z+3*alpha.*M(3))-4*S.*Z).*S_1_2;
% |grad|B||
gradBm = sqrt(gradBx.^2+gradBy.^2+gradBz.^2);
% |B|/|grad|B||
LB = Bm./gradBm;

gradB = {gradBx, gradBy, gradBz, gradBm, LB};
varargout{1} = gradB;

%fprintf(1,'gradB(%.2e,%.2e,%.2e)=(%.2e,%.2e,%.2e)\n',X,Y,Z,gradB{1:3})
%pause

if 0, % test grad B calculation for M along z in polar coordinates
EPS = 1e-6;
rp = R*(1+EPS);
rm = R*(1-EPS);
dr = 2*R*EPS;
% B  = M(3)./R3.*sqrt(1+3*cost.^2)
Bm = @(r,mu) M(3)./r.^3.*sqrt(1+3*mu.^2);
% gradient with respect to r
dBdr = (Bm(rp,cost)-Bm(rm,cost))./dr;
%
mup = cost*(1+EPS)+EPS; 
mum = cost*(1-EPS)-EPS;
rdmu = R.*(mup-mum);
dmudt = -sqrt(1-cost.^2);
% gradient with respect to theta 1/r dB/dmu dmu/dt
dBdt = (Bm(R,mup)-Bm(R,mum))./rdmu.*dmudt;
gradBx1 = (dBdr.*coslat+dBdt.*sinlat).*coslon;
gradBy1 = (dBdr.*coslat+dBdt.*sinlat).*sinlon;
gradBz1 = (dBdr.*sinlat-dBdt.*coslat);
if 0,
disp('grad B components')
[gradBx;gradBy;gradBz]
[gradBx1;gradBy1;gradBz1]
[gradBx;gradBy;gradBz]-[gradBx1;gradBy1;gradBz1]
subplot(211), plot(Z,[gradBx;gradBz]','-o'),
legend({'dB/dx','dB/dz'})
subplot(212), plot(Z,[gradBx1;gradBz1]','-o'),
legend({'dB/dx','dB/dz'})
pause
end
end

if nargout>2,

% Q = 3 (M.R) R - (R.R) M
Qx = 3*alpha.*X-R2*M(1);
Qy = 3*alpha.*Y-R2*M(2);
Qz = 3*alpha.*Z-R2*M(3);

% X = (18 \alpha^2-3 R^2 M^2)R -7 alpha R^2 M
Xx = (18*alpha2-3*R2*M2).*X-7*alpha.*R2*M(1);
Xy = (18*alpha2-3*R2*M2).*Y-7*alpha.*R2*M(2);
Xz = (18*alpha2-3*R2*M2).*Z-7*alpha.*R2*M(3);

% Q.X = \alpha R^2 (15 \alpha^2 + R^2 M^2)
QdotX = alpha.*R2.*(15*alpha2+R2*M2);

%Q2 = Q.Q = Qx.^2+Qy.^2+Qz.^2;
Q2 = R2.*(R2*M2+3*alpha2);
% (Q.Q)^2
Q4 = Q2.^2;

fac1 = 1./Q2;
fac2 = -QdotX./Q4;
% K = 1/Q^2 * (X - (Q.X)/Q^2 Q) = 1/Q^2 * (I-bb) X
curvBx = Xx.*fac1+fac2.*Qx;
curvBy = Xy.*fac1+fac2.*Qy;
curvBz = Xz.*fac1+fac2.*Qz;

% K = (9\alpha^2+6\alpha^2 R^2 M^2 -3 R^4 M^4)/R^2/(R^2 M^2+3\alpha^2 R
%     -6\alpha R^2(R^2 M^2+\alpha^2)/R^2/(R^2 M^2+3\alpha^2 M
%fac1 = (9*alpha2.^2+6*alpha2.*R2*M2-3*R4*M2^2)./R2./(R2*M2+3*alpha2).^2;
%fac2 = -6*alpha.*R2.*(R2*M2+alpha2)./R2./(R2*M2+3*alpha2).^2;
%curvBx = fac1.*X+fac2.*M(1);
%curvBy = fac1.*Y+fac2.*M(2);
%curvBz = fac1.*Z+fac2.*M(3);

% |K|
curvBm = sqrt(curvBx.^2+curvBy.^2+curvBz.^2);

% R_c = 1/|K|
Rc = 1./curvBm;

curvB = {curvBx, curvBy, curvBz, curvBm, Rc};
varargout{2} = curvB;

if 0, % test curvature calculation for M along z in polar coordinates
EPS = 1e-6;
rp = R*(1+EPS);
rm = R*(1-EPS);
dr = 2*R*EPS;
br = @(mu) 2*mu./sqrt(1+3*mu.^2);
bt  = @(mu) sqrt(1-mu.^2)./sqrt(1+3*mu.^2);
% gradient with respect to r
dbrdr = zeros(size(R));
dbtdr = zeros(size(R));

%
mup = cost*(1+EPS)+EPS;
mum = cost*(1-EPS)-EPS;
rdmu = R.*(mup-mum);
dmudt = -sqrt(1-cost.^2);
% gradient with respect to theta
dbrdt = (br(mup)-br(mum))./rdmu.*dmudt;
dbtdt = (bt(mup)-bt(mum))./rdmu.*dmudt;
% curvature
Kr = br(cost).*dbrdr+bt(cost).*dbrdt-bt(cost).^2./R;
Kt = br(cost).*dbtdr+bt(cost).*dbtdt+br(cost).*bt(cost)./R;

curvBx1 = (Kr.*coslat+Kt.*sinlat).*coslon;
curvBy1 = (Kr.*coslat+Kt.*sinlat).*sinlon;
curvBz1 = (Kr.*sinlat-Kt.*coslat);
if 0,
disp('curvature components')
[curvBx;curvBy;curvBz]
[curvBx1;curvBy1;curvBz1]
[curvBx;curvBy;curvBz]-[curvBx1;curvBy1;curvBz1]
subplot(211), plot(Z,[curvBx;curvBz]','-o'),
legend({'K_x','K_z'})
subplot(212), plot(Z,[curvBx1;curvBz1]','-o'),
legend({'K_x','K_z'})
pause
end
end

end

end
