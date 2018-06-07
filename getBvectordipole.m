function B = getBvectordipole(m, rm, r, RP)
% function B = dipoleMagneticField3D(m, rm, r)
%
%    m  : magnetic moment vector
%    rm : position of the magnetic moment
%    r  : cartesion coordinates where to calculate magnetic field
%    B  : dipole magnetic field (cell array of length(m), each cell is
%         an array of size(r))
%

% $Id: getBvectordipole.m,v 1.1 2017/10/06 11:15:22 patrick Exp $
%
% Copyright (c) 2009-2017 Patrick Guio <patrick.guio@gmail.com>
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

% Subtract moment position (which is likely to be at 0,0,0 anyway)
X = X-rm(1);
Y = Y-rm(2);
Z = Z-rm(3);

R2 = X.^2 + Y.^2 + Z.^2; % Modulus square
R3 = R2.^(3/2); % Cube

mDotR = m(1)*X + m(2)*Y + m(3)*Z; % scalar product M.R

%EarthEQ = 0.306;
% Bx, By and Bz
B{1} = 10000./R3.*(3*mDotR.*X./R2 - m(1))/RP.^3;
B{2} = 10000./R3.*(3*mDotR.*Y./R2 - m(2))/RP.^3;
B{3} = 10000./R3.*(3*mDotR.*Z./R2 - m(3))/RP.^3;


