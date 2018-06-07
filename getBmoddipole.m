function Bmod = getBmoddipole(m, rm, r, RP)

%    m  : magnetic moment vector
%    rm : position (origin) of the magnetic moment
%    r  : cartesion coordinates where to calculate magnetic field
%    B  : dipole magnetic field (cell array of length(m), each cell is
%         an array of size(r))
%    RP : PLanetary radius / Earth's
%

% $Id: getBmoddipole.m,v 1.1 2017/10/06 11:15:22 patrick Exp $
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
% Subtract moment position (which is likely to be at 0,0,0 anyway)
R = r - rm;

% Bx, By and Bz
Bx = 10000./norm(R).^3*(3*dot(m,R).*R(1)./norm(R).^2 - m(1))/RP.^3;
By = 10000./norm(R).^3*(3*dot(m,R).*R(2)./norm(R).^2 - m(2))/RP.^3;
Bz = 10000./norm(R).^3*(3*dot(m,R).*R(3)./norm(R).^2 - m(3))/RP.^3;

Bmod = sqrt(Bx^2 + By^2 + Bz^2);

