function [Br, Bt, varargout] = MDiscField(r, theta)
% function [Br, Bt, gradB, curvB] = MDiscField(r, theta)
%
% r      : radial distance [Rp]
% theta  : colatitude [rad]
%
% Br     : radial magnetic field component  [T]
% Bt     : meridional magnetic field component [T]
%
% gradB = {dBm/dr, 1/r dBm/dt, |grad Bm|, LB=Bm/|grad Bm|}
%        : optional magnetic field amplitude gradient parameters 
%          in polar coordinates
%         {[T/m], [T/m], [T/m], [m]}
%
% curvB = {Kr, Kt, Km, Rc=1/Km, dbr/dr, 1/r dbr/dt, dbt/dr, 1/r dbt/dt} 
%        : optional magnetic field curvature parameters 
%          in polar coordinates
%         {[m-1], [m-1], [m-1], [m], [T/m], [T/m], [T/m], [T/m]}
%
% The magnetodisc model filename is set with the global variable MDFile.
% For instance
%
%   global MDFile
%   MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
% 
% The magnetodisc data can be removed by clearing the MDiscField function
%
%   clear MDiscField
% 
% By default the magnetodisc data is selected, to select the dipole data set
% the global variable DipoleOn to true
% For instance
% 
%   global DipoleOn
%   DipoleOn = true;

% $Id: MDiscField.m,v 1.15 2026/07/08 13:03:24 patrick Exp $
%
% Copyright (c) 2016 Patrick Guio <patrick.guio@gmail.com>
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

global MDFile DipoleOn
persistent MD Fr Ft Fm fr ft

if isempty(MD),
  fprintf(1,'Loading %s ... ', MDFile);
  eval(['load ', MDFile]);
  fprintf(1,'done.\n');

  %opts = {'linear', 'linear'};
  %opts = {'spline', 'spline'};
  opts = {'cubic', 'cubic'}; % for regular grid only
  %opts = {'makima', 'makima'};

  R = MD.c2d.r'; 
  Mu = MD.c2d.mu';

  fprintf(1,'Generating field functions (%s,%s) from ',opts{:});
  if isempty(DipoleOn) || ~DipoleOn, % magnetodisk field

    fprintf(1,'magnetodisk field ... '); 
    Br = MD.v2d.Br'*MD.scales.bfield;
    Fr = griddedInterpolant(R,Mu,Br,opts{:});

    Bt = MD.v2d.Bth'*MD.scales.bfield;
    Ft = griddedInterpolant(R,Mu,Bt,opts{:});

    B = sqrt(MD.v2d.Br.^2+MD.v2d.Bth.^2)'*MD.scales.bfield;
    Fm = griddedInterpolant(R,Mu,B,opts{:});

		br = Br./B;
    fr = griddedInterpolant(R,Mu,br,opts{:});

		bt = Bt./B;
    ft = griddedInterpolant(R,Mu,bt,opts{:});

  else, % dipole field

    fprintf(1,'dipole field ... '); 
    Br = MD.v2d.Brdip'*MD.scales.bfield;
    Fr = griddedInterpolant(R,Mu,Br,opts{:});

    Bt = MD.v2d.Bthdip'*MD.scales.bfield;
    Ft = griddedInterpolant(R,Mu,Bt,opts{:});

    B = sqrt(MD.v2d.Brdip.^2+MD.v2d.Bthdip.^2)'*MD.scales.bfield;
    Fm = griddedInterpolant(R,Mu,B,opts{:});

		br = Br./B;
    fr = griddedInterpolant(R,Mu,br,opts{:});

		bt = Bt./B;
    ft = griddedInterpolant(R,Mu,bt,opts{:});

	end
  fprintf(1,'done.\n');

end

mu = cos(theta);

% Field components in model (tilted, rotating) frame

% fast gridded data interpolant
rs = size(r);
% one argument X seems faster than two arguments r(:), mu(:)
X = [r(:), mu(:)];
Br = reshape(Fr(X), rs);
Bt = reshape(Ft(X), rs);
%Br = reshape(Fr(r(:), mu(:)), rs);
%Bt = reshape(Ft(r(:), mu(:)), rs);

if nargout>2,

EPS = 1e-6;

% radial component d/dr
dr = 2*r(:)*EPS * MD.scales.length; % [m]
rp = r(:)*(1+EPS); % normalised
rm = r(:)*(1-EPS); % normalised
% one argument seems faster than two arguments
Xp = [rp, mu(:)];
Xm = [rm, mu(:)];
dBdr = reshape((Fm(Xp)-Fm(Xm))./dr,rs); % [T/m]
%dBdr = reshape((Fm(rp,mu(:))-Fm(rm,mu(:)))./dr,rs);

if nargout>2, % B for gradB 

B = reshape(Fm(X), rs);

if nargout>3, % radial components for curvB

br = reshape(fr(X), rs);
bt = reshape(ft(X), rs);

dbrdr = reshape((fr(Xp)-fr(Xm))./dr,rs); % [m-1]
dbtdr = reshape((ft(Xp)-ft(Xm))./dr,rs); % [m-1]

end
end

% tangential component 1/r d/dt
dmudt = -sqrt(1-mu(:).^2);
% add EPS to deal with mu=0
mup = mu(:)*(1+EPS)+EPS;
mum = mu(:)*(1-EPS)-EPS;
rdmu = r(:).*(mup-mum) * MD.scales.length; % [m]
% one argument seems faster than two arguments
Xp = [r(:), mup];
Xm = [r(:), mum];
dBdt = reshape((Fm(Xp)-Fm(Xm))./rdmu.*dmudt,rs); % [T/m]
%dBdt = reshape((Fm(r(:),mup)-Fm(r(:),mum))./rdmu.*dmudt,rs);

gradBm = sqrt(dBdr.^2+dBdt.^2);
LB = B./gradBm;

gradB = {dBdr, dBdt, gradBm, LB};

varargout{1} = gradB;

if nargout>3, % tangential components for curvB

dbrdt = reshape((fr(Xp)-fr(Xm))./rdmu.*dmudt,rs); % [m-1]
dbtdt = reshape((ft(Xp)-ft(Xm))./rdmu.*dmudt,rs); % [m-1]

Kr = br.*dbrdr + bt.*dbrdt - bt.^2./(r*MD.scales.length);
Kt = br.*dbtdr + bt.*dbtdt + br.*bt./(r*MD.scales.length);

Km = sqrt(Kr.^2+Kt.^2);

Rc = 1./Km;

curvB = {Kr, Kt, Km, Rc, dbrdr, dbrdt, dbtdr, dbtdt};

varargout{2} = curvB;

end

end

return

end
