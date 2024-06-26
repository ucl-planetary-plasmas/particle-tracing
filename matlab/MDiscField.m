function [Br, Bt, varargout] = MDiscField(r, theta)
% function [Br, Bt, gradB] = MDiscField(r, theta)
%
% r = radial distance in RJ
% theta = colatitude in radians
%
% Br = radial field in Tesla
% Bt = meridional field in Tesla
% gradB = gradient (r,theta) of the amplitude of the magnetic field in
%         Tesla/m
%
% The disc model filename is set with the global variable MDFile, for
% instance
%
%   global MDFile
%   MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
% 
% The disc file can be cleared by clearing the MDiscField function
%
%   clear MDiscField
% 

% $Id: MDiscField.m,v 1.10 2020/11/17 13:36:45 patrick Exp $
%
% Copyright (c) 2016-2017 Patrick Guio <patrick.guio@gmail.com>
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
persistent MD Fr Ft Fm

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

  else, % dipole field

    fprintf(1,'dipole field ... '); 
    Br = MD.v2d.Brdip'*MD.scales.bfield;
    Fr = griddedInterpolant(R,Mu,Br,opts{:});

    Bt = MD.v2d.Bthdip'*MD.scales.bfield;
    Ft = griddedInterpolant(R,Mu,Bt,opts{:});

    B = sqrt(MD.v2d.Brdip.^2+MD.v2d.Bthdip.^2)'*MD.scales.bfield;
    Fm = griddedInterpolant(R,Mu,B,opts{:});

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
rs = size(r);

% radial component
dr = 2*r(:)*EPS * MD.scales.length;
rp = r(:)*(1+EPS);
rm = r(:)*(1-EPS);
% one argument seems faster than two arguments
Xp = [rp, mu(:)];
Xm = [rm, mu(:)];
dBdr = reshape((Fm(Xp)-Fm(Xm))./dr,rs);
%dBdr = reshape((Fm(rp,mu(:))-Fm(rm,mu(:)))./dr,rs);

% tangential component
dmudt = -sqrt(1-mu(:).^2);
% add EPS to deal with mu=0
mup = mu(:)*(1+EPS)+EPS;
mum = mu(:)*(1-EPS)-EPS;
rdmu = r(:).*(mup-mum) * MD.scales.length;
% one argument seems faster than two arguments
Xp = [r(:), mup];
Xm = [r(:), mum];
dBdt = reshape((Fm(Xp)-Fm(Xm))./rdmu.*dmudt,rs);
%dBdt = reshape((Fm(r(:),mup)-Fm(r(:),mum))./rdmu.*dmudt,rs);

varargout{1} = {dBdr, dBdt};

end

return

end
