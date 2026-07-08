function check_length_scales(plot_id, mdfile, dipole_on)
% function check_length_scales(plot_id, mdfile, dipole_on)
%
% Magnetic field length scale parameters along dipole L-shells
%
% plot_id  : id number for set of parameters to plot 
%            0: magnetic field gradient length scale
%            1: magnetic curvature 
%            2: derivative terms for directional derivative in curvatures
%
% mdfile   : magnetodisc mat-file filename 
%            default file is jup_mdisc_kh3e7_rmp90.mat
%
% dipole_on: flag for dipole or magnetodisc field
%            default value is true

%
% $Id: check_length_scales.m,v 1.1 2026/07/08 21:33:12 patrick Exp $
%
% Copyright (c) 2026 Patrick Guio <patrick.guio@gmail.com>
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

if ~exist('mdfile','var') | isempty(mdfile), % define MDfile if not given
  MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
else
  MDFile = mdfile;
end

if ~exist('dipole_on','var') | isempty(dipole_on), % define DipoleOn if not given
  DipoleOn = true;
else
  DipoleOn = dipole_on;
end

D = @(t) 1+3*cos(t).^2;

% dipole magnetic field length scale
Bm = @(r,t) sqrt(D(t))./r.^3;
dBmdr = @(r,t) -3./r.*Bm(r,t);
dBmdt = @(r,t) -3./r.*Bm(r,t).*sin(t).*cos(t)./D(t);
gradBm = @(r,t) sqrt(dBmdr(r,t).^2+dBmdt(r,t).^2);

LB = @(r,t) r/3.*D(t)./sqrt(1+7*cos(t).^2+10*cos(t).^4);
Lpar = @(r,t) r/3.*D(t).^1.5./abs(cos(t))./(3+cos(t).^2);

% dipole curvature
br = @(r,t) 2*cos(t)/D(t).^.5;
dbrdr = @(r,t) zeros(size(t));
dbrdt = @(r,t) -2*sin(t)./r./D(t).^1.5;

bt = @(r,t) sin(t)/D(t).^.5;
dbtdr = @(r,t) zeros(size(t));
dbtdt = @(r,t) 4*cos(t)./r./D(t).^1.5;

Kr = @(r,t) -3./r.*sin(t).^2.*(1+cos(t).^2)./D(t).^2;
Kr = @(r,t) -6./r.*sin(t).^2.*(3+cos(2*t))./(5+3*cos(2*t)).^2;
Kt = @(r,t) 3./r*sin(t).*cos(t).*(3+cos(t).^2)./D(t).^2;
Kt = @(r,t) 6./r.*sin(2*t).*(3+cos(2*t))./(5+3*cos(2*t)).^2;
K = @(r,t) 3./r.*sin(t)./D(t).^1.5;
K = @(r,t) 3*sqrt(2)./r.*sin(t).*(3+cos(2*t))./(5+3*cos(2*t)).^1.5;

Rc = @(r,t) r/3.*D(t).^1.5./sin(t);
Rc = @(r,t) r/3/sqrt(2).*(5+3*cos(2*t)).^1.5./sin(t)./(3+cos(2*t));


clear MDiscField
load(MDFile,'MD')

Bp = MD.scales.bfield;
Rp = MD.scales.length;

L = 5:5:30;
% dipole mirror point at Rp provides loss cone
lm = acos(1./sqrt(L));

cols = get(gca,'ColorOrder');
c1 = cols(1,:);
c2 = cols(2,:);
c3 = cols(3,:);
c4 = cols(4,:);

clf

for i=1:length(L),
  t{i} = linspace(pi/2,pi/2+lm(i),100);
	% dipole field line
	r{i} = L(i)*sin(t{i}).^2;
	if 0,
	x = (t{i}-pi/2)*180/pi; % colatitude
	xlbl = 'colatitude';
	xLim = [0,90];
	else
  x = r{i}; % radial distance
	xlbl = 'radial distance';
	xLim = [0,max(L)];
	end
	[Br, Bt, gradB, curvB] = MDiscField(r{i}, t{i});
	Rcs{i}=Rc(r{i},t{i});
	LBs{i}=LB(r{i},t{i});
	%subplot(length(L),1,i)
	switch plot_id,
	case 0,
	  yn = {gradB{1}*Rp/Bp, gradB{2}*Rp/Bp,... 
          gradB{3}*Rp/Bp, gradB{4}/Rp};
    ya = {dBmdr(r{i},t{i}), dBmdt(r{i},t{i}),... 
          gradBm(r{i},t{i}), LB(r{i},t{i})};
    ylbl = {'\partial{B}/\partial{r}', '\partial{B}/\partial\theta',...
            '|\nabla{B}|','L_B'};
    myPlot(x,yn,ya,xLim,xlbl,ylbl,i==1,DipoleOn)

	case 1,
	  yn = {curvB{1}*Rp, curvB{2}*Rp,... 
          curvB{3}*Rp, curvB{4}/Rp};
    ya = {Kr(r{i},t{i}), Kt(r{i},t{i}),... 
          K(r{i},t{i}), Rc(r{i},t{i})};
    ylbl = {'\kappa_r', '\kappa_t',...
            '\kappa','R_c'};
    myPlot(x,yn,ya,xLim,xlbl,ylbl,i==1,DipoleOn)

  case 2,
	  yn = {curvB{5}*Rp, curvB{6}*Rp,... 
		      curvB{7}*Rp, curvB{8}*Rp};
		ya = {dbrdr(r{i},t{i}), dbrdt(r{i},t{i}),... 
		      dbtdr(r{i},t{i}), dbtdt(r{i},t{i})};
		ylbl = {'\partial{b_r}/\partial{r}', '1/r\partial{b_r}/\partial\theta',...
		        '\partial{b_\theta}/\partial{r}','1/r\partial{b_\theta}/\partial\theta'};
		myPlot(x,yn,ya,xLim,xlbl,ylbl,i==1,DipoleOn)
	end
end


function myPlot(x,yn,ya,xLim,xlbl,ylbl,lblOn,dipoleOn)

cols = get(gca,'ColorOrder');
c1 = cols(1,:);
c2 = cols(2,:);
c3 = cols(3,:);
c4 = cols(4,:);

np = length(yn);

for i=1:np,

  subplot(np,1,i),
  hold on
  plot(x,yn{i},'Color', c1)
  plot(x,ya{i},'Color', c2)
  if lblOn, 
    ylabel(ylbl{i}), 
    if i==np, 
      xlabel(xlbl), 
		end
  end
  xlim(xLim), 
  hold off

end


