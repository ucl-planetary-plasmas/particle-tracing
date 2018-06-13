function [td,dtd]=getdriftperiod(t,x,omb,a4td)
% function [td,dtd]=getdriftperiod(t,x,omb,a4td)
%
% t,x: time series of longitudes
% omb: estimate for bounce frequency \omega_b = 2\pi/\tau_b in rad/s
% a4td: angle for 1 drift period (360 for longitude in degrees, 2*pi 
%       fpr longitude in radians)



%
% $Id: getperiod.m,v 1.3 2017/01/12 15:40:52 patrick Exp $
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

% fit straight line
pol = polyfit(t,x,1);
% slope is estimate for drift angular frequency
omd = pol(1);

if 0,
clf
plot(t,x-polyval(pol,t)), pause
end

% initial parameters 
% [modulation amplitude, modulation period, phase, drift frequency]
pin = [max(abs((x-polyval(pol,t))/omd));2*omb;0;omd];

% what to fit
dp = [1; 1; 1; 1];
% fixed omb
%dp = [1; 0; 1; 1];


stol = 1e-6; 
niter = 200;
minstep = [0; 0; 0; 0];
maxstep = [Inf; Inf; Inf; Inf];
options = [minstep, maxstep];

t = t(:);
x = x(:);
wt = ones(size(t));


global verbose;
verbose = [0 0];

[f, p, kvg, iter, corp, covp, covr, stdresid, Z, r2] = ...
    leasqr (t, x, pin, @F, stol, niter, wt, dp, @dFdp, options);
% standard deviation for estimated parameters
psd = zeros(size(p));
psd(dp==1) = sqrt(diag(covp));

% 2*omb this 2*tb
tb = 2*(2*pi/p(2));
dtb = 2*(2*pi*psd(2)/p(2)^2);

td = a4td/p(4);
dtd = a4td*psd(4)/p(4)^2;

fprintf(1,'kvg=%d iter=%2d tb=%.4g (%.2g) A=%.2g (%.2g) td=%.4g (%.2g)\n',...
        kvg, iter, tb, dtb, p(1), psd(1),td,dtd);

if 0,
clf
plot(t,x,t,F(t,pin),t,F(t,p)), pause
end

function y = F(t,p)
% fprintf(1, 'called F(t,[%e,%e,%e]\n', p(1),p(2),p(3))
y = p(4)*(t + p(1)*sin(p(2)*t+p(3)));


function y = dFdp(t,f,p,dp,func)
y = [p(4)*sin(p(2)*t+p(3)),...
     p(4)*p(1)*t.*cos(p(2)*t+p(3)),...
     p(4)*p(1)*cos(p(2)*t+p(3)),...
		 t+p(1)*sin(p(2)*t+p(3))];

