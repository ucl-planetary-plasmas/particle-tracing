function [tb,dtb] = getbounceperiod(t,x,omb)
% function [tb,dtb] = getbounceperiod(t,x,omb)
%
% t,x: time series of latitudes or other periodic time series with bounce
%      period
% omb: estimate for bounce frequency omega_b = 2\pi/\tau_b

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

% initial parameters [amplitude, angular frequency, phase]
pin = [max(abs(x));omb;0];

% what to fit
dp = [1; 1; 1];
%dp = [1; 1; 0];

stol = 1e-6; 
niter = 200;
minstep = [0; 0; 0];
maxstep = [Inf; Inf; Inf];
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

tb = 2*pi/p(2);
dtb = 2*pi*psd(2)/p(2)^2;

fprintf(1,'kvg=%d iter=%2d tb=%.4g (%.2g) A=%.2g (%.2g)\n',...
        kvg, iter, tb, dtb, p(1), psd(1));

if 0,
clf
plot(t,x,t,F(t,pin),t,F(t,p)), pause
end

function y = F(t,p)
% fprintf(1, 'called F(t,[%e,%e,%e]\n', p(1),p(2),p(3))
y = p(1)*sin(p(2)*t+p(3));

function y = dFdp(t,f,p,dp,func)
y = [sin(p(2)*t+p(3)),...
     p(1)*p(2)*t.*cos(p(2)*t+p(3)),...
     p(1)*cos(p(2)*t+p(3))];

