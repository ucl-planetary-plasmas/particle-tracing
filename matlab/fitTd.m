function fit = fitTd(t,x,tb,a4td)
% function fit = fitTd(t,x,tb,a4td)
%
% t, x : time series of longitudes, time t is in time unit,
%        latitude x is in angle unit
%
% tb   : estimate for bounce period \tau_b [time unit]
%
% a4td : angle [unit angle] for one drift period
%         * 360 for longitude in degrees
%         * 2\pi for longitude in radians
%
% fit: structure returned containing the following fields
%      kvg: fit status
%      r2 : coefficient of multiple determination
%      tb : bounce period [time unit]
%      dtb: standard deviation estimate for tb [time unit]
%      td : drift period [time unit]
%      dtd: standard deviation estimate for td [time unit]
%      f  : time series of fitted longitudes
%      nlm: class returned by fitnlm 

%
% $Id: fitTd.m,v 1.3 2019/06/14 14:09:32 patrick Exp $
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
td = 1/pol(1);

if 0,
clf
plot(t,x-polyval(pol,t)), pause
end

% initial parameters 
% [modulation amplitude, modulation period, phase, drift period]
% note that modulation period is equal to half the bounce period
p0 = [max(abs((x-polyval(pol,t)))); tb/2  ;   0; td];
if td>0,
lb = [0                           ; tb/2/3; -pi; td/3];  % lower boundary
ub = [2*p0(1)                     ; 3*tb/2;  pi; 3*td];  % upper boundary
else
lb = [0                           ; tb/2/3; -pi; 3*td];  % lower boundary
ub = [2*p0(2)                     ; 3*tb/2;  pi; td/3];  % upper boundary
end
np = 4;                     % number of parameters

% model function
model = @(p,t) t/p(4) + p(1)*sin(2*pi/p(2)*t+p(3));

% objective function for global optimisation
objfun = @(p) norm(x-model(p,t));

% various inlined print functions
pp = @(text,b)fprintf(1,'%s: %+8f,%+8f,%+8f,%+8f\n',text,b);
ps = @(text,b)fprintf(1,'%s: %9f,%9f,%9f,%9f\n',text,b);
pgof = @(text,r2,r2adj,mse)fprintf(1,'%s: %9f %9f %9f\n',text,r2,r2adj,mse);

if 0,
% Genetic algorithm global optimisation
opts = optimoptions('ga','Display','off',...
                    'PopulationSize',300,...
                    'UseParallel',false);
[p1,fval,exitflag,output,population,scores] = ...
    ga(objfun,np,[],[],[],[],lb,ub,[],opts);
%pp('p1 from ga             ',p1)
end

% pattern search global optimimsation 
opts = optimoptions('patternsearch','Display','off',...
                    'MaxIterations',600,'MeshTolerance',1e-6,...
                    'UseParallel',false);
%pp('p0                     ',p0)
[p1,fval,exitflag,output] = ...
   patternsearch(objfun,p0,[],[],[],[],lb,ub,[],opts);
%pp('p1 from patternearch   ',p1)

% nonlinear fit optimisation
display = 'iter';
display = 'final';
display = 'off';
stol = 1e-12;
robust = 'off';

prec = {'TolFun',stol,'TolX',stol,'MaxIter',400};
deriv = {'Jacobian','on','GradObj','on'};
opts = statset(prec{:},deriv{:},'Robust',robust,'UseParallel',true);

%fprintf(1,'Robust %s\n', robust);
tbl = table(t(:),x(:));

% For information on the class returned by fitnlm see
% https://uk.mathworks.com/help/stats/nonlinearmodel-class.html
nlm = fitnlm(tbl,model,p1,'Options',opts);

if 0,
pp('fitnlm                 ',nlm.Coefficients.Estimate);
ps('fitnlm std             ',nlm.Coefficients.SE);
pgof('fitnlm r2 r2adj mse    ',nlm.Rsquared.Ordinary,nlm.Rsquared.Adjusted,nlm.MSE);
end

r2 = nlm.Rsquared.Adjusted;

% fitted parameters
p = nlm.Coefficients.Estimate;
% standard deviation for fitted parameters
psd = nlm.Coefficients.SE;
f = nlm.Fitted;

% omfit=2*omb thus tb=2*tfit
tb = 2*p(2);
dtb = 2*psd(2);

td = a4td*p(4);
dtd = a4td*psd(4);

fprintf(1,'r2=%.2f tb=%6.2f s (%.2g) td=%6.2f s (%.2g)\n',...
        r2, tb, dtb, td,dtd);

fit = struct('r2',r2,'tb',tb,'dtb',dtb,'td',td,'dtd',dtd,'f',f,'nlm',nlm);

if 0,
clf
plot(t,x,t,F(t,p0),t,F(t,p)), pause
end

