function fit = fitTb(t,x,tb)
% function fit = fitTb(t,x,tb)
%
% t, x : time series of latitudes, time t is in time unit, 
%        latitude x is in angle unit
%
% tb   : estimate for bounce period \tau_b [time unit]
%
% fit  : returned stucture containing the following fields
%        r2   : coefficient of multiple determination
%        tb   : bounce period [time unit]
%        dtb  : std estimate for tb [time unit]
%        mplat: mirror point latitude,  max(abs(x) [angle unit]
%        f    : time series of fitted latitudes
%        nlm  : class returned by fitnlm

%
% $Id: fitTb.m,v 1.2 2019/06/12 17:40:57 patrick Exp $
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

% initial parameters [amplitude, bounce period, phase offset]
p0 = [max(abs(x)); tb  ;   0];
lb = [0          ; tb/3; -pi]; % lower boundary
ub = [2*p0(1)    ; 3*tb;  pi]; % upper boundary
np = 3;                        % number of parameters

% model function
model = @(p,t) p(1)*sin(2*pi/p(2)*t+p(3));
% alternate model triangle wave
%model = @(p,t) p(1)*sawtooth(2*pi/p(2)*t+p(3)+pi/2,.5);

% objective function for global optimisation
objfun = @(p) norm(x-model(p,t));

% various inlined print functions
pp = @(text,b)fprintf(1,'%s: %+8f,%+8f,%+8f\n',text,b);
ps = @(text,b)fprintf(1,'%s: %9f,%9f,%9f\n',text,b);
pgof = @(text,r2,r2adj,mse)fprintf(1,'%s: %9f %9f %9f\n',text,r2,r2adj,mse);

if 0,
% Genetic algorithm global optimisation
opts = optimoptions('ga','Display','off',...
                    'PopulationSize',300,...
                    'UseParallel',false);
[p1,fval,exitflag,output,population,scores] = ...
    ga(objfun,np,[],[],[],[],lb,ub,[],opts);
pp('p1 from ga             ',p1)
end

% pattern search global optimimsation 
opts = optimoptions('patternsearch','Display','off',...
                    'MaxIterations',600,'MeshTolerance',1e-6,...
                    'UseParallel',false);
%pp('p0                     ',p0)
[p1,fval,exitflag,output] = ...
   patternsearch(objfun,p0,[],[],[],[],[],[],[],opts);
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

% for information on the class nlm see
% https://uk.mathworks.com/help/stats/nonlinearmodel-class.html
nlm = fitnlm(tbl,model,p1,'Options',opts);

if 0,
pp('fitnlm                 ',nlm.Coefficients.Estimate);
ps('fitnlm std             ',nlm.Coefficients.SE);
pgof('fitnlm r2 r2adj mse    ',nlm.Rsquared.Ordinary,nlm.Rsquared.Adjusted,nlm.MSE);
end

iter = 0;
% 
r2 = nlm.Rsquared.Adjusted;
% fitted parameters
p = nlm.Coefficients.Estimate;
% standard deviation for fitted parameters
psd = nlm.Coefficients.SE;
f = nlm.Fitted;

tb = p(2);
dtb = psd(2);
mplat = max(abs(x));

fprintf(1,'r2=%.2f tb=%6.2f s (%.2g) lm=%.2f deg A=%.2f deg (%.2g)\n',...
        r2,  tb, dtb, mplat, p(1), psd(1));

fit = struct('r2',r2,'tb',tb,'dtb',dtb,'mplat',mplat,'f',f,'nlm',nlm);


if 0,
clf
plot(t,x,t,F(t,p0),t,F(t,p)), pause
end

