function example_length_scale
% function example_length_scale

%
% $Id: example_length_scales.m,v 1.2 2026/07/08 21:28:23 patrick Exp $
%
% Copyright (c) 2025 Patrick Guio <patrick.guio@gmail.com>
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

mdfile = 'jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
partype = 'p';                    % proton
Ep = .25;                         % energy [MeV]
Ri = 30;                          % initial equatorial distance [Rp]
ai = 30;                          % initial equatorial pitch angle [deg]
timespec = [0,0,2,0];             % run for 2 dipole bounce periods
savefile = [];                    % name for result mat-file
pauseOn = false;                  % flag to pause on/off.
plotOn = false;                   % flag to plot on/off.
npertc = 20;                      % number of iterations/equatorial gyroperiod

savefile='jup_scale';

mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc);

track = load(savefile);

clf

Re = track.Re;

% simulation 1 (no centrifugal force)
t1 = track.s1.t/track.s1.latfit.tb;
rho1 = track.s1.rhog;
Lb1 = track.s1.Lb;
Eb1 = track.s1.Eb;
Rc1 = track.s1.Rc;
Ec1 = track.s1.Ec;

% simulation 2 (with centrifugal force)
t2 = track.s2.t/track.s2.latfit.tb;
rho2 = track.s2.rhog;
Lb2 = track.s2.Lb;
Eb2 = track.s2.Eb;
Rc2 = track.s2.Rc;
Ec2 = track.s2.Ec;

subplot(411)
semilogy(t1,rho1/Re,t1,Lb1/Re,t1,Rc1/Re)
title('no centrifugal force');
legend({'\rho','L_B','R_c'})

subplot(412)
plot(t1,Eb1,t1,Ec1)
legend({'E_b=\rho/L_B','E_c=\rho/R_c'})

subplot(413)
semilogy(t2,rho2/Re,t2,Lb2/Re,t2,Rc2/Re)
title('with centrifugal force');

subplot(414)
plot(t2,Eb2,t2,Ec2)
xlabel('time t/t_b')

return

x1 = (rho1-rho2)./rho1*100;
x2 = (Lb1-Lb2)./Lb1*100;
x3 = (Rc1-Rc2)./Rc1*100;
min(x1(:)),max(x1(:))
plot(t1,x1,t1,x2,t1,x3)
