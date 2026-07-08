function example_trajectory
% function example_trajectory

%
% $Id: example_trajectory.m,v 1.2 2026/07/08 18:52:48 patrick Exp $
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
partype = 'S+';                           % proton
Ep = 25;                                  % energy [MeV]
Ri = 70;                                  % initial equatorial distance [Rp]
ai = 10;                                  % initial equatorial pitch angle [deg]
timespec = [0,0,20,0];                    % run for 20 dipole bounce periods
savefile = [];                            % name for result mat-file
pauseOn = true;                           % flag to pause on/off.
plotOn = true;                            % flag to plot on/off.
npertc = 20;                              % number of iterations/equatorial gyroperiod


mdbtracer(mdfile,partype,Ep,Ri,ai,timespec,savefile,pauseOn,plotOn,npertc);

