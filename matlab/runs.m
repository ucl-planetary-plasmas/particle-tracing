%
% $Id: runs.m,v 1.7 2026/07/08 21:46:39 patrick Exp $
%
% Copyright (c) 2017 Patrick Guio <patrick.guio@gmail.com>
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

% Earth p+ energy 10 MeV 
% initial condition Req=1.5 R_E, pitch angle=80 deg
% simulation time 3 bounce period
%trajectory_main(1.5,0,0,1000,80,'p','E',0,'N',-1,1000,1,'A',3);
planet = 'earth';                 % Earth
partype = 'p';                    % proton
Ep = 10;                          % energy [MeV]
Ri = 1.5;                         % initial equatorial distance [Rp]
ai = 80;                          % initial equatorial pitch angle [deg]
timespec = [0,0,3,0];             % run for 3 dipole bounce periods
dipbtracer(planet,partype,Ep,Ri,ai,timespec);

% Dipolar Jupiter p+ energy 100 MeV 
% initial condition Req=1.5 R_J, pitch angle 80=deg
% simulation time 2 bounce period
%trajectory_main(1.5,0,0,10000,80,'p','J',0,'N',-1,1000,1,'A',2);
planet = 'jupiter';               % Jupiter
partype = 'p';                    % proton
Ep = 100;                         % energy [MeV]
Ri = 1.5;                         % initial equatorial distance [Rp]
ai = 80;                          % initial equatorial pitch angle [deg]
timespec = [0,0,2,0];             % run for 2 dipole bounce periods
dipbtracer(planet,partype,Ep,Ri,ai,timespec);

% Magnetodisc Jupiter p+ energy 100 MeV
% initial condition Req=1.5 R_J, pitch angle=80 deg
% simulation time 2 bounce period
%trajectory_main(1.5,0,0,10000,80,'p','J',0,'Y',-1,1000,1,'A',2);
mdfile = 'jup_mdisc_kh3e7_rmp90'; % magnetodisc mat-file
partype = 'p';                    % proton
Ep = 100;                         % energy [MeV]
Ri = 1.6;                         % initial equatorial distance [Rp]
ai = 80;                          % initial equatorial pitch angle [deg]
timespec = [0,0,2,0];             % run for 2 dipole bounce periods
mdbtracer(mdfile,partype,Ep,Ri,ai,timespec);


% Compare algorithm to integrate equation of motion for Earth
% simulation time 2 bounce periods
cmptraj_earth([0,0,2,0])

% Compare algorithm to integrate equation of motion for Jupiter
% simulation time 1 bounce period
cmptraj_jupiter([0,0,1,0])

% Compare algorithm to integrate equation of motion for Saturn
% simulation time 2 bounce periods
cmptraj_saturn([0,0,2,0])

