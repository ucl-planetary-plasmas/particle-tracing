function trajectory_main(R0, mlat0, mlong0, KE0, pitcheq, particle, planet, llimit, mdisc, initial_angle, outputthinning, solver_choice, Method, bounces)
%Input parameters for particle initially in Equatorial plane
%R0         Initial distance from planet centre / planet's diameter.
%mlat0      Initial magnetic latitude.
%mlong0     Initial magnetic longitude.
%KE0        Initial kinetic energy of particle (gamma-1)*mc2
%pitcheq    Initial pitch angle in degrees.
%particle  Type 'p' for proton. Anything else is electron.
%planet   Type 'j' or 'J' for Jupiter, s or 'S' for Saturn. Anything else
%           for Earth.
%llimit     Enter maximum order of Legendre polynomial expansion for
%           calculating field potential. Enter 0 to use pure dipole
%           equation.
%mdisc
%initial_angle  This defines where in the bounce cycle the trajectory
%           starts.  If set to 0 or 180 degrees, the initial veloity will
%           be perpendicular to the guiding centre position vector R, which 
%           =R0(R0, mlat0, mlong0). If set to 90, the initial v will be 
%           anti-parallel to R0 and if set to 270 it will be parallel to R0.  
%           Any angle between 0 and 360 degrees may be entered. 
%           If set to -1, the angle is set so that the particle
%           starts at the same distance from the planet as is its
%           (predicted) guided centre.
%outputthininning  During each simulation, up-to-date information is
%           recorded to a text file at each integration step, which can used later for detailed
%           analysis.  This parameter limits the amount of data recorded to
%           prevent the file size from becoming overwhelming.  It is
%           recommended that this be set to 100 for shorter simulations.
%           For longer sims, set it to 1000 or 10000 even.
%solver_choice, enter either 1 or 2:  1 - ode45 (non-stiff), 2 - ode15s (stiff)
%Method     Mirror point detection method. Enter 'A' if doing < 2 bounces or you 
%           believe that the bounce period remians constant. Enter 'B' if period is inclined to drift.  
%bounces    Number of bounce periods to simulate (based on prediction of
%           bounce period).
%     
%           Typical input example: trajectory_main(1.5,0,0,1000,80,'p','E',0,'N',-1,1000,1,'A',3)

%
% $Id: trajectory_main.m,v 1.4 2017/10/18 09:40:00 patrick Exp $
%
% Copyright (c) 2015-2017 Patrick Guio <patrick.guio@gmail.com>
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
if solver_choice < 1 || solver_choice > 2
    solver_choice = 1;
end
if (upper(Method) ~= 'B')
    Method = 'A';
else Method = 'B';
end
stepspergyro = 10;
if lower(mdisc) == 'y'
    llimit = 0;
end
if llimit > 13
    llimit = 13; % llimit cannot be larger than 13.
end
Electroncharge = -1.6022E-19;
Electronmass = 9.1094E-31;
restmassenergy = 511.99891; %keV
Protonmass = Electronmass*1836;
Eartheqradius = 12736/2; %Earth equatorial radius in km
Earthaverageradius = 6371; %Earth equatorial radius in km
Earthmagmoment = 7940000;  %Mag momnet in Tkm3

magmoment = Earthmagmoment; 
planetname = 'Earth';%If 'planet' in argument list is anything other than J or S
                            %we make it the Earth.
planetnameposessive = 'Earth''s';
planeteqradius = Eartheqradius;

planet = upper(planet); %converts lowercase to uppercase.
switch (planet)
    case 'J',
        planetname = 'Jupiter';
        planetnameposessive = 'Jovian';
        magmoment = 1.564e11;
        planeteqradius = 71492;

    case 'S',
        planetname = 'Saturn';
        planetnameposessive = 'Kronian';
        magmoment = 4.2e9;%532;
        planeteqradius = 120536/2;
end

particlemass = Electronmass; %default particle is electron
particlecharge = Electroncharge;
particlename = 'Electron';
if particle == 'p' 
    particlemass = Protonmass;
    particlecharge = -Electroncharge;
    particlename = 'Proton';
    %outputthinning = 100;
    restmassenergy = 938272.0813;
end

[x0, y0, z0] = polartocart(R0, mlat0, mlong0);
initial_GCposition_cartesian = [x0, y0, z0]; %Initial position of guiding centre.
L0 = R0/(cosd(mlat0)).^2;  %inital L-value.

c = 299792.458; %speed of light in km/s
Gamma0 = KE0/restmassenergy + 1;
fprintf(1,'Initial gamma = %4g.\n',Gamma0);
speed0c = sqrt(1-1/(Gamma0^2));
speed0kps = speed0c*c;  %Initial speed of particle in km/s
speed0 = speed0kps/planeteqradius;  %  Initial speed of particle converted to planetary radii per second.
B0 = 0;
dr = 1; %distance length for measuring potential gradient in km.
xplus = [x0+dr/planeteqradius y0 z0];
xminus = [x0-dr/planeteqradius y0 z0];
yplus = [x0 y0+dr/planeteqradius z0];
yminus = [x0 y0-dr/planeteqradius z0];
zplus = [x0 y0 z0+dr/planeteqradius];
zminus = [x0 y0 z0-dr/planeteqradius];
dipole_vector = [0 0 0];
if lower(mdisc) == 'y'
    global MDFile
		switch (planet)
		  case 'J', MDFile = 'jup_mdisc_kh3e7_rmp90.mat';
			case 'S', MDFile = 'sat_mdisc_kh2e6_rmp25.mat';
      otherwise, error('no MDFile for that planet');
		end
    % clear function with persistent variables
    clear MDiscField
    [Brad, Blat] = MDiscField(R0,(90-mlat0)*pi/180);
    B0 = sqrt(Brad.^2 + Blat.^2)*10000;
elseif llimit < 1
    B0 = getBmoddipole([0, 0,magmoment],dipole_vector , initial_GCposition_cartesian , planeteqradius); %calculates modulus of B-field at initial (guiding centre) position in gauss.
      %This is only used for calculating gyroperiod prediction.
else   %Finds magnetic potenitial at points close to desired postion.
    Vxplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(xplus,mlong0),Earthaverageradius, llimit);
    Vxminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(xminus,mlong0),Earthaverageradius, llimit);
    B_x = -(Vxplus - Vxminus)/dr/2; %gradient gives B-field components
    Vyplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(yplus,mlong0),Earthaverageradius, llimit);
    Vyminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(yminus,mlong0),Earthaverageradius, llimit);
    B_y = -(Vyplus - Vyminus)/dr/2;
    Vzplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(zplus,mlong0),Earthaverageradius, llimit);
    Vzminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(zminus,mlong0),Earthaverageradius, llimit);
    B_z = -(Vzplus - Vzminus)/dr/2;
    B0 = sqrt(B_x.^2 + B_y.^2 + B_z.^2);
end   
      
fprintf(1,'Magnetic field strength at initial guiding cente position = %4g G.\n',B0);
Bounceperiodpred = L0*planeteqradius*sqrt(2)/speed0kps*(3.7-1.6*sind(pitcheq)); %Bounce period prediction for dipole field.
EquatorialGyroperiodpred = 2 * pi * particlemass*Gamma0/(abs(particlecharge)*B0/10^4); %Inital gyroperiod prediction  2 pi m / q B
                                               %Drift period prediction.
Driftorbitperiodpred = pi*abs(particlecharge)*getBmoddipole([0, 0,magmoment],dipole_vector , [1,0,0] , planeteqradius)/10000*planeteqradius.^2*1e6/3/L0/(0.5*particlemass*speed0kps.^2*1e6*Gamma0)/(0.35+0.15*sind(pitcheq)); 
Lossconeangle = asind(1/(4*L0.^6-3*L0.^5)^0.25); %Loss cone angle prediction.

EquatorialGyroradiuspred = abs(EquatorialGyroperiodpred*speed0kps*sind(pitcheq)/(2*pi)); %Initial Gyroradius in km.
EquatorialGyroradiuspredRP = EquatorialGyroradiuspred/planeteqradius; %Initial Gyroradius in planetary radii.
if initial_angle == -1
    if particlecharge >= 0
        initial_angle = 270-atand(EquatorialGyroradiuspredRP/2/R0);
    else 
        initial_angle = 90 + atand(EquatorialGyroradiuspredRP/2/R0);
    end
elseif particlecharge < 0
    initial_angle = initial_angle + 180; %Electron start postion correctionn is reversed.
end
if initial_angle > 360
   initial_angle = initial_angle - 360;
end
speed0g = speed0*EquatorialGyroperiodpred; %speed0g is inital speed of particle in planetary radii per gyroperiod.
output_filename = strcat(particlename,'_',planetname,'_L',num2str(L0,3),'_lat',num2str(mlat0,3),'_long',num2str(mlong0,3),'_',num2str(KE0,4),'keV_p',num2str(pitcheq,3),'_pa',num2str(initial_angle,3));
if lower(mdisc) == 'y'
    output_filename = strcat(output_filename,'mdisc_', num2str(bounces,3),'bounces');
elseif llimit < 1
    output_filename = strcat(output_filename,'dipole_', num2str(bounces,3),'bounces');
else 
    output_filename = strcat(output_filename,'IGRF_llimit_',num2str(llimit), num2str(bounces,3),'_bounces');
end
trajectory_filename = strcat('trajectory_',output_filename,'.txt');
gyroperiods_filename = strcat('gyroperiods_',output_filename,'.txt');
bounceperiods_filename = strcat('bounceperiods_',output_filename,'.txt');
mirrorpoints_filename = strcat('mirrorpoints_',output_filename,'.txt');
results_filename = strcat('results_log_',output_filename,'.txt');
mat_filename = strcat(output_filename,'.mat');
if exist('results log.txt','file') ~= 2  %If 'results log.txt' does not exist it has to be created.
    resultsfile = fopen(results_filename,'a');
    fprintf(resultsfile, 'Date & time\tPlanet\tParticle\tInitial L\tInitial magnetic latitude (%s)\tInitial magnetic longitude (%s)\tK.E. (keV)\tGamma\tEquatorial pitch angle (%s)\tPlanet equatorial radius (km)\tPlanet magnetic moment (Tkm3)\tEqu. Gyroradius / R0\tInitial phase angle (%s)\tMethod\tNo. of bounces\tSteps per predicted gyro-orbit\tl limit\tsolver\tAbsolute tolerance\tRelative tolerance\tB-field strength at initial GC poistion (G)\tinitial speed (km/s)\t', char(176), char(176), char(176), char(176));
    fprintf(resultsfile, 'Predicted equatorial gyro-radius at GC (km)\tEqu. Gyroradius / R0\tLoss cone angle (%s)\tGreatest magnetic latitude (%s)\tMax R (Rp)\tMin R (Rp)\tMax Rcyl (Rp)\tMin Rcyl (Rp)\tMax z (Rp)\tMin z (Rp)\tMax magnetic latitude (%s)\tMin magnetic latitude (%s)\tMax magnetic longitude (%s)\tMin magnetic longitude (%s)\tMax L\tMin L\tGreatest gyroperiod, calculated from local magnetic field (s)\t', char(176), char(176), char(176), char(176),	char(176), char(176));
    fprintf(resultsfile, 'Smallest gyroperiod, calculated from local magnetic field (s)\tGreatest gyro-radius, calculated from local magnetic field (s)\tSmallest gyro-radius,calculated from local magnetic field (s)\tSimulation time (s)\tEquatorial gyroperiod predicted (s)\tEquatorial gyroperiod from simulation (s)\tMax gyroperiod from simulation (s)\tMin gyroperiod from simulation (s)\tEquatorial gyroperiod error (s)\t');
    fprintf(resultsfile, 'Equatorial gyroperiod discrepancy (%%)\tBounce period predicted (s)\tBounce period from simulation (s)\tBounce period error (s)\tBounce period discrepancy (%%)\tBounce period from sim. average (s)\tBounce period from sim. max (s)\tBounce period from sim. min (s)\tBounce period from sim. spread (s)\tBounce period from sim. SD (s)\tDrift orbital period predicted (s)\tDrift orbital period from simulation, polyfit method (s)\t');
    fprintf(resultsfile, 'Drift orbital period error, polyfit method (s)\tDrift orbital period discrepancy, polyfit method (%%)\tDrift orbital period from simulation, simple gradient method (s)\tDrift orbital period from simulation, simple gradient method error (s)\tDrift orbital period  simple gradient method discrepancy (%%)\tDrift orbit period delta long / tmax\tMirror point latitude predicted (%s)\tMirror point latitude from simulation  (%s)\t', char(176), char(176));
    fprintf(resultsfile, 'Mirror point latitude error (%s)\tnumber of mirror points\tMirror point latitudes standard deviation (%s)\tMirror point latitudes max (%s)\tMirror point latitudes min (%s)\tMirror point latitude discrepancy  alternative (%%)\tMu mean(Tm3)\tMu SD(Tm3)\tMax magnetic moment (Tm3)\tMin magnetic moment (Tm3)\tMu spread /Mu mean\tPeriod from FFT 1 (s)\tPeriod from FFT 2 (s)\tPeriod from FFT 3 (s)\tPeriod from FFT 4 (s)\t', char(176), char(176), char(176), char(176));
    fprintf(resultsfile, 'Time resolution (s)\tNumber of time divisions\tParticle flight time (s)\tFate of particle\n');
else    resultsfile = fopen(results_filename,'a');
end
trajectoryfile = fopen(trajectory_filename,'w');
logfile = fopen('trajectories log.txt','a');
fprintf(trajectoryfile,'Time (s)\tx (Rp)\ty (Rp)\tz (Rp)\tR (Rp)\tRcyl (Rp)\tmag lat (%c)\tmag long (%c)\tL\tvx (Rp/s)\tvy (Rp/s)\tvz (Rp/s)\tmu\tpitch (%c)\tBmod (G)\tGamma\tspeed (km/s)\tgyroperiod (s)\tgyroradius (km)\t', char(176), char(176), char(176));
fprintf(trajectoryfile, '%s\t%s\t%s\tinitial L =\t%.4g\tinitial magnetic latitude =\t%4g\tinitial magnetic longitude =\t%4g\tenergy  =\t%.4g\tkeV\tequatorial pitch angle =\t%.4g\tdeg\tmagnetic moment =\t%.4g\tTkm3\tinitial phase angle =\t%4g\tdeg\t%3g\tbounce periods\t%d\tsteps per gyro.\t', datestr(now), planetname,particlename, L0, mlat0, mlong0, KE0, pitcheq, magmoment, initial_angle, bounces, stepspergyro);
fprintf(logfile, '%s: %s, %s, initial L =\t%.4g\tinitial magnetic latitude =\t%4g\t%c, initial magnetic longitude =\t%4g\t%c, initial KE  =\t%.4g\tkeV, equatorial pitch angle =\t%.4g\tdeg, magnetic moment =\t%.4g\tTkm3,\tinitial phase angle = \t%4g\t%c,\tMirror point detection method =\t%d,\t%3g bounce periods,\t%d steps per gyro.\n', datestr(now), planetname,particlename, L0, mlat0, char(176), mlong0, char(176), KE0, pitcheq, magmoment, initial_angle, char(176), Method, bounces, stepspergyro);
fprintf(resultsfile, '%s\t%s\t%s\t%.4g\t%5g\t%5g\t%.4g\t%.6g\t%4g\t%4g\t%4g\t%4g\t%4g\t%c\t%3g\t%d\t', datestr(now), planetname,particlename, L0, mlat0, mlong0, KE0, Gamma0, pitcheq, planeteqradius, magmoment, EquatorialGyroradiuspredRP/R0,initial_angle, Method, bounces, stepspergyro);
gyroperiodsfile = fopen(gyroperiods_filename,'w');
mirrorpointsfile = fopen(mirrorpoints_filename,'w');
fprintf(gyroperiodsfile, '%s\t%s\t%s\tinitial L =\t%.4g\tinitial magnetic latitude =\t%4g\tinitial magnetic longitude =\t%4g\tenergy  =\t%.4g\tkeV\tequatorial pitch angle =\t%.4g\tdeg\tmagnetic moment =\t%.4g\tTkm3\t%3g\tinitial phase angle = %4g\tdeg\tbounce periods\t%d\tsteps per gyro.\t', datestr(now), planetname,particlename, L0, mlat0, mlong0, KE0, pitcheq, magmoment, initial_angle, bounces, stepspergyro);
fprintf(mirrorpointsfile, '%s\t%s\t%s\tinitial L =\t%.4g\tinitial magnetic latitude =\t%4g\tinitial magnetic longitude =\t%4g\tenergy  =\t%.4g\tkeV\tequatorial pitch angle =\t%.4g\tdeg\tmagnetic moment =\t%.4g\tTkm3\t%3g\tinitial phase angle = %4g\tdeg\tbounce periods\t%d\tsteps per gyro.\t', datestr(now), planetname,particlename, L0, mlat0, mlong0, KE0, pitcheq, magmoment, initial_angle, bounces, stepspergyro);
if lower(mdisc) == 'y'
   fprintf(trajectoryfile, 'Magnetodisc.\t\t');
   fprintf(gyroperiodsfile, 'Magnetodisc.\n');
   fprintf(mirrorpointsfile, 'Magnetodisc.\n');
   fprintf(resultsfile, 'Magnetodisc\t');
elseif llimit < 1
   fprintf(trajectoryfile, 'Dipole.\t\t');
   fprintf(gyroperiodsfile, 'Dipole.\n');
   fprintf(mirrorpointsfile, 'Dipole.\n');
   fprintf(resultsfile, 'dipole\t');
else 
    fprintf(trajectoryfile, 'IGRF, l limit =\t%d\n', llimit);
    fprintf(gyroperiodsfile, 'IGRF, l limit =\t%d\n', llimit);
    fprintf(mirrorpointsfile, 'IGRF, l limit =\t%d\n', llimit);
    fprintf(resultsfile, '%d\t', llimit);
end
fprintf(gyroperiodsfile,'Time (s)\tGyroperiod (s)');
fprintf(mirrorpointsfile,'Time (s)\tMirror point latitude (s)\tMirror point latitude uncertanty(s)');

fprintf(1,'Initial cyclotron period prediction = %4g s.  Initial bounce period prediction = %4g s.\nInitial drift orbital period prediction = %4g s.  Loss cone angle prediction = %3g degrees.\n',EquatorialGyroperiodpred, Bounceperiodpred, Driftorbitperiodpred, Lossconeangle);
fprintf(logfile,'Initial cyclotron period prediction =\t%4g\ts. Initial bounce period prediction =\t%4g\ts.\nInitial drift orbital period prediction =\t%4g\ts.  Loss cone angle =\t%3g\tdegrees.\n',EquatorialGyroperiodpred, Bounceperiodpred, Driftorbitperiodpred, Lossconeangle);

fprintf(1,'Initial Gyroradius = %4g km or %4g %s radii.\n', EquatorialGyroradiuspred, EquatorialGyroradiuspredRP, planetnameposessive);
fprintf(1,'Initial speed = %4gc = %4g %s radii per second or %4g %s radii per gyroperiod.\n',speed0c, speed0, planetnameposessive, speed0g, planetnameposessive);
fprintf(logfile,'Initial Gyroradius =\t%4g\tkm or\t%4g\t%s radii.\n', EquatorialGyroradiuspred, EquatorialGyroradiuspredRP, planetnameposessive);
fprintf(logfile,'Initial speed =\t%4g\tc =\t%4g\t%s radii per second or\t%4g\t%s radii per gyroperiod.\n',speed0c, speed0, planetnameposessive, speed0g, planetnameposessive);

fprintf(trajectoryfile, 'Initial cyclotron period prediction =\t%4g\ts\tInitial gyroradius prediction =\t%4g\tkm\tInitial bounce period prediction =\t%4g\ts\tInitial drift period prediction =\t%4g\ts\tLoss cone angle =\t%3g\tdegrees,\t', EquatorialGyroperiodpred, EquatorialGyroradiuspred, Bounceperiodpred, Driftorbitperiodpred,Lossconeangle);
M = sym('M','positive');  %M represents cos(magnetic latitude of mirror point)
eqn = M.^6 + 3*M*(sind(pitcheq)).^4 == 4*(sind(pitcheq)).^4; 
assume(M <= 1);
assume(M >= 0);
assume(M == 'Real'); %These assumes don't appear to be working!
COS2M = vpasolve(eqn,M,0.5);% Finds mirror point latitude prediction from wequatorial pitch angle
mirrormaglatitudeprediction = double( 180/pi*acos(sqrt(COS2M(2))) );  %Solution #2 is only real positive solution.
fprintf(1,'Mirror point angle prediction = %4g %c. ', mirrormaglatitudeprediction, char(176));
fprintf(logfile,'Mirror point angle prediction =\t%4g\t%c. ', mirrormaglatitudeprediction, char(176));
fprintf(trajectoryfile,'Mirror point angle prediction =\t%4g\t%c. ', mirrormaglatitudeprediction, char(176));
if lower(mdisc) == 'y'
    fprintf(logfile, 'Axisymmetric magnetodisc.\n');
elseif llimit < 1
   fprintf(logfile, 'Pure dipole field.\n');
else 
    fprintf(logfile, 'IGRF, l limit =\t%d\t\n', llimit);
end

numberofgyros = ceil (bounces*Bounceperiodpred/EquatorialGyroperiodpred); %approximate number of gyroperiods in 'bounces'
tmax = numberofgyros*EquatorialGyroperiodpred;   % max time. Should be a large number of gyroperiods
tsteps = numberofgyros*stepspergyro;% Overall number of time steps
fprintf(logfile,'Tmax =\t%4g\ts.\n',tmax);
fprintf(1,'Tmax =\t%4g\ts.\t%d\tsteps per gyro-orbit.\nSimulation started at %s\t',tmax, stepspergyro, datestr(now));

tic; %simulation timer started.
tspans = linspace(0,tmax,tsteps); % time: start, limit, divisions.

Rinitial = sqrt(R0*(R0+2*EquatorialGyroradiuspredRP*cosd(initial_angle))+EquatorialGyroradiuspredRP.^2);
deltaz = EquatorialGyroradiuspredRP*cosd(initial_angle)*sind(mlat0);
if initial_angle > 180
    deltaz = -deltaz;
end
zinitial = z0 + deltaz;
maglatinitial = asind(zinitial/Rinitial);
rcylinitial = sqrt(Rinitial.^2 - zinitial.^2);
rcyl0 = sqrt(R0.^2 - z0.^2);
rgcyl = sqrt(EquatorialGyroradiuspredRP.^2-deltaz.^2);
%deltalongitudeold = abs(atand(EquatorialGyroradiuspredRP*sind(initial_angle)/(R0+EquatorialGyroradiuspredRP*cosd(initial_angle))))
deltalongitude = real(acosd((rcylinitial.^2+rcyl0.^2-rgcyl.^2)/rcyl0/rcylinitial/2));
maglonginitial = mlong0+deltalongitude;
previousmaglongitude = maglonginitial;

[xinitial, yinitial, zinitial] = polartocart(Rinitial, maglatinitial, maglonginitial);
Vr = speed0*sind(pitcheq)*(-sind(initial_angle)) ;
Vlat = speed0*cosd(pitcheq);
Vlong = speed0*sind(pitcheq)*(cosd(initial_angle));%*360/pi
vzinitial = Vr*sind(mlat0)+Vlat*cosd(mlat0);
vxinitial = (R0*Vr-z0*vzinitial)/rcyl0*cosd(mlong0)-Vlong*sind(mlong0);
vyinitial = (R0*Vr-z0*vzinitial)/rcyl0*sind(mlong0)+Vlong*cosd(mlong0);
                                   
        %initial position and velocity components
posvelmatrix = [xinitial; yinitial; zinitial;vxinitial;vyinitial; vzinitial];
rmax = sqrt(xinitial.^2 + yinitial.^2 + zinitial.^2); %maximum & minimum values of variou metrics initialised.
rmin = sqrt(xinitial.^2 + yinitial.^2 + zinitial.^2);
rcylmax = rcylinitial;
rcylmin = rcylinitial;
xmax = xinitial;
ymax = yinitial;
zmax = zinitial;
xmin = xinitial;
ymin = yinitial;
zmin = zinitial;
maglatbiggest = atand(zinitial/rcylinitial);
maglatmax = atand(zinitial/rcylinitial);
maglatmin = atand(zinitial/rcylinitial);
maglongmin = atand(yinitial/xinitial);
maglongmax = atand(yinitial/xinitial);
Lmin = L0;
Lmax = L0;
Gyroradiusmax = 0;
Gyroradiusmin = 2*EquatorialGyroradiuspred;
Gyroperiodmax = 0;
Gyroperiodmin = 2*EquatorialGyroperiodpred ;
mutotal = 0;
musquarestotal = 0;
mumax = 0;
mumin = 1000;
framecount = 0;

% ODE solver options
absolute_tolerance =  1e-6;
relative_tolerance =  1e-6;
options = odeset('AbsTol',absolute_tolerance,'RelTol',relative_tolerance,'Stats','on','vectorized','on');
switch solver_choice,
    case 1, solvername = 'ode45';
    case 2, solvername = 'ode15s';
end
fprintf(resultsfile, '%s\t%e\t%e\t', solvername, absolute_tolerance, relative_tolerance);
fprintf(logfile, 'Solver =\t%s\t, relative tolerance =\t%e\t, absolute tolerance =\t%e\t\n', solvername, absolute_tolerance, relative_tolerance);
fprintf(trajectoryfile, 'Solver =\t%s\t, relative tolerance =\t%e\t, absolute tolerance =\t%e\t. %s\t%s\n', solvername, absolute_tolerance, relative_tolerance, trajectory_filename, gyroperiods_filename);
fprintf(gyroperiodsfile, 'Solver =\t%s\t, relative tolerance =\t%e\t, absolute tolerance =\t%e\t. %s\t%s\n', solvername, absolute_tolerance, relative_tolerance, trajectory_filename, gyroperiods_filename);
fprintf(1, 'Solver = %s\t, relative tolerance = %e\t, absolute tolerance = %e\n', solvername, absolute_tolerance, relative_tolerance);
fprintf(1,'\nPlease wait.\n\n');
done20percent = false;
done40percent = false;
done60percent = false;
done80percent = false;
switch solver_choice,
    case 1, [t1,X1] = ode45(@dynamic3d,tspans,posvelmatrix,options); solver=solvername; %non-stiff
    case 2, [t1,X1] = ode15s(@dynamic3d,tspans,posvelmatrix,options); solver=solvername; %stiff
end
% save trajectory t1,X1 and input arguments (metadata)
save(mat_filename,'t1','X1','R0','mlat0','mlong0','KE0','pitcheq','particle','planet','llimit','mdisc','initial_angle','outputthinning','solver_choice','Method','bounces')
mumean = mutotal/framecount;
musd = sqrt(musquarestotal/framecount-mumean.^2);
rcylplotstart = rcylmin; %These values set the plotting ranges.
rcylplotend = rcylmax;
rplotstart = rmin;
rplotend = rmax;
xplotstart = xmin;
xplotend = xmax;
yplotstart = ymin;
yplotend = ymax;
zplotstart = zmin;
zplotend = zmax+0.000001;
Lplotstart = Lmin;
Lplotend = Lmax+0.000001;

[numberoftimeslots, n] = size(X1);
MLongplot = zeros(numberoftimeslots,1);  %sets up new array for longitude vs. time plot.
Lplot = zeros(1,numberoftimeslots);  %sets up new array for longitude vs. L plot.
MLatplot = zeros(numberoftimeslots,1);  %sets up new array for longitude vs. time plot.
MLongplot(1) = maglong(X1(1,1), X1(1,2), maglonginitial);
for S = 2: 1: numberoftimeslots
    MLongplot(S) = maglong(X1(S,1), X1(S,2), MLongplot(S-1));
end
%for S = 1: 1: numberoftimeslots
    MLatplot(:) = atand(X1(:,3)./(sqrt(X1(:,1).^2+X1(:,2).^2)));
    Lplot(:) = sqrt(X1(:,1).^2+X1(:,2).^2+X1(:,3).^2)./(cosd(atand(X1(:,3)./sqrt(X1(:,1).^2+X1(:,2).^2)))).^2;
%end
                      %title for figures
figuretitle = strcat(num2str(KE0,3),' keV',' ',particlename,' in',' ',planetnameposessive,' magnetosphere. initial R = ',num2str(R0,3),', initial mag latitude = ',num2str(mlat0,3), char(176),' initial mag longitude = ',num2str(mlong0,3), char(176),', equatorial pitch angle = ',num2str(pitcheq,3),char(176),', initial phase angle = ',num2str(initial_angle,3),char(176), '. ',num2str(bounces,2),' bounce periods, ', num2str(stepspergyro,2),' steps per gyro-orbit.');
if lower(mdisc) == 'y'
    figuretitle = strcat(figuretitle,' Magnetodisc.');
elseif llimit < 1
    figuretitle = strcat(figuretitle,' Pure dipole field.');
else 
    figuretitle = strcat(figuretitle,' IGRF, l limit ', num2str(llimit));
end    
[Bounceperiod, sineamplitude, Bounceperiod_error, sineamplitude_error,  convergence] = sinefit(t1, atand(X1(:,3)./(sqrt(X1(:,1).^2+X1(:,2).^2))), 1/Bounceperiodpred, figuretitle);
if (convergence == 0)
    Bounceperiod = Bounceperiodpred;
    sineamplitude = mirrormaglatitudeprediction;
    Bounceperioderror_old = 0;
else
    Bounceperioderror_old = 1/ceil(bounces*2)*0.006*Bounceperiod; %Bounceperiod precision is imporoved by increading the number of bounces.
    if Bounceperioderror_old < (t1(numberoftimeslots)  - t1(1))/(numberoftimeslots-1)/ceil(bounces*2) %Bounceperiod precision is limited by the time resolution.
        Bounceperioderror_old = (t1(numberoftimeslots)  - t1(1))/(numberoftimeslots-1)/ceil(bounces*2);
    end
end
Bounceperiod_discrepancy =  (Bounceperiod- Bounceperiodpred)/Bounceperiodpred*100;

margin = ceil(Bounceperiod/40/Gyroperiodmin)/2;  % *Marker C*
Bounceperiodaverage = Bounceperiod; %initialising.
Bounceperiodspread = 0;
BounceperiodSD = 0;
Bounceperiodmax = Bounceperiod; 
Bounceperiodmin = Bounceperiod;
numberofmirrorpoints = ceil((tmax-Bounceperiod/4)/Bounceperiod*2);

if numberofmirrorpoints > 0
    if upper(Method) == 'A'
        mirrorpoints = zeros(numberofmirrorpoints,3);
        for W = 1:1:numberofmirrorpoints
            counter = 0;
            mlatitudefluctuationsaroundmirrorpoint = 0;
            maglattempmax = -2*maglatbiggest;
            maglattempmin = 2*maglatbiggest;
            for S = 2:1:numberoftimeslots
                if abs(t1(S) - Bounceperiod*(2*(W-1)+1)/4) <= margin*Gyroperiodmin  %this code adds magnetic latitude values of points very close
                    counter = counter + 1;               % to first mirror point and calculates the mean.
                    maglat = atand(X1(S,3)./(sqrt(X1(S,1).^2+X1(S,2).^2)));
                    mlatitudefluctuationsaroundmirrorpoint = mlatitudefluctuationsaroundmirrorpoint + maglat;
                    if maglat > maglattempmax
                        maglattempmax = maglat;
                    end
                    if maglat < maglattempmin
                        maglattempmin = maglat;
                    end
                end
            end
            mirrorpoints(W,1) = Bounceperiod*(2*(W-1)+1)/4;
            mirrorpoints(W,2) = abs(mlatitudefluctuationsaroundmirrorpoint/counter); % mean
            mirrorpoints(W,3) = abs(maglattempmax- maglattempmin)/4;
            fprintf(mirrorpointsfile,'\n%4g\t%4g\t%4g', mirrorpoints(W,1), mirrorpoints(W,2), mirrorpoints(W,3));
        end
        mirrormaglatitude = mean(mirrorpoints(:,2)); %alternative determination of mirror point latitude.
        mirrormaglatitude_error = mean(mirrorpoints(:,3));
        mirrormaglatitude_sd = std(mirrorpoints(:,2));
        mirrormaglatitudemin = min(mirrorpoints(:,2));
        mirrormaglatitudemax = max(mirrorpoints(:,2));
    else  %Method B
        counter = 0;
        bounceperiodsfile = fopen(bounceperiods_filename,'w');
        fprintf(bounceperiodsfile, '%s\t%s\t%s\tinitial L =\t%.4g\tinitial magnetic latitude =\t%4g\tinitial magnetic longitude =\t%4g\tenergy  =\t%.4g\tkeV\tequatorial pitch angle =\t%.4g\tdeg\tmagnetic moment =\t%.4g\tTkm3\t%3g\tinitial phase angle = %4g\tdeg\tbounce periods\t%d\tsteps per gyro.\t', datestr(now), planetname,particlename, L0, mlat0, mlong0, KE0, pitcheq, magmoment, initial_angle, bounces, stepspergyro);
        if lower(mdisc) == 'y'
           fprintf(bounceperiodsfile, 'Magnetodisc.\t\n');
        elseif llimit < 1
            fprintf(bounceperiodsfile, 'Dipole.\t\n');
        else 
            fprintf(bounceperiodsfile, 'IGRF, l limit =\t%d\n', llimit);
        end
        fprintf(bounceperiodsfile,'Time (s)\tBounceperiod (s)');
        lastmaxormintime = -Bounceperiod/4; %Roughly where the previous mirror point would be if went backwards in time.
        padding = 10;
        for S = padding-1:1:numberoftimeslots-padding
            if abs(MLatplot(S)) > abs(MLatplot(S-1)) && abs(MLatplot(S)) > abs(MLatplot(S+1)) && t1(S)- lastmaxormintime > Bounceperiod/4%Maximum or minimum found.
                if abs(MLatplot(S)) > abs(mean(MLatplot(S-padding:S-1))) && abs(MLatplot(S)) > abs(mean(MLatplot(S+1:S+padding))) %Further conditons to exclude 'mini' maxima or minima.
                    counter = counter + 1;               % Counting the number of minima.
                end
            end
        end
        numberofmirrorpoints = counter;
        bounceperiods = zeros(numberofmirrorpoints-2,3);
        mirrorpoints = zeros(numberofmirrorpoints,3);
        counter = 0;
        for S = padding-1:1:numberoftimeslots-padding
            if abs(MLatplot(S)) > abs(MLatplot(S-1)) && abs(MLatplot(S)) > abs(MLatplot(S+1)) && t1(S)- lastmaxormintime > Bounceperiod/4%Maximum or minimum found.
                if abs(MLatplot(S)) > abs(mean(MLatplot(S-padding:S-1))) && abs(MLatplot(S)) > abs(mean(MLatplot(S+1:S+padding))) %Further conditons to exclude 'mini' maxima or minima.
                    counter = counter + 1;               % to first mirror point and calculates the mean.
                    mirrorpoints(counter,1) = t1(S);
                    mirrorpoints(counter,2) = abs(MLatplot(S));
                    mirrorpoints(counter,3) = 0;
                    fprintf(mirrorpointsfile,'\n%4g\t%4g\t%4g', mirrorpoints(counter,1), mirrorpoints(counter,2), mirrorpoints(counter,3));
                    lastmaxormintime = t1(S);
                end
            end
        end
        for D = 2:1:counter-1
            bounceperiods(D-1,1)= mirrorpoints(D,1);
            bounceperiods(D-1,2)= mirrorpoints(D+1,1)-mirrorpoints(D-1,1);
            fprintf(bounceperiodsfile,'\n%4g\t%4g', mirrorpoints(D,1), mirrorpoints(D+1,1)-mirrorpoints(D-1,1));
        end
        Bounceperiodaverage = mean(bounceperiods(:,2));
        Bounceperiodspread = max(bounceperiods(:,2))-min(bounceperiods(:,2));
        BounceperiodSD = std(bounceperiods(:,2));
        Bounceperiodmin = min(bounceperiods(:,2));
        Bounceperiodmax = max(bounceperiods(:,2));
        mirrormaglatitude = mean(mirrorpoints(:,2)); %alternative determination of mirror point latitude.
        mirrormaglatitude_error = mean(mirrorpoints(:,3));
        mirrormaglatitudemin = min(mirrorpoints(:,2));
        mirrormaglatitudemax = max(mirrorpoints(:,2));
        mirrormaglatitude_sd = std(mirrorpoints(:,2));
    end
else
    mirrormaglatitude = maglatbiggest; %Just to prevent error messages if no maxima or minima exist.
    mirrormaglatitude_error = 0;
    mirrormaglatitude_sd = 0;
    mirrormaglatitudemin = maglatbiggest;
    mirrormaglatitudemax = maglatbiggest;
end
mirrormaglatitude_discrepancy =  (mirrormaglatitude - mirrormaglatitudeprediction)/mirrormaglatitudeprediction*100;

Lminimacount = 0;
firstLminimumtime = 0;
firstLminimummarker = 1;
lastLminimummarker = numberoftimeslots;
if tmax < Bounceperiod/10
    numberofminsforgyroperiodcalc = floor(tmax/EquatorialGyroperiodpred); % *Marker A*
else  numberofminsforgyroperiodcalc = floor(Bounceperiod/EquatorialGyroperiodpred*0.03);
end
if numberofminsforgyroperiodcalc < 2
    numberofminsforgyroperiodcalc = 2;
elseif numberofminsforgyroperiodcalc > 11
    numberofminsforgyroperiodcalc = 11;
end
nthLminimumtime = 0;
nthLminimummarker = 1;
J = 0;
for S = 2: 1: numberoftimeslots-1  %This section of code is used to calculate the number of minima.
    if Lplot(S) <= Lplot(S-1) && Lplot(S) <= Lplot(S+1) %local minimum of L found.
        Lminimacount = Lminimacount +1;
        if firstLminimumtime == 0
            firstLminimumtime = t1(S);
            firstLminimummarker = S;
        end
        if Lminimacount < numberofminsforgyroperiodcalc+1
            nthLminimumtime = t1(S);
            nthLminimummarker = S;
            J = Lminimacount; 
        end
        lastLminimummarker = S;
    end
end
if Lminimacount > 1 
    gyroperiods = zeros(Lminimacount-1,2);
else  gyroperiods = zeros(2,2);
end
minimumtime = 0;
if Lminimacount > 1 %we need at lest 2 minima to measure gyroperiod.
    Lminimacount = 0;
    for S = 2: 1: numberoftimeslots-1 %This section of code is used to calculate gyroperiods from the minimum positions.
        if Lplot(S) <= Lplot(S-1) && Lplot(S) <= Lplot(S+1) %local minimum of L found.
            Lminimacount = Lminimacount +1;
            if Lminimacount > 1
                gyroperiods(Lminimacount-1,2) = t1(S)-minimumtime;
            gyroperiods(Lminimacount-1,1) = (t1(S)+minimumtime)/2; %First column in gyroperiods is the time of the midpoint between the first two minima.
            fprintf(gyroperiodsfile, '%6g\t%5g\n',gyroperiods(Lminimacount-1,1), gyroperiods(Lminimacount-1,2));
            end
        minimumtime = t1(S); %used for calculating next gyroperiod.
        end
    end
end
time_resolution = (t1(numberoftimeslots) - t1(1))/(numberoftimeslots-1);%average time interval.
EquatorialGyroperiod = 0;
if (Method =='B') 
    firstbounceperiod = bounceperiods(1,2);
else
    firstbounceperiod = Bounceperiod;
end
if tmax > firstbounceperiod/2 + numberofminsforgyroperiodcalc/2*Gyroperiodmin 
                %This checks whether are likely to be sufficient minima close to the... 
    tstart = firstbounceperiod/2 - numberofminsforgyroperiodcalc/2*Gyroperiodmin;  %first return to the equatorial plane.
    tstartmarker = 1;
    for D = 2:1:Lminimacount-1
        if abs(gyroperiods(D,1)- tstart) <= abs(gyroperiods(D-1,1)- tstart)
            tstartmarker = D;
        end
    end
    for E = tstartmarker : 1 : tstartmarker+numberofminsforgyroperiodcalc-2 %takes average of several gyro periods...
        EquatorialGyroperiod = EquatorialGyroperiod + gyroperiods(E,2)/(numberofminsforgyroperiodcalc-1); %close to particle's first return to equatorial plane.
    end
    EquatorialGyroperiod_error = time_resolution/(numberofminsforgyroperiodcalc-1)/2;
elseif nthLminimumtime > 0 && firstLminimumtime > 0 %Takes average of first J-1 gyroperiods...
    EquatorialGyroperiod = (nthLminimumtime-firstLminimumtime)/(J-1); %to calculate gyroperiod if tmax is less than one half a bounce cycle.
    EquatorialGyroperiod_error = (nthLminimumtime-firstLminimumtime)/(nthLminimummarker-firstLminimummarker)/(J-1);
else EquatorialGyroperiod = EquatorialGyroperiodpred; 
    EquatorialGyroperiod_error = 0;%Shouldn't happen if things working correctly.
end
Gyroperiod_discrepancy = (EquatorialGyroperiod-EquatorialGyroperiodpred)/EquatorialGyroperiodpred*100;

lastxyplanecrossingtimeadj = -1;
if tmax >= firstbounceperiod/2
    if Method == 'A'
        lastxyplanecrossingtime = tmax-rem(tmax,Bounceperiod/2);
           %lastxyplanecrossingtime is the most appropriate endpoint for
           %finding the drift gradient and subsequently the drift period.
    else
        for G=numberoftimeslots-1:-1:1
            if (MLatplot(G)< 0 && MLatplot(G+1)> 0)
                lastxyplanecrossingtime = t1(G);
                break;
            end
        end
    end
    for T = 2: 1: Lminimacount-1
        if abs(gyroperiods(T,1)- lastxyplanecrossingtime)< abs(gyroperiods(T-1)- lastxyplanecrossingtime)
            lastxyplanecrossingtimeadj = gyroperiods(T,1); %time of L minimum that closest to time of last x-y plane crossing point is used to determine drift period.
        end
    end
    if lastxyplanecrossingtimeadj >= 0
        for S = 2: 1: numberoftimeslots
            if abs(t1(S)- lastxyplanecrossingtimeadj)< abs(t1(S-1)- lastxyplanecrossingtimeadj)
                lastxyplanecrossingtimemarker = S;
            end
        end
        [maglonggradient, G] = polyfit(t1(firstLminimummarker:lastxyplanecrossingtimemarker), MLongplot(firstLminimummarker:lastxyplanecrossingtimemarker),1); %Calculates gradient of mag lontitude vs time plot in seconds/degree.
        Driftorbitperiod = abs(360/maglonggradient(1)); 
        deltat = t1(lastLminimummarker)-t1(firstLminimummarker);
            
            %Drift orbit period can also be determined by calculating a
            %simple between the start and end points.
        Driftorbitperiodalt = abs(360/(MLongplot(lastxyplanecrossingtimemarker) -MLongplot(firstLminimummarker))*(t1(lastxyplanecrossingtimemarker)-t1(firstLminimummarker)));
        discrepancyplus = abs(360/(MLongplot(lastxyplanecrossingtimemarker-1) -MLongplot(firstLminimummarker+1))*(t1(lastxyplanecrossingtimemarker-1)-t1(firstLminimummarker+1))- Driftorbitperiodalt );
        discrepancyminus = abs(360/(MLongplot(lastxyplanecrossingtimemarker+1) -MLongplot(firstLminimummarker-1))*(t1(lastxyplanecrossingtimemarker+1)-t1(firstLminimummarker-1))- Driftorbitperiodalt );
        Driftorbitperiodalt_error = (discrepancyplus + discrepancyminus)/2;
    else
        [maglonggradient, G] = polyfit(t1(1:numberoftimeslots), MLongplot(1:numberoftimeslots),1); %Last resort if program cannot find any L-minima
        Driftorbitperiod = abs(360/maglonggradient(1));
        deltat = tmax-t1(1); %should = tmax.
        Driftorbitperiodalt = 0;
        Driftorbitperiodalt_error = 0;
    end
else % if tmax is less than half a bounce period.
    [maglonggradient, G] = polyfit(t1(firstLminimummarker:lastLminimummarker), MLongplot(firstLminimummarker:lastLminimummarker),1); %Calculates gradient of mag lontitude vs time plot in seconds/degree.
    Driftorbitperiod = abs(360/maglonggradient(1));
    deltat = t1(lastLminimummarker)-t1(firstLminimummarker);
    Driftorbitperiodalt = 0;
    Driftorbitperiodalt_error = 0;
end
[a_longitude,deltalong] = polyval(maglonggradient,deltat,G) ; %a_longtiude should be very close to MLongplot(lastLminimummarker);
      %deltalong is uncertianity on a_longitude from polyval calculation.
Driftorbitperiod_error = Driftorbitperiod - abs(360*deltat/(a_longitude- maglonggradient(2)+deltalong));
      %converting deltalong into uncertainity on drift orbit period.

fprintf(1,'Intercept from magnetic longitiude fit = %4g, Polyval results = %4g, period error = %4g\n', maglonggradient(2),a_longitude,Driftorbitperiod_error);
Driftorbitperiod_discrepancy = (Driftorbitperiod-Driftorbitperiodpred)/Driftorbitperiodpred*100;
Driftorbitperiodalt_discrepancy = (Driftorbitperiodalt-Driftorbitperiodpred)/Driftorbitperiodpred*100;

figure('Name', figuretitle)
  
subplot(241), % plot trajectory rcyl vs. z
plot(sqrt(X1(:,1).^2+X1(:,2).^2),X1(:,3),'k-'), 
xlabel('Rcyl (Rp)','FontName','Cambria'),
ylabel('z (Rp)','FontName','Cambria'),
title(sprintf('Rcyl vs. z'));
uicontrol('Style', 'text',...
       'String', figuretitle,... 
       'Units','normalized',...
       'Position', [0.075 0.98 0.85 0.02]); %4 paras are text boc start point x, strart point y, width, height.
set(gca,'xlim',[rcylplotstart rcylplotend])
set(gca,'zlim',[zplotstart zplotend])

subplot(242),
plot(X1(:,1),X1(:,2)), % plot trajectory x vs y
xlabel('x (Rp)','FontName','Cambria'),
ylabel('y (Rp)','FontName','Cambria'),
title(sprintf('x vs. y'));
set(gca,'xlim',[xplotstart xplotend])
set(gca,'ylim',[yplotstart yplotend])
  
subplot(243),
plot(sqrt(X1(:,1).^2+X1(:,2).^2+X1(:,3).^2),atand(X1(:,3)./sqrt(X1(:,1).^2+X1(:,2).^2))), % plotr vs. mag lat 
xlabel('R (Rp)','FontName','Cambria'),
ylabel('mag lat (deg)'),
title(sprintf('R vs. latitude'));
set(gca,'xlim',[rplotstart rplotend])
set(gca,'ylim',[maglatmin maglatmax+0.01])

subplot(244),     
plot(t1,atand(X1(:,3)./(sqrt(X1(:,1).^2+X1(:,2).^2)))), % plot magnetic latitude vs t
xlabel('time (s)','FontName','Cambria'),
ylabel('mag lat (deg)','FontName','Cambria'),
title(sprintf('mag latitude vs. time'));
set(gca,'xlim',[0 tmax])
set(gca,'ylim',[maglatmin maglatmax+0.01])
  
subplot(245),
plot(t1,MLongplot), % plot magnetic longitude vs t
xlabel('time (s)','FontName','Cambria'),
ylabel('mag long (deg)'),
title(sprintf('mag longitude vs. time'));
set(gca,'xlim',[0 tmax])
set(gca,'ylim',[maglongmin maglongmax])
  
subplot(246),
plot(t1,Lplot);
xlabel('time (s)','FontName','Cambria'),
ylabel('L'),
title(sprintf('L vs. time'));
set(gca,'xlim',[0 tmax])
set(gca,'ylim',[Lplotstart Lplotend])

  subplot(247),
plot(gyroperiods(:,1),gyroperiods(:,2),'k.');
xlabel('Time (s)','FontSize',10),
ylabel('Gyroperiod (s)'),
title(sprintf('Gyroperiod vs. time'));
set(gca,'xlim',[0 tmax])

Harmonics = fourier(t1,atand(X1(:,3)./(sqrt(X1(:,1).^2+X1(:,2).^2))));

figure('pos',[0 40 1280 900],'Name', figuretitle);
sphere(40)
hold on
uicontrol('Style', 'text',...
       'String', figuretitle,... 
       'Units','normalized',...
       'Position', [0.075 0.98 0.85 0.02]);
line(X1(:,1),X1(:,2),X1(:,3), 'color', 'black');
view([45 60])
set(gca,'xlim',[-rcylmax rcylmax]);
set(gca,'ylim',[-rcylmax rcylmax]);
xlabel('X (Rp)','FontSize',10,'FontName','Cambria');
ylabel('Y (Rp)','FontSize',10,'FontName','Cambria');
zlabel('Z (Rp)','FontSize',10,'FontName','Cambria');
zpluslimit = zmax;
if zmax < 1
    zpluslimit = 1;
end
zminuslimit = zmin;
if zmin > -1
    zminuslimit = -1;
end
set(gca,'zlim',[zminuslimit zpluslimit]);
hold off

figure('pos',[0 40 1280 900],'Name', figuretitle);
hold on
uicontrol('Style', 'text',...
       'String', figuretitle,... 
       'Units','normalized',...
       'Position', [0.075 0.98 0.85 0.02]);
line(sqrt(X1(:,1).^2+X1(:,2).^2),X1(:,3), 'color', 'black');
set(gca,'ylim',[zminuslimit zpluslimit], 'Fontsize', 14);
set(gca,'xlim',[0 rcylmax], 'Fontsize', 14);
xlabel('Rcyl (Rp)','FontSize',18,'FontName','Cambria');
ylabel('Z (Rp)','FontSize',18,'FontName','Cambria');
phi = linspace(0,360,361);
cx = cosd(phi);  %draws a cricle (planet's surface).
cy = sind(phi);
plot(cx,cy,'b-','LineWidth',2);
hold off

simtime = toc; %timer stopped. 
  
fprintf(logfile, 'Max magnetic latitude =\t%3g\t%c, min R =\t%4g\tRp, Max rcyl =\t%.6g\tRp, min rcyl =\t%.6g\tRp,\n',maglatbiggest, char(176),rmin, rcylmax, rcylmin);
fprintf(logfile, 'max z =\t%.6g\tRp, min z =\t%.6g\tRp, min mag lat =\t%.3g\t%c, max mag long =\t%.3g\t%c, min mag long =\t%.3g\t%c, max L =\t%5g\t, min L =\t%5g,\nmax gyroperiod =\t%5g\t, min gyroperiod =\t%5g\t, max gyroradius =\t%5g\t, min gyroradius =\t%5g\t.\nTime elapsed =\t%3g\ts.\n',zmax, zmin, maglatmin, char(176),maglongmax,char(176), maglongmin, char(176),Lmax, Lmin, Gyroperiodmax, Gyroperiodmin, Gyroradiusmax, Gyroradiusmin, simtime);
fprintf(logfile, 'Gyroperiod from plot =\t%5g\ts, gyroperiod discrepancy =\t%5g\t%%, error on Gyroperiod =\t%5g\ts, drift orbit period from plot =\t%5g\ts, drift orbit period discrepancy =\t%5g\t%%, error on drift orbit period =\t%5g\tdrift orbit period from plot alternative = \t%5g\ts, drift orbit period alt discrepancy =\t%5g\t%%, error on drift orbit period alt =\t%5g\n', EquatorialGyroperiod, Gyroperiod_discrepancy, EquatorialGyroperiod_error, Driftorbitperiod, Driftorbitperiod_discrepancy, Driftorbitperiod_error, Driftorbitperiodalt, Driftorbitperiodalt_discrepancy, Driftorbitperiodalt_error);
fprintf(logfile, 'bounceperiod from plot =\t%5g\ts, bounceperiod discrepancy =\t%5g\t%%, error on bounceperiod from plot =\t%5g\ts.\n', Bounceperiod, Bounceperiod_discrepancy, Bounceperiod_error);
fprintf(logfile, 'mirror point latitude from plot alternative =\t%3g\t%c, mirror point latitude alternative discrepancy =\t%3g\t%%, error on mirrorpointlatitude alternative =\t%3g\t%c\tnumber of mirror points\t%d\tmirror point latitudes sd\t%4g\t%s.\n', mirrormaglatitude, char(176), mirrormaglatitude_discrepancy, mirrormaglatitude_error, char(176), numberofmirrorpoints, mirrormaglatitude_sd, char(176));
fprintf(logfile, 'Mean magnetic moment = \t%4g\t, standard devitation =\t%4g\t, max magnetic moment = \t%4g\t, min magnetic medium =\t%4g\t. Periods from FFT (Hz) - \t%5g \t%5g\t%5g\t%5g\n', mumean, musd, mumax, mumin, Harmonics(1), Harmonics(2),Harmonics(3),Harmonics(4));
fprintf(1, 'Max magnetic latitude = %3g degrees, min R = %4g  Rp, Max rcyl = %.6g Rp, min rcyl = %.6g Rp,\n',maglatbiggest, rmin, rcylmax, rcylmin);
fprintf(1, 'max z =\t%.6g\tRp, min z =\t%.6g\tRp, min mag lat =\t%.3g\t%c, max mag long =\t%.3g\t%c, min mag long =\t%.3g\t%c, max L =\t%5g\t, min L =\t%5g\nmax gyroperiod =\t%5g\t, min gyroperiod =\t%5g\t, max gyroradius =\t%5g\t, min gyroradius =\t%5g\t.\nTime elapsed =\t%3g\ts.\n',zmax, zmin, maglatmin, char(176),maglongmax,char(176), maglongmin, char(176),Lmax, Lmin, Gyroperiodmax, Gyroperiodmin, Gyroradiusmax, Gyroradiusmin, simtime);
fprintf(1, 'Gyroperiod from plot =\t%5g\ts, gyroperiod discrepancy =\t%5g\t%%, error on Gyroperiod =\t%5g\ts, drift orbit period from plot =\t%5g\ts, drift orbit period discrepancy =\t%5g\t%%, error on drift orbit period =\t%5g\tdrift orbit period from plot alternative = \t%5g\ts, drift orbit period alt discrepancy =\t%5g\t%%, error on drift orbit period alt =\t%5g\n', EquatorialGyroperiod, Gyroperiod_discrepancy, EquatorialGyroperiod_error, Driftorbitperiod, Driftorbitperiod_discrepancy, Driftorbitperiod_error, Driftorbitperiodalt, Driftorbitperiodalt_discrepancy, Driftorbitperiodalt_error);
fprintf(1, 'bounceperiod from plot =\t%5g\ts, bounceperiod discrepancy =\t%5g\t%%, error on bounceperiod from plot =\t%5g\ts.\n', Bounceperiod, Bounceperiod_discrepancy, Bounceperiod_error);
fprintf(1, 'mirror point latitude from plot alternative =\t%3g\t%c, mirror point latitude alternative discrepancy =\t%3g\t%%, error on mirrorpointlatitude alternative =\t%3g\t%c\tnumber of mirror points\t%d\tmirror point latitudes sd\t%4g\t%s.\n', mirrormaglatitude, char(176), mirrormaglatitude_discrepancy, mirrormaglatitude_error, char(176), numberofmirrorpoints, mirrormaglatitude_sd, char(176));
fprintf(1, 'Mean magnetic moment = \t%4g\t, standard devitation =\t%4g\t, max magnetic moment = \t%4g\t, min magnetic medium =\t%4g\t. Periods from FFT (Hz) - \t%5g \t%5g\t%5g\t%5g\n', mumean, musd, mumax, mumin, Harmonics(1), Harmonics(2),Harmonics(3),Harmonics(4));
if upper(Method) == 'B' 
    fprintf(1, 'Mean bounce period =\t%4g\tMax bounce period =\t%4g\tMin bounce period =\t%4g\tBounce period SD =\t%4g\n', Bounceperiodaverage, Bounceperiodmax, Bounceperiodmin, BounceperiodSD);
    fprintf(logfile, 'Mean bounce period =\t%4g\tMax bounce period =\t%4g\tMin bounce period =\t%4g\tBounce period SD =\t%4g\n', Bounceperiodaverage, Bounceperiodmax, Bounceperiodmin, BounceperiodSD);
end
fprintf(1, '%4g Bounces per orbit.\n', Driftorbitperiodalt/Bounceperiod);
fprintf(1, 'Framecount = %i\n', framecount);
fprintf(resultsfile, '%5g\t%5g\t%5g\t%4g\t%4g\t%5g\t%.6g\t%.6g\t%.6g\t%.6g\t',B0, speed0kps,EquatorialGyroradiuspred, EquatorialGyroradiuspred/(R0*planeteqradius),Lossconeangle, maglatbiggest, rmax, rmin, rcylmax, rcylmin);
fprintf(resultsfile, '%.6g\t%.6g\t%5g\t%5g\t%5g\t%5g\t%.6g\t%.6g\t%5g\t%5g\t%5g\t%5g\t%5g\t',zmax, zmin, maglatmax, maglatmin, maglongmax,maglongmin,Lmax, Lmin, Gyroperiodmax, Gyroperiodmin, Gyroradiusmax, Gyroradiusmin,simtime);
fprintf(resultsfile, '%5g\t%5g\t%5g\t%5g\t%4g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t', EquatorialGyroperiodpred, EquatorialGyroperiod,max(gyroperiods(:,2)),min(gyroperiods(:,2)), EquatorialGyroperiod_error, Gyroperiod_discrepancy, Bounceperiodpred, Bounceperiod, Bounceperiod_error, Bounceperiod_discrepancy, Bounceperiodaverage, Bounceperiodmax, Bounceperiodmin, Bounceperiodspread, BounceperiodSD);
fprintf(resultsfile, '%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%4g\t', Driftorbitperiodpred,Driftorbitperiod, Driftorbitperiod_error, Driftorbitperiod_discrepancy, Driftorbitperiodalt, Driftorbitperiodalt_error, Driftorbitperiodalt_discrepancy, 360*tmax/(maglongmax-maglongmin));
fprintf(resultsfile, '%5g\t%5g\t%5g\t%5g\t%5g\t%5g\t%4g\t%4g\t',mirrormaglatitudeprediction, mirrormaglatitude, mirrormaglatitude_error, numberofmirrorpoints, mirrormaglatitude_sd, mirrormaglatitudemax, mirrormaglatitudemin, mirrormaglatitude_discrepancy);
fprintf(resultsfile, '%4g\t%4g\t%4g\t%4g\t%4g\t%5g\t%5g\t%5g\t%5g\t', mumean, musd, mumax, mumin, musd/mumean, Harmonics(1),Harmonics(2),Harmonics(3),Harmonics(4));
fprintf(resultsfile, '%5g\t%5g\t%4g\t',time_resolution,numberoftimeslots, tmax);
if rmin < 1;% + 100/planetpolarradius lost = true;
   fprintf(1, '\aParticle is lost to atmosphere!\n'); %This means that lowest radius reached is inside planet's surface.
   fprintf(logfile, 'Particle is lost to atmosphere!\n');
   fprintf(resultsfile, 'Lost to atmosphere\t');
elseif rmax > R0+20*EquatorialGyroperiodpred %This means that particle is energetic enough to escape magnetosphere altogether.
   fprintf(1, '\aParticle escaped to infinity!\n');
   fprintf(logfile, 'Particle escaped to infinity!\n');
   fprintf(resultsfile, 'Escaped to infinity\t');
else
   fprintf(resultsfile, 'Trapped\t'); 
   fprintf(logfile, 'Particle is trapped.\n');
   fprintf(1, 'Particle is trapped.\n');
end
fprintf(resultsfile, '%s\t%s\n', trajectory_filename, gyroperiods_filename);
fprintf(logfile, '%s\t%s\t%s\t%s\n', trajectory_filename, gyroperiods_filename, mirrorpoints_filename, bounceperiods_filename);

% 3D ODE function dr/dt=f(r) where r=(x,y,z,vx,vy,vz)
function dr_dt = dynamic3d(t,r)
    if t > 0.2*tmax && done20percent == false
        done20percent = true;
        fprintf(1,'%s 20%% done,\t',datestr(now)); %Output gives a guide as to how far into the simulation the computer has got.
    end
    if t > 0.4*tmax && done40percent == false
        done40percent = true;
        fprintf(1,'%s 40%% done,\t',datestr(now));
    end
    if t > 0.6*tmax && done60percent == false
        done60percent = true;
        fprintf(1,'%s 60%% done,\t',datestr(now));
    end
    if t > 0.8*tmax && done80percent == false
        done80percent = true;
        fprintf(1,'%s 80%% done.\n',datestr(now));
    end
    
    [x,y,z,vx,vy,vz] = deal(r(1,:),r(2,:),r(3,:),r(4,:),r(5,:),r(6,:)) ;% split coordinates
    rcyl = sqrt(x.^2 + y.^2);
    
    R = sqrt(x.^2 + y.^2 + z.^2);
    vmod = sqrt(vx.^2 + vy.^2+ vz.^2);
    vmodkps = vmod*planeteqradius;
    Gamma = 1./sqrt(1-(vmod.^2).*(planeteqradius/c).^2);
    
    maglatitude = atand(z/rcyl);
    maglongitude = maglong(x,y,previousmaglongitude);
    previousmaglongitude = maglongitude;
    L = R / (cosd(maglatitude)).^2;
    dr = EquatorialGyroradiuspred/10;
    xplus = [x+dr/planeteqradius y z];
    xminus = [x-dr/planeteqradius y z];
    yplus = [x y+dr/planeteqradius z];
    yminus = [x y-dr/planeteqradius z];
    zplus = [x y z+dr/planeteqradius];
    zminus = [x y z-dr/planeteqradius];
    if lower(mdisc) == 'y'
        [Br, Bcolat] = MDiscField(R0,(90-maglatitude)*pi/180);
        Br = 10000*Br;
        Blat = -10000*Bcolat;
        B_z = Br*sind(maglatitude)+Blat*cosd(maglatitude);
        B_x = (R*Br-z*B_z)/rcyl*cosd(maglongitude);
        B_y = (R*Br-z*B_z)/rcyl*sind(maglongitude);
    elseif llimit < 1  % calculates magnetic field components assuming pure dipole.
        Bvector = getBvectordipole([0,0,magmoment], dipole_vector , {x,y,z},planeteqradius);
        [B_x,B_y,B_z] = deal(Bvector{:});
    else % calculates magnetic field components using IGRF coefficients.
        Vxplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(xplus,maglongitude), Earthaverageradius, llimit);
        Vxminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(xminus,maglongitude),Earthaverageradius, llimit);
        B_x = (Vxplus - Vxminus)/dr/2;
        Vyplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(yplus,maglongitude),Earthaverageradius, llimit);
        Vyminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(yminus,maglongitude),Earthaverageradius, llimit);
        B_y = (Vyplus - Vyminus)/dr/2;
        Vzplus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(zplus,maglongitude),Earthaverageradius, llimit);
        Vzminus = getMpot_REM(magmoment/Earthmagmoment,cart_to_polar(zminus,maglongitude),Earthaverageradius, llimit);
        B_z = (Vzplus - Vzminus)/dr/2;
    end   
    Bmod = sqrt(B_x.^2 + B_y.^2+ B_z.^2);
    pitch = 180-acosd((vx.*B_x+vy.*B_y+vz.*B_z)./(vmod.*Bmod));  %*Marker D* instantaneous pitch angle calculated.
    vperpkps = vmodkps*sind(pitch);
    Mu = Gamma*particlemass*(vperpkps*1000).^2/Bmod/2*10000; %*Marker E* instantaneous mu (magnetic moment) calculated.
    %Mu = kinetic_energy_perp./Bmod*abs(Electroncharge)*10000*1e6; 
    mutotal = mutotal + Mu;
    musquarestotal = musquarestotal + Mu.^2;
    Gyroperiod = 2 * pi * particlemass*Gamma/(abs(particlecharge)*Bmod/10^4); %Instantaneous gyroperiod calculated.
    Gyroradius = abs(Gyroperiod*vmodkps*sind(pitch)/2/pi); %Instataneous gyro-radius in km calculated.

%    if rem(framecount,outputthinning) == 0  %data written to trajectory file. Not all framses recorded (determined by outputthinning.
%        fprintf(trajectoryfile,'%5g\t%5g\t%5g\t%5g\t%5g\t%4g\t%4g\t%4g\t%5g\t%4g\t%4g\t%4g\t%4g\t%3g\t%4g\t%.5g\t%5g\t%6g\t%5g\n',t,x,y,z,R,rcyl, maglatitude,maglongitude,L, vx,vy,vz,Mu, pitch, Bmod, Gamma, vmodkps, Gyroperiod, Gyroradius);
%    end
    % Lorentz Force Z V x B / m
    ax = particlecharge*(vy.*B_z - vz.*B_y)./10000/particlemass./Gamma;
    ay = particlecharge*(vz.*B_x - vx.*B_z)./10000/particlemass./Gamma;
    az = particlecharge*(vx.*B_y - vy.*B_x)./10000/particlemass./Gamma;

    if abs(maglatitude) > maglatbiggest  %updates maglatbiggest if necessary.
        maglatbiggest = abs(maglatitude);
    end
    if Mu > mumax  
        mumax = Mu;
    end
    if Mu < mumin  
        mumin = Mu;
    end
    if maglatitude > maglatmax  
        maglatmax = maglatitude;
    end
    if maglatitude < maglatmin  
        maglatmin = maglatitude;
    end
    if  maglongitude > maglongmax
        maglongmax = maglongitude;
    end
    if  maglongitude < maglongmin
        maglongmin = maglongitude;
    end
    if rcyl > rcylmax
        rcylmax = rcyl;
    end
    if rcyl < rcylmin
        rcylmin = rcyl;
    end
    if R > rmax
        rmax = R;
    end
    if R < rmin
        rmin = R;
    end
    if x > xmax
        xmax = x;
    end
    if x < xmin
        xmin = x;
    end
    if y > ymax
        ymax = y;
    end
    if y < ymin
        ymin = y;
    end
    if z > zmax
        zmax = z;
    end
    if z < zmin
        zmin = z;
    end
    if L > Lmax
        Lmax = L;
    end
    if L < Lmin
        Lmin = L;
    end
    if Gyroradius > Gyroradiusmax
        Gyroradiusmax = Gyroradius;
    end
    if Gyroradius < Gyroradiusmin
        Gyroradiusmin = Gyroradius;
    end
    if Gyroperiod > Gyroperiodmax
        Gyroperiodmax = Gyroperiod;
    end
    if Gyroperiod < Gyroperiodmin
        Gyroperiodmin = Gyroperiod;
    end
    framecount = framecount+1; %framecoutn incremented.
    dr_dt = [vx;vy;vz;ax;ay;az]; % return dr/dt
end

function magnetic_longitude = maglong(XX, YY, previous_maglong)%
    if XX == 0 && YY == 0 %If point lies on polar axis.
        ml=0;
    else   
        ml = atan2d(YY,XX);
        if particlecharge > 0
            while abs(previous_maglong - ml) > 270
                ml = ml + 360;
            end
        else 
            while abs(previous_maglong - ml) > 270
                ml = ml - 360;
            end
        end
    end
   magnetic_longitude = ml;
end

function polar = cart_to_polar(cart, amaglongitude)
    R = norm(cart); % R is modulus of position vector R2 = x2 + y2 + z2.
    phi = maglong(cart(1),cart(2),amaglongitude);
    theta = asind(cart(3)/R);
    polar = [R theta phi];
end
function [X, Y, Z] = polartocart(R, latitude, longitude)
    Z = R*sind(latitude);  %initial potition in cartesian, x0, y0 & zo.
    X = sqrt(R.^2 - Z.^2)*cosd(longitude);
    Y = X*tand(longitude);
    if longitude == 90
        Y = sqrt(R.^2 - Z.^2);
    elseif longitude == 270
        Y = -sqrt(R.^2 - Z.^2);
    end
end
function [period, amplitude, period_error, amplitude_error, convergence] = sinefit(t,maglat,frequency_guess, figuretitle)
    p0 = [max(abs(maglat));frequency_guess;0]; % initial guess
    dp = [1; 1; 1]; % fit 3 parameters
    % do not fit phase parameter p(3)
    stol=1e-8;  % tolerance
    niter=100;  % max iter
    minstep = [0; 0; 0];
    maxstep = [Inf; Inf; Inf];
    options = [minstep, maxstep];

    t = t(:);
    maglat = maglat(:);
    weight = ones(size(t)); % weight (=1 everywhere by default.

    global verbose;
    verbose = [0 0];
    [f, p, kvg, iter, corp, covp, covr, stdresid, Z, r2, ss] = ...
        leasqr(t, maglat, p0, @F, stol, niter, weight, dp, @dFdp, options);
%    fprintf(1,'method3 fi  %10g / kvg %d / iter %d / r2 %g / ss %g\n',...
%        frequency_guess, iter, kvg, r2, ss)

    figure ('Name', figuretitle);
    plot(t, maglat, '.', t, F(t, p0), 'g:', t, f,'k-');
    uicontrol('Style', 'text',...
       'String', figuretitle,... 
       'Units','normalized',...
       'Position', [0.075 0.98 0.85 0.02]); 
   xlabel('Time (s)','FontName','Cambria','FontSize',11),
    ylabel('Magnetic latitude (degrees)','FontName','Cambria','FontSize',11),
    title(sprintf('Magnetic latitude vs. time with fitted sinusoid'));
    
    period = 1/p(2);
    amplitude = p(1);
    frequency_error = sqrt(covp(2,2)); % *Marker B* Error calculated from covariance matrix.
    period_error = period.^2*frequency_error;
    amplitude_error = sqrt(covp(1,1)); % Error calculated from covariance matrix.
    convergence = kvg;
end
% function to fit
function y = F(t,p)
    y = p(1)*sin(2*pi*p(2)*t+p(3));
end

% and derivative
function y = dFdp(t,f,p,dp,func)
    %fprintf(1,'called dFdp(t,[%e,%e,%e]\n', p(1),p(2),p(3))
    y = [sin(2*pi*p(2)*t+p(3)), ...
         p(1)*2*pi*t.*cos(2*pi*p(2)*t+p(3)), ...
        p(1)*cos(2*pi*p(2)*t+p(3))];
end

function periods = fourier(t,maglat)
    ns = length(t); % Number of samples
    dt = t(2)-t(1); % sampling time
    fs = 1/dt; % sampling frequency
    fmax = fs/2; % Nyquist frequency
    fex = fmax/300;

    if mod(ns,2)==1, % zero-centred frequency axis: odd ns
        frq = [-floor(ns/2):floor(ns/2)]/floor(ns/2)*fmax;
    else % even ns
        frq = [-floor(ns/2):floor(ns/2)-1]/floor(ns/2)*fmax;
    end

    % alternative formulation without testing ns parity 
    % total time interval
    T = max(t) - min(t);
    % nearest lower even number
    neven = 2*floor(ns/2);
    % more accurate formulation
    frq1 = [floor(-(ns-1)/2):floor((ns-1)/2)]/T*((ns-1)/neven);
    % check they are equal to machine precision
    frq = frq1;
    Spc = fftshift(abs(fft(maglat))); % zero-centred spectra amplitude
    numberofmaxima = 0;
    for S = 2:1:size(Spc)-1%ceil(size(Spc)/2)-1
        if Spc(S) > Spc(S-1) && Spc(S) > Spc(S+1) && frq(S) >= 0
            numberofmaxima = numberofmaxima+1; %finds how many maxima exist in plot.
        end
    end
    frequencies = zeros(numberofmaxima,2);
    periods = zeros(4,1);
    maximacounter = 0;
    for S = 2:1:size(Spc)-1%ceil(size(Spc)/2)-1
        if Spc(S) > Spc(S-1) && Spc(S) > Spc(S+1) && frq(S) >= 0
            maximacounter = maximacounter+1;
            frequencies(maximacounter,2)= frq(S); 
            frequencies(maximacounter,1)= Spc(S);
                %extracts frequency associated with each maximum.
        end
    end
   frequencies = sortrows(frequencies);
    periods(1) = 1./frequencies(numberofmaxima,2); %Frequency that contains the most power (bounce period).
    if size(frequencies) > 1, periods(2)=  1./frequencies(numberofmaxima-1,2); %Frequency that contains the second most power.
    end
    if size(frequencies) > 2, periods(3)=  1./frequencies(numberofmaxima-2,2); %Frequency that contains the third most power. 
    end
    if size(frequencies) > 3, periods(4)=  1./frequencies(numberofmaxima-3,2); %Frequency that contains the third most power. 
    end
    subplot(248), % plot spectra and signal frequency
    plot(frq, Spc,'-'), 
   set(gca,'xlim',[-fmax fmax])
    xlabel('frequency (Hz)','FontName','Cambria'),
    ylabel('FT(mag lat)','FontName','Cambria'),
    title(sprintf('FFT(mag latitude vs. time)'));
end

end
