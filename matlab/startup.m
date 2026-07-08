%startup.m

% settings *needed* to avoid error below when quitting
% Error using quit
% Unrecognized function or variable 'Settings'.
%
% workarounf for bug in R2021b, fixed in update 2
% settings;

if exist('groot') == 5,
  h = groot;
else, % prior to 2014b
  h = 0;
end

% from R2025a needed to get previous release behaviour
set(h, 'DefaultFigureToolBar', 'figure');
set(h, 'DefaultFigureMenuBar', 'figure');
set(h, 'DefaultFigurePosition', [360 58 560 420]);

%set(h,'DefaultFigurePosition',[550 500 505 405]);
set(h,'DefaultFigurePaperType','a4letter');
set(h,'DefaultFigurePaperUnits','centimeters');
%set(h,'DefaultFigurePaperPosition',[2.5 2.5 16.5 24.2]/2.53807);
set(h,'DefaultFigurePaperPosition',[2.5 2.5 16.5 24.2]);
%set(h,'DefaultFigureNumberTitle','off');

set(h,'DefaultAxesFontName','times');
set(h,'DefaultAxesFontWeight','normal');
set(h,'DefaultAxesFontSize',14);
set(h,'DefaultAxesTickDir','in');
set(h,'DefaultAxesXGrid','on');
set(h,'DefaultAxesYGrid','on');
set(h,'DefaultAxesZGrid','on');

% trick for R2014b onward to get same behaviour as earlier
set(h,'DefaultAxesGridLineStyle',':');
set(h,'DefaultAxesGridAlpha',1);
set(h,'DefaultAxesGridColor',[0.1 0.1 0.1]);

set(h,'DefaultTextFontName','times');
%set(h,'DefaultTextFontWeight','demi');
set(h,'DefaultTextFontSize',16);

set(h,'DefaultLineMarkersize',2);
set(h,'DefaultLineLineWidth',.2);

set(h,'DefaultSurfaceLineWidth',.65);

%p = {[getenv('HOME') '/research/codes']};
p = {'../data'};

for i=1:length(p),
  if exist(p{i},'dir'),
    addpath(p{i});
  end
end

clear h i p

