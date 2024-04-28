%#######################################################################
%
%        * Femur REGionS from BioMarkers Female Knees Program *
%
%          M-File which reads the standard grid and plots the trochlea
%     and condyle regions from the MS-Excel spreadsheet
%     "Partitions6_f3.xlsx".
%
%     NOTES:  1.  MAT file fb_rngf.mat and fcthkf.mat must be in the
%             current directory.
%
%             2.  The MS-Excel spreadsheet, Partitions6_f3.xlsx, must
%             be in the current path or directory.
%
%     28-Jan-2021 * Mack Gardner-Morse
%

%#######################################################################
%
% Get "Standard" Grid
%
load fb_rngf.mat quad tq zq;
%
% Get Radius
%
load fcthkf.mat rq;
idn = isnan(rq);
rq(idn) = 0;
rq = sum(rq)./sum(~idn);     % Get average radius
rq = rq';
idn = isnan(rq);
rq(idn) = 0;
%
% Get XYZ Coordinates from Cylindrical Coordinates
%
deg2rad = pi/180;
%
trz = [tq rq zq];
%
idb = find(trz(:,1)>180);
trz(:,1) = trz(:,1)*deg2rad;
trz(idb,1) = trz(idb,1)-2*pi;
%
% Convert to Cartesian Coordinates
%
[x,y,z] = pol2cart(trz(:,1),trz(:,2),trz(:,3));
x(idn) = NaN;
y(idn) = NaN;
z(idn) = NaN;
%
% Get Regions
%
xnam = 'Partitions6_f3.xlsx';
idreg = xlsread(xnam);
%
% Plot Regions
%
figure;
orient landscape;
colormap([0.2 0.2 1; 1 0 0]);
%
patch(x(quad'),y(quad'),z(quad'),idreg(quad'));
%
view(180,-50);
axis equal;
axis off;
axis tight;
%
colorbar;
%
title({'Biomarkers Females'; ['Trochlea (1) and Condyle (0) ', ...
       'Regions']},'FontSize',20,'FontWeight','bold');
%
print -dpsc2 -r600 -fillpage fregs_bmf3.ps
%
return