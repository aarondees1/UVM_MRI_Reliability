%#######################################################################
%
%             * Femur Cartilage THicKnesses Male Program *
%
%          M-File which reads femur cartilage and subchondral bone MAT
%     files and calculates the cartilage thicknesses.  The cartilage
%     thickness data is saved to the MAT file:  fcthkm.mat.
%
%     NOTES:  1.  For male biomarkers femur data.
%
%             2.  The femur data is scaled based on the male tibia
%             coordinate system data in the "\bm male tibia cs"
%             directory in MAT files ending with "*CS.mat".
%
%             3.  The MAT data files must be in the "Biomarker MAT
%             files" directory.
%
%             4.  The cartilage and subchondral bone MAT files must end
%             with "EAbc7tsc.mat".
%
%             5.  MAT file fb_rngm.mat must be in the current directory.
%
%             6.  The M-files car_thk5.m, in_tri2d.m. nod2tri.m,
%             nod_norm.m, plane_fit.m, tri_fix2.m, tri_norm.m, tsect4.m
%             and xprod.m must be in the current path or directory.
%
%     20-Aug-2015 * Mack Gardner-Morse
%

%#######################################################################
%close all; clear all; clc;
%
% Get "Standard" Grid
%
load fb_rngm.mat;
%
% Make Figure Handles
%
hf1 = figure;
hf2 = figure;
% autoArrangeFigures(1,2);
%
% Data Directory and File Names
%
ddir = 'Femur';          % Data directory
d = dir(fullfile(ddir,'*CS.mat'));
fnams = deblank(char(d.name));
nf = size(fnams,1);     % Number of files
%
% Initialize Variables
%
% iprt = false;
iprt = true;
rad2deg = 180/pi;       % Radians to degrees
%
% Get Tibia CS File Directory
%
%
if ismac
  tdir = 'Tibia';
else
  tdir = 'Tibia';
end
tnams = dir(fullfile(tdir,'tcart08_1.mat'));
tnams = char(tnams.name);
%
% Get Tibia Widths
%
nft = size(tnams,1);
tpw = zeros(nft,1);
for k = 1:nft
   fnamt = tnams(k,:);
   load(fullfile(tdir,fnamt));
   tpw(k) = max(xyzpt(:,2))-min(xyzpt(:,2));
end
tnams = tnams(:,1:5);
tpw_mn = mean(tpw);
%
if exist('fcthkm.mat','file')
    load fcthkm;
else
    cthk = zeros(nf,nq);               % Cartilage thicknesses
    xyzic = zeros(nf,nq,3);            % Cartilage intersection
    rq = zeros(nf,nq);                 % Radial coordinates
    zqs = zeros(nf,nq);                % Polar Z coordinates (scaled)
end
%
% Loop through MAT Files
%
for kf = 1:nf           % Loop through MAT files
%
   fnam = fnams(kf,:);                 % File name
   idot = strfind(fnam,'.');
   fstr = fnam(1:idot-1);              % File name w/o extension
%
% Load Data
%
%    d = dir([fstr(1:8) '*.mat']);    % Check for modified data
%    fn = sortrows(char(d.name));
%    nfsz = size(fn,1);
%    fn = fn(nfsz,:);     % Use the most modified data
%
   load(fullfile(ddir,fnam));
   
   trifb = mk_tri4f(datfb);       %taken from fcomba.m for scaled data to solve 'xyzsb' and 'trisb' AD 3/3/23
   xyzfb = cell2mat(datfb);
   trifb = tri_fix2(trifb,xyzfb);

   datfc = comb_dat(datlct,datmct,dattct);
   trifc = mk_tri4f(datfc);
   xyzfc = cell2mat(datfc);
   trifc = tri_fix2(trifc,xyzfc);

   
%
% Get Cylindrical Coordinates
%
   [t,r,z] = cart2pol(xyzfb(:,1),xyzfb(:,3),xyzfb(:,2));
   id = find(t>pi/2);
   t(id) = t(id)-2*pi;
   td = t*rad2deg;      % Angles in degrees
%
% Check Triangle Normals and Make Orientations +R
% Not Necessary For Thickness Calculation, Add this Step to MK_TRI4F.M?
%
%    [nt,nz,nr] = tri_norm(trisb,[td z r]);
%
%    idn = find(nr<0);
%    if ~isempty(idn)
%      trisb(idn,:) = trisb(idn,[1 3 2]);
%    end
%
% Scale Polar "Z" Coordinates
%
   iz = 1;
   zsc = tpw(iz)./tpw_mn;
   zs = zsc*zq;
%
% Interpolate to Scaled "Standard" Grid
%
   F = TriScatteredInterp(td,z,r);
   r = F(tq,zs);
   tzrg = [tq zs r];
   it = in_tri2d(trifb,[td z],tzrg(:,1:2));
   tzrg(~it,:) = NaN;
   r(~it) = NaN;
   rq(kf,:) = r';
   zqs(kf,:) = zs';
%
% Back to Cartesian Coordinates
%
   [xg,zg,yg] = pol2cart(tzrg(:,1)/rad2deg,tzrg(:,3),tzrg(:,2));
   xyzg = [xg yg zg];
%
% Improve the Mesh
%
   trig = tri_fix2(trig,xyzg);
%
% Calculate Cartilage Thicknesses
%
   [ct,bi] = car_thk5(trifc,xyzfc,trig,xyzg,iprt,fstr,[],hf1);
%car_thk5(trib,xyzb,tric,xyzc,iplt,fstr,cmprt,hf1);
% Save Thicknesses and Print Plot
%
   cthk(kf,:) = ct';
   xyzic(kf,:,:) = bi;
   if iprt
     if kf==1
       print('-dpsc2','-r300','-fillpage',hf1,'fcthkm.ps');
     else
       print('-dpsc2','-r300','-fillpage','-append',hf1,'fcthkm.ps');
     end
     if kf==1
       print('-dpsc2','-r300','-fillpage',hf2,'fcthkm2.ps');
     else
       print('-dpsc2','-r300','-fillpage','-append',hf2,'fcthkm2.ps');
     end
%      pause;
     print('-dpdf','-r300','-fillpage',hf1,['fcthkm1_' fstr(1:5)]);
     print('-dpdf','-r300','-fillpage',hf2,['fcthkm2_' fstr(1:5)]);
     clf(hf1);
     clf(hf2);
   end
%
   fprintf(1,'File number %i completed!\n',kf);
%
end                     % File loop kf
%
% Save Data in MAT File
%
save fcthkm cthk fnams rq xyzic zqs
%
return