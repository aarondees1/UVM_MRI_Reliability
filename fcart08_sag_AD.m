clear;
close all;
clc;


csnam = '_femurCS.mat';
cnam = '_femurCart.mat';
%
%
% Get All Subject Subdirectories
%
div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');
cdir = fullfile(rdir, 'Cartilage');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
load(fullfile(rdir, 'Full Tibia Full Femur Subject Files.mat'));
ns=size(sd,1);
%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1);    % 0 (false) for left/1 (true) for right
kids = repmat(blanks(ns)',1,5);        % Knee IDs
tribfs = cell(ns,1);    % Sagittal bone triangular mesh connectivity
xyzbfs = cell(ns,1);    % Sagittal bone point coordinates
xyzmnl = zeros(ns,3);   % Minimum individual lateral femur coordinates
xyzmxl = zeros(ns,3);   % Maximum individual lateral femur coordinates
xyzmnm = zeros(ns,3);   % Minimum individual medial femur coordinates
xyzmxm = zeros(ns,3);   % Maximum individual medial femur coordinates
xyzs = zeros(ns,3);     % Range of proximal femur outline
rf = zeros(ns,3);       % Radii of femur condyles

scan(1,:)='FFE';
scan(2,:)='RHO';
scan(3,:)='T2S';
%

% Initialize Overall Minimums and Maximums
%
tzrmn = zeros(1,3);     % Minimums
tzrmx = zeros(1,3);     % Maximums
%
hpi = pi/2;             % Half pi
dpi = 2*pi;             % Double pi
%
% Loop through MAT Files
%
%% Load Bone Data and Generate min/max and range of Grid
for j=1:3
    for i = 1:ns
    i_str=int2str(i);

        if j==1
            fstr=sd(i).FFE.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end
    
        bone = load(fullfile(rdir,'Bone',[fstr csnam]));


   eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.trimesh = bone.trifb;']);
   eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.points = bone.xyzfb;']);

   cart = load(fullfile(rdir,'Cartilage',[fstr cnam]));

   eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.trimesh = cart.trifc;']);
   eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.points = cart.xyzfc;']);
%
% Plot Bone Data
%
   figure;
   view(3);
   hold on;
   plt_datsl(bone.datfb,'b.-');
   axis equal;
   xlabel('X','FontSize',12,'FontWeight','bold');
   ylabel('Y','FontSize',12,'FontWeight','bold');
   zlabel('Z','FontSize',12,'FontWeight','bold');
   title(fstr,'FontSize',24,'FontWeight','bold', ...
          'Interpreter','none');

% Transform to Cylindrical Coordinate System
%
   xyz = cell2mat(bone.datfb);
%
   [th,r,z] = cart2pol(xyz(:,1),xyz(:,3),xyz(:,2));
   id = find(th>hpi);                  % Greater than pi/2
   th(id) = th(id)-dpi;                % Minus 2*pi
   tzr = [th,z,r];
%
% Get Minimums and Maximums
%
   tzrn = min(tzr);
   tzrx = max(tzr);

      itst = tzrn<tzrmn;
   tzrmn(itst) = tzrn(itst);
   itst = tzrx>tzrmx;
   tzrmx(itst) = tzrx(itst);

    end
    for i=1:ns
        i_str=int2str(i);
        gnam = [scan(j,:) '_fgrid08_.mat'];
        eval(['sd(' i_str ').' scan(j,:) '.femur.grid = gnam;']);
        gnam = fullfile(gdir,gnam);   % Grid MAT file
    end

    save gnam nq quad tq zq trig tzrmn tzrmx;

end