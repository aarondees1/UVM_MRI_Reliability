% %#######################################################################
%
%                  * Femur Bone RaNGe Male Program *
%
%          M-File which reads the male femur data to get the femur
%     coordinate range and generate a "standard" quadrilateral mesh in
%     the cylindrical coordinate system.  The resulting "standard" mesh
%     is saved to MAT file fb_rngm.mat.
%
%     NOTES:  1.  The MAT data files must be in the
%             "BM_Male_Femur_MatFiles" directory.
%
%     14-Aug-2015 * Mack Gardner-Morse
%

%#######################################################################
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
for i = 1:ns
    i_str=int2str(i);

    for j=1:3
        if j==1
            fstr=sd(i).FFE.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end
    end

        bone = load(fullfile(rdir,'Bone',[fstr csnam]));
%
% Load Data
%
   load(fullfile(pnams{i},fnam));
%
% Plot Bone Data
%
   figure;
   view(3);
   hold on;
   plt_datsl(datfb,'b.-');
   axis equal;
   xlabel('X','FontSize',12,'FontWeight','bold');
   ylabel('Y','FontSize',12,'FontWeight','bold');
   zlabel('Z','FontSize',12,'FontWeight','bold');
   title(fstr,'FontSize',24,'FontWeight','bold', ...
          'Interpreter','none');

% Transform to Cylindrical Coordinate System
%
   xyz = cell2mat(datfb);
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
%
% Find Overall Minimums and Maximums
%
   itst = tzrn<tzrmn;
   tzrmn(itst) = tzrn(itst);
   itst = tzrx>tzrmx;
   tzrmx(itst) = tzrx(itst);
%
    
end
%
% Get Range and Make a Rectangular Grid
%
tzrmn(1) = tzrmn(1)*180/pi;            % Convert to degrees
tzrmx(1) = tzrmx(1)*180/pi;            % Convert to degrees
tzmn = floor(tzrmn(1:2));              % Get degrees/millimeters closest to -infinity
tzmx = ceil(tzrmx(1:2));               % Get degrees/millimeters closest to +infinity
%
t = (tzmn(1):2:tzmx(1))';              % Theta range in 2 degrees increments
nt = size(t,1);
% z = (tzmn(2):2:tzmx(2))';              % Z range in 2 mm increments
z = (tzmn(2):tzmx(2)+2)';              % Z range in 1 mm increments
nz = size(z,1);
[T,Z] = meshgrid(t,z);  % Rectangle grid of square quadrilaterals
tq = T(:);              % Grid points as a vector
zq = Z(:);              % Grid points as a vector
nq = size(tq,1);
%
% FROM QUAD CONNECT
quad = (1:nz-1)';
quad = [quad quad+nz quad+nz+1 quad+1];
quad = repmat(quad,nt-1,1);
addcol = repmat(nz*(0:nt-2),(nz-1)*4,1);
addcol = reshape(addcol,4,(nt-1)*(nz-1))';
quad = quad+addcol;
%
% Triangular Grid
%
ne = size(quad,1);      % Number of elements
%
trig = [quad(:,1:3) quad(:,1) quad(:,3:4)]';
trig = reshape(trig,3,2*ne)';
%
% Save Data
%
save fb_rngm.mat nq quad tq zq trig tzrmn tzrmx;
%
return