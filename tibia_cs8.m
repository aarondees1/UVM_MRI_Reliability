function [xyzo,xyzax,aspect,widt,height,xyzpt,prox_tib_x,prox_tib_y,prox_tib_ml,prox_tib_ap] = tibia_cs8(faxial,leg,iplt); % Define function to calculate tibia coordinate system
%TIBIA_CS8 Reads Osirix regions-of-interest comma separated files to
%          calculate tibia based coordinate systems (CS).
%
%          [XYZO,XYZAX] = TIBIA_CS8(FAXIAL,LEG) given the filename of
%          an Osirix regions-of-interest comma separated of transverse
%          sections of the tibia, FAXIAL, and an integer 1 (or logical
%          true) for a right knee and 0 (or logical false) for a left
%          knee, LEG, calculates tibia coordinate system based on the
%          centroids from the most distal section of the tibia visible
%          on the MRI and the most proximal section of the tibia.  The
%          most posterior points on the tibia at the transverse section
%          through the PCL form the second axis.  The function returns 
%          the X, Y and Z center of the tibia coordinate system, XYZO,
%          and unit vectors for the X, Y and Z axes in the columns of
%          matrix XYZAX.
%
%          [XYZO,XYZAX] = TIBIA_CS8(FAXIAL,LEG,IPLT) generates plots as
%          it calculates the axes.
%
%          [XYZO,XYZAX,ASPECT,WIDT,HEIGHT] = TIBIA_CS8(FAXIAL,LEG)
%          returns the aspect ratio, ASPECT, of the width of the tibia
%          plateau, WIDT, to the height between the distal and proximal
%          outlines of the tibia plateau, HEIGHT.
%
%          [XYZO,XYZAX,ASPECT,WIDT,HEIGHT,XYZPT,PROX_TIB_ML,PROX_TIB_AP] = 
%          TIBIA_CS8(FAXIAL,LEG)
%          returns the aspect ratio, ASPECT, of the width of the tibia
%          plateau, WIDT, to the height between the distal and proximal
%          outlines of the tibia plateau, HEIGHT. Also returns med/lat 
%          dimension PROX_TIB_ML and ap dimension PROX_TIB_AP. This data
%          will be used for scaling the grids used in cartilage thickness
%          measurements 
%
%
%          NOTES:  1.  The filename must be in single row character
%                  array.
%
%                  2.  The names of the regions of interest (ROIs) must
%                  follow a naming convention.
%
%                  3.  The name of the transverse ROI of the posterior
%                  tibia must be "postaxis".
%
%                  4.  The names of the distal transverse ROIs of the
%                  tibias must be "dist_tibia" and the proximal
%                  sections must be named "prox_tibia".
%
%                  5.  The Matlab files li_clos.m, rd_roi4.m, rotxyz.m
%                  and tri_area.m must be in the current directory or
%                  path.
%
%                  6.  This is a cleaned up version of tibia_cs7.m with
%                  a loop to find improved posterior points of the
%                  tibia.
%
%          07-May-2019 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<2) % Check if fewer than 2 inputs are provided
  error([' *** ERROR in tibia_cs8:  At least an input', ...
         ' filename and left/right knee identifier are required!']); % Throw error for insufficient inputs
end
%
if (nargin<3) % Check if plot flag is missing
  iplt = false; % Set default plot flag to false
end
%
% File Name
%
nrow = size(faxial,1); % Get number of rows in filename
if nrow~=1 % Check if filename is not a single row
  error([' *** ERROR in TIBIA_CS8:  Filename must be in single', ...
         ' row character array!']); % Throw error for invalid filename format
end
%
% Get Current Figure
%
if iplt % Check if plotting is enabled
  gcf; % Get current figure
  gca; % Get current axes
  hold on; % Enable hold for multiple plots
end
%
% Posterior Points Parameters
%
ptol = pi/180; % Define one-degree tolerance in radians
ang = zeros(2,1); % Initialize rotation angle array
dang = 10*ptol; % Set initial angle difference
%
% Initialize Variables
%
xme = zeros(2,1); % Initialize X coordinates for posterior points
yme = zeros(2,1); % Initialize Y coordinates for posterior points
zme = zeros(2,1); % Initialize Z coordinates for posterior points
%
% Read File
%
filenam = deblank(faxial); % Remove trailing blanks from filename
%
% Get ROI Data
%
roi = rd_roi6(filenam); % Read ROI data from file
nroi = size(roi,1); % Get number of ROIs
%
% Check for DIST_TIBIA ROI
%
iuse = false; % Initialize flag for using anataxis ROI
%
roinams = lower(char(roi.name)); % Convert ROI names to lowercase
idd = strmatch('dist_tibia',roinams); % Find 'dist_tibia' ROI
if isempty(idd) % Check if 'dist_tibia' is not found
  idd = strmatch('anataxis',roinams); % Find 'anataxis' ROI
  if isempty(idd) % Check if 'anataxis' is not found
    error([' *** ERROR in tibia_cs8:  No distal tibia ROI found in', ...
           ' data file!']); % Throw error for missing distal tibia ROI
  end
  iuse = true; % Set flag to use anataxis ROI
end
%
% Loop through ROIs
%
for l = 1:nroi % Loop over ROIs
%
   roiname = roi(l).name; % Get ROI name
   roinaml = lower(roiname); % Convert ROI name to lowercase
%
   dat = roi(l).data'; % Transpose ROI data
   nslice = size(dat,1); % Get number of slices
%
% Get Distal Tibia Centroids
%
   if strcmp(roinaml,'anataxis')&iuse % Check for 'anataxis' ROI and flag
%
% Get Digitized Data
%
     xyzc = zeros(nslice,3); % Initialize centroid array
%
     for n = 1:nslice % Loop over slices
        xyz = dat{n}; % Get slice data
        npts = size(xyz,1); % Get number of points
%
        xyzm = mean(xyz); % Compute mean of points
        xyzp = [xyz; xyzm]; % Append mean to points
%
% Get Centroids
%
        tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)']; % Create triangle connectivity
        tri = [tri; [npts+1 npts 1]]; % Add final triangle
%
        [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri); % Compute triangle areas and centroids
        xyzc(n,:) = sum(repmat(at,1,3).*cgt)./sum(at); % Compute weighted centroid
%
        if iplt % Check if plotting is enabled
          plot3(xyzc(n,1),xyzc(n,2),xyzc(n,3),'bs', ...
                'LineWidth',1,'MarkerSize',8); % Plot centroid
        end
     end
%
% Get and Plot Most Distal Centroid
%
     [~,idz] = min(xyzc(:,3)); % Find index of most distal centroid
     xyzc1 = xyzc(idz,:); % Get most distal centroid
     if iplt % Check if plotting is enabled
       xyz = dat{idz}; % Get corresponding slice data
     end
%
     if iplt % Check if plotting is enabled
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','MarkerSize',8, ...
             'LineWidth',1); % Plot slice points
%
       plot3(xyzc1(:,1),xyzc1(:,2),xyzc1(:,3),'bs','MarkerSize',8, ...
             'LineWidth',1); % Plot distal centroid
     end
%
   end
%
% Get Most Distal Tibia Centroid
%
   if strcmp(roinaml,'dist_tibia') % Check for 'dist_tibia' ROI
%
     xyz = cell2mat(dat); % Convert data to matrix
     npts = size(xyz,1); % Get number of points
%
     xyzm = mean(xyz); % Compute mean of points
     xyzp = [xyz; xyzm]; % Append mean to points
%
% Get and Plot Centroid
%
     tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)']; % Create triangle connectivity
     tri = [tri; [npts+1 npts 1]]; % Add final triangle
     [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri); % Compute triangle areas and centroids
     xyzc1 = sum(repmat(at,1,3).*cgt)./sum(at); % Compute weighted centroid
%
     if iplt % Check if plotting is enabled
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-', ...
             'LineWidth',1); % Plot slice points
%
       plot3(xyzc1(:,1),xyzc1(:,2),xyzc1(:,3),'bs', ...
             'LineWidth',1); % Plot centroid
     end
%
   end
%
% Posterior PCL
%
   if strcmp(roinaml,'postaxis')||strcmp(roinaml,'post_axis') % Check for 'postaxis' or 'post_axis' ROI
%
% Plot Posterior PCL
%
     if nslice>2 % Check for excessive slices
       warning([' *** WARNING in tibia_cs8:  Too many POSTAXIS ', ...
                'slices in ' faxial '!']); % Issue warning
     end
%
     nsz = zeros(2,1); % Initialize slice size array
     xyz_post = zeros(2,3); % Initialize posterior points array
     for n = 1:2 % Loop over slices
         xyz = dat{n}; % Get slice data
         nsz(n) = size(xyz,1); % Get number of points
         if iplt % Check if plotting is enabled
           plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k^-', ...
                 'MarkerSize',8,'LineWidth',1); % Plot slice points
         end
%
% Find Most Posterior Points on the Posterior Tibia Plateau at the PCL
%
         [~,id] = max(xyz(:,2)); % Find index of most posterior point
         xyz_post(n,:) = xyz(id,:); % Store posterior point
     end
%
% Combine Both Posterior Sides
%
     ns = sum(nsz); % Compute total number of points
     idp{2} = nsz(1)+1:ns; % Set indices for second side
     idp{1} = 1:nsz(1); % Set indices for first side
%
     xyz = cell2mat(dat); % Convert data to matrix
     xyzo = mean(xyz); % Compute mean center
     xyzro = xyz-repmat(xyzo,ns,1); % Center data
%
% Loop until Y-Axis Changes are Less Than One (1) Degree
%
     while dang>ptol % Loop until angle change is small
%
          if leg==1 % Check for right knee
            if xyz_post(1,1)>xyz_post(2,1) % Compare X coordinates
              xme(1) = xyz_post(1,1); % Assign first point X
              yme(1) = xyz_post(1,2); % Assign first point Y
              zme(1) = xyz_post(1,3); % Assign first point Z
              xme(2) = xyz_post(2,1); % Assign second point X
              yme(2) = xyz_post(2,2); % Assign second point Y
              zme(2) = xyz_post(2,3); % Assign second point Z
            else
              xme(1) = xyz_post(2,1); % Assign second point X
              yme(1) = xyz_post(2,2); % Assign second point Y
              zme(1) = xyz_post(2,3); % Assign second point Z
              xme(2) = xyz_post(1,1); % Assign first point X
              yme(2) = xyz_post(1,2); % Assign first point Y
              zme(2) = xyz_post(1,3); % Assign first point Z
            end
          else % Handle left knee
            if xyz_post(1,1)<xyz_post(2,1) % Compare X coordinates
              xme(1) = xyz_post(1,1); % Assign first point X
              yme(1) = xyz_post(1,2); % Assign first point Y
              zme(1) = xyz_post(1,3); % Assign first point Z
              xme(2) = xyz_post(2,1); % Assign second point X
              yme(2) = xyz_post(2,2); % Assign second point Y
              zme(2) = xyz_post(2,3); % Assign second point Z
            else
              xme(1) = xyz_post(2,1); % Assign second point X
              yme(1) = xyz_post(2,2); % Assign second point Y
              zme(1) = xyz_post(2,3); % Assign second point Z
              xme(2) = xyz_post(1,1); % Assign first point X
              yme(2) = xyz_post(1,2); % Assign first point Y
              zme(2) = xyz_post(1,3); % Assign first point Z
            end
          end
%
% Get Y-Axis
%
          yax = [xme(2)-xme(1) yme(2)-yme(1) zme(2)-zme(1)]'; % Compute Y-axis vector
%
% Get Rotation Angle and Rotate about the Z-Axis
%
          ang(1) = ang(2); % Store previous angle
          ang(2) = -atan(yax(2)/yax(1)); % Compute new rotation angle
          dang = diff(abs(ang)); % Compute angle difference
%
          r = rotxyz([0 0 ang(2)]); % Compute rotation matrix
          xyzr = xyzro*r'; % Rotate centered points
%
% Find Most Posterior Points on the Posterior Tibia Plateau at the PCL
%
          ids = zeros(2,1); % Initialize indices array
          for n = 1:2 % Loop over sides
             [~,id] = max(xyzr(idp{n},2)); % Find max posterior Y
             ids(n) = idp{n}(id); % Store index
             xyz_post(n,:) = xyz(ids(n),:); % Update posterior point
          end
     end
%
% Plot Posterior Points and Line
%
     if iplt % Check if plotting is enabled
       plot3(xyz_post(:,1),xyz_post(:,2),xyz_post(:,3), ...
             'bo','LineWidth',2); % Plot posterior points
     end
   end
%
% Get Proximal Tibia Section
%
   if strcmp(roinaml,'prox_tibia') % Check for 'prox_tibia' ROI
%
     if nslice>1 % Check for excessive slices
       warning([' *** WARNING in tibia_cs8:  Too many PROX_TIBIA ', ...
                'slices in ' faxial '!']); % Issue warning
     end
%
     xyz = dat{1}; % Get first slice data
     npts = size(xyz,1); % Get number of points
%
     xyzm = mean(xyz); % Compute mean of points
     xyzp = [xyz; xyzm]; % Append mean to points
%
     if iplt % Check if plotting is enabled
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-', ... 
             'LineWidth',1); % Plot slice
     end
%
% Get and Plot Centroid of Proximal Tibia Section
%
     tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)']; % Create triangle connectivity
     tri = [tri; [npts+1 npts 1]]; % Add final triangle
%
     [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri); % Compute triangle areas and centroids
     xyzc2 = sum(repmat(at,1,3).*cgt)./sum(at); % Compute weighted centroid
     maxheight = xyzc2(3); % Store maximum height
     if iplt % Check if plotting is enabled
       plot3(xyzc2(1),xyzc2(2),xyzc2(3),'bs', ...
             'LineWidth',1); % Plot centroid
     end
   end
%
end
%
% Draw a Line Connecting Posterior Tibia Plateau Points
%
if iplt % Check if plotting is enabled
  plot3(xme,yme,zme,'b-','LineWidth',2); % Plot line between posterior points
end
%
% Get Coordinate System
%
rvec = xyzc2 - xyzc1; % Compute vertical Z-axis vector
%
% Get Coordinate System Origin
%
xyzm = [xme yme zme]; % Combine posterior points
xyzo = li_clos([xyzc1;xyzc2],xyzm); % Find closest point on centroid line
xc = xyzo(1); % Extract X coordinate
yc = xyzo(2); % Extract Y coordinate
zc = xyzo(3); % Extract Z coordinate
%
% Get Coordinate System
%
yax = [xme(2)-xme(1) yme(2)-yme(1) zme(2)-zme(1)]'; % Compute Y-axis vector
yax = yax./norm(yax); % Normalize Y-axis
zax = sign(rvec(3))*rvec'; % Compute Z-axis vector
zax = zax./norm(zax); % Normalize Z-axis
xax = cross(yax,zax); % Compute X-axis vector
xax = xax./norm(xax); % Normalize X-axis
yax = cross(zax,xax); % Recompute Y-axis vector
yax = yax./norm(yax); % Normalize Y-axis
%
% Plot Coordinate System
%
if iplt % Check if plotting is enabled
  h = quiver3(xc,yc,zc,xax(1),xax(2),xax(3),50,'k'); % Plot X-axis
  set(h,'LineWidth',2); % Set line width
  h = quiver3(xc,yc,zc,yax(1),yax(2),yax(3),50,'k'); % Plot Y-axis
  set(h,'LineWidth',2); % Set line width
  h = quiver3(xc,yc,zc,zax(1),zax(2),zax(3),50,'k'); % Plot Z-axis
  set(h,'LineWidth',2); % Set line width
  axis equal; % Set equal axis scaling
end
%
% Return Axis Vectors in Matrices
%
xyzax = [xax yax zax]; % Combine axis vectors
%
% Calculate Aspect Ratio Based on Maximum Width of the Proximal Tibia
% in the Tibial Coordinate System and Height Between the Distal and
% Proximal Tibia Outlines
%
npts = size(xyzp,1); % Get number of points
xyzp = xyzp - repmat(xyzo,npts,1); % Center points
xyzp = xyzp*xyzax; % Transform to tibia coordinate system
widt = abs(max(xyzp(:,1))-min(xyzp(:,1))); % Compute width
minheight = xyzc1(3); % Set minimum height
maxheight = xyzc2(3); % Set maximum height
height = maxheight-minheight; % Compute height
aspect = widt/height; % Compute aspect ratio
%
% Proximal Tibia Outline in Tibia Coordinate System
%
xyzpt = xyzp(1:end-1,:); % Extract outline points
prox_tib_x = sortrows(xyzpt,1); % Sort by X coordinate
prox_tib_y = sortrows(xyzpt,2); % Sort by Y coordinate
prox_tib_ap = abs(prox_tib_x(1,1) - prox_tib_x(end,1)); % Compute AP dimension
prox_tib_ml = abs(prox_tib_y(1,2) - prox_tib_y(end,2)); % Compute ML dimension
%
if nargout==0 % Check if no output arguments
  fprintf(1,'\n Aspect ratio = %5.3f\n\n'); % Print aspect ratio
end
% %%
% %
% % Proximal Tibial ML and AP dimensions by min/max
% %
% prox_tib = roi(3).data{1, 1}; % Proximal Tibial Data 
% prox_tib_x = sortrows(prox_tib,1); % Proximal tibial Data Sorted for X
% prox_tib_y = sortrows(prox_tib,2); % Proximal Tibial data sorted for y
% prox_tib_ml = abs(prox_tib_x(1,1) - prox_tib_x(end,1)); % ML Dimension
% prox_tib_ap = abs(prox_tib_y(1,2) - prox_tib_y(end,2)); % AP Dimesnion
%%
%
return % Exit the function