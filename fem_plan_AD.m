function [lmp,lmv,ppc,ppn,ppl,pnl,ppm,pnm] = fem_plan_AD(f3,tmin,iplt,ltxt,fsav) % Define function to calculate femur dividing planes
%FEM_PLAN  Uses femur bone sagittal plane segmentations to calculate
%          planes dividing the femur into regions of interest (ROIs).
%
%          [LMP,LMV,PPC,PPN,PPL,PNL,PPM,PNM] = FEM_PLAN(F3) Given a
%          sagittal knee segmentation cell array, F3, calculates the
%          lateral-medial dividing plane defined by a point, LMP,
%          and a normal vector, LMV, the posterior dividing plane
%          defined by a point, PPC, and a normal vector, PPN, the
%          lateral plane dividing the lateral condyle from the trochlea
%          defined by a point, PPL, and a normal vector, PNL, and a
%          medial plane dividing the medial condyle from the trochlea
%          defined by a point, PPM, and a normal vector, PNM.
%
%          [...] = FEM_PLAN(F3,IPLT,LTXT,FSAV) given a logical true,
%          IPLT, creates plots of the femoral dividing planes. LTXT is
%          a string with the name of the leg to be used in the plot
%          title.  If the file name, FSAV, is provided, the plots are
%          saved to the file, FSAV.
%
%          NOTES:  1.  The femur bone sagittal plane segmentations
%                  should be in the femoral coordinate system.
%
%                  2.  The Matlab files plnorm.m and plt_datsl.m must
%                  be in the current directory or path.
%
%                  3.  The origin, [0 0 0], is a point in the lateral-
%                  medial plane and the lateral and medial trochlear
%                  planes.
%
%          16-Dec-2022 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<2) % Check if fewer than 2 inputs are provided
  error([' *** ERROR in FEM_PLAN:  An input structure with', ...
         ' sagittal femoral segmentations is required!']); % Throw error for insufficient inputs
end
%
if (nargin<3)||isempty(iplt) % Check if plot flag is missing or empty
  iplt = false; % Set default plot flag to false
elseif iplt == true % Check if plot flag is explicitly true
  iplt = true; % Ensure plot flag is true
end
%
if (nargin<4)||isempty(ltxt) % Check if leg text is missing or empty
  ltxt = []; % Set default leg text to empty
end
%
if (nargin<5)||isempty(fsav) % Check if save file name is missing or empty
  iprt = false; % Set default print flag to false
else
  iprt = true; % Set print flag to true
end
%
% Posterior and Medial/Lateral Cutoffs
%
% tmin = -145/rad2deg;    % Angle (theta) cutoff (-150, -145, or -140)
tmin = deg2rad(tmin); % Convert angle cutoff to radians
y0 = 0; % Set Y cutoff
hpi = pi/2; % Set half pi constant
%
% Lateral-Medial Plane Normal Vector
%
lmp = [0 0 0]; % Define lateral-medial plane point
lmv = [0 1 0]; % Define lateral-medial plane normal vector
%
if iplt % Check if plotting is enabled
  hf2 = figure; % Create new figure
  hold on; % Enable hold for multiple plots
  plt_datsl(f3,'b.-',0.5,[0 0 0.7]); % Plot femoral segmentation data
  axis equal; % Set equal axis scaling
  view(3); % Set 3D view
  grid on; % Enable grid
  xlabel('X (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
  title({'Dividing Planes in Femur CS'; ltxt},'FontSize',16,'FontWeight','bold'); % Set plot title
  view(40,20); % Set view angle
  axis equal; % Re-set equal axis scaling
  %
  axlim = axis; % Get axis limits
  pxyz = [axlim(1) 0 axlim(5); axlim(2) 0 axlim(5); axlim(2) 0 axlim(6); axlim(1) 0 axlim(6); axlim(1) 0 axlim(5)]; % Define lateral-medial plane vertices
  patch(pxyz(:,1),pxyz(:,2),pxyz(:,3),pxyz(:,3),'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33); % Plot lateral-medial plane
end
%
% Convert XZ to Polar Coordinates
%
xyzlbf = cell2mat(f3); % Convert femoral data to matrix
%
[thlf,rlf,zlf] = cart2pol(xyzlbf(:,1),xyzlbf(:,3),xyzlbf(:,2)); % Convert to polar coordinates
idc = thlf>hpi; % Identify angles greater than pi/2
thlf(idc) = thlf(idc)-2*pi; % Adjust angles to range [-3pi/2, pi/2]
%
if iplt % Check if plotting is enabled
  figure; % Create new figure
  plot3(thlf,zlf,rlf,'.'); % Plot theta-Z-R data
  axlim = axis; % Get axis limits
  view(2); % Set 2D view
  hold on; % Enable hold for additional plots
  plot([tmin tmin],axlim(3:4),'r-','LineWidth',1); % Plot angle cutoff line
  xlabel('Theta (radians)','FontSize',12,'FontWeight','bold'); % Label X-axis
  ylabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
  title({'Theta-Z CS'; ltxt},'FontSize',16,'FontWeight','bold'); % Set plot title
else
  tzrmn = min([thlf,rlf,zlf]); % Compute minimum of theta, R, Z
  tzrmx = max([thlf,rlf,zlf]); % Compute maximum of theta, R, Z
  axlim = [tzrmn; tzrmx]; % Combine min/max into axis limits
  axlim = axlim(:)'; % Reshape to row vector
end
%
dr = axlim(6)-axlim(5); % Compute R range
rmn = mean(axlim(5:6)); % Compute mean R
axlim(5) = rmn-0.75*dr; % Adjust minimum R limit
axlim(6) = rmn+1.25*dr; % Adjust maximum R limit
%
ppxyz = [tmin axlim(3) axlim(5); tmin axlim(4) axlim(5); tmin axlim(4) axlim(6); tmin axlim(3) axlim(6); tmin axlim(3) axlim(5)]; % Define posterior plane vertices
[ppxyz(:,1),ppxyz(:,3),ppxyz(:,2)] = pol2cart(ppxyz(:,1),ppxyz(:,3),ppxyz(:,2)); % Convert posterior plane to Cartesian coordinates
%
ppc = mean(ppxyz(1:4,:)); % Compute center of posterior plane
[nx,ny,nz] = plnorm(ppxyz(1:3,1),ppxyz(1:3,2),ppxyz(1:3,3)); % Compute posterior plane normal
ppn = [nx,ny,nz]; % Store posterior plane normal
%
if (iplt) % Check if plotting is enabled
  figure(hf2); % Activate main figure
  patch(ppxyz(:,1),ppxyz(:,2),ppxyz(:,3),ppxyz(:,3),'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33); % Plot posterior plane
end
%
% Get Maximum "Peaks" of Femur in X-Direction
%
nsl = size(f3,1); % Get number of slices
%
mx = zeros(nsl,1); % Initialize X maxima array
my = zeros(nsl,1); % Initialize Y values at X maxima
%
for ks = 1:nsl % Loop through slices
   [mx(ks),idx] = max(f3{ks}(:,1)); % Find maximum X and its index
   my(ks) = f3{ks}(idx,2); % Get corresponding Y value
end
%
[~,ids] = sort(mx); % Sort slices by X maxima
idp = find(my(ids)>0); % Identify lateral slices (Y > 0)
idp = ids(idp(end-2:end)); % Select top three lateral slices
idn = find(my(ids)<0); % Identify medial slices (Y < 0)
idn = ids(idn(end-2:end)); % Select top three medial slices
mxp = mean(mx(idp)); % Compute mean lateral peak X
myp = mean(my(idp)); % Compute mean lateral peak Y
mxn = mean(mx(idn)); % Compute mean medial peak X
myn = mean(my(idn)); % Compute mean medial peak Y
%
% Get Trochlea Groove
%
idxmn = [round(mean(idp)); round(mean(idn))]; % Compute mean indices for lateral/medial peaks
idxmn = sort(idxmn); % Sort indices
idmn = find(ids>idxmn(1)&ids<idxmn(2)); % Identify trochlea groove slices
idmn = ids(idmn(1:3)); % Select first three trochlea slices
mxm = mean(mx(idmn)); % Compute mean trochlea X
mym = mean(my(idmn)); % Compute mean trochlea Y
%
% Get Line Directions
%
xyc = [mxm mym]; % Define trochlea groove center
xyl = [mxp myp]; % Define lateral peak
xym = [mxn myn]; % Define medial peak
vm = xym-xyc; % Compute medial direction vector
vl = xyl-xyc; % Compute lateral direction vector
%
% Get Center of Trochlea (Center of Coordinate System)
%
% xyc = [0 0 0];          % Use origin
%
% Get Lateral and Medial Trochlear Planes
%
vl = [vl 0]; % Extend lateral vector to 3D
vl = vl./norm(vl); % Normalize lateral vector
vm = [vm 0]; % Extend medial vector to 3D
vm = vm./norm(vm); % Normalize medial vector
xv = [1 0 0]; % Define X-axis vector
%
zvl = cross(xv,vl); % Compute lateral intermediate vector
pnl = cross(vl,zvl); % Compute lateral plane normal
pnl = pnl./norm(pnl); % Normalize lateral plane normal
ppl = [0 0 0]; % Define lateral plane point
%
zvm = cross(xv,vm); % Compute medial intermediate vector
pnm = cross(vm,zvm); % Compute medial plane normal
pnm = pnm./norm(pnm); % Normalize medial plane normal
ppm = [0 0 0]; % Define medial plane point
%
% Plot Lateral and Medial Trochlear Planes
%
if iplt % Check if plotting is enabled
%
  figure(hf2); % Activate main figure
  orient landscape; % Set figure orientation to landscape
  axis equal; % Set equal axis scaling
  axlim = axis; % Get axis limits
%
  xlmx = -pnl(2)*axlim(4)/pnl(1); % Compute X limit for lateral plane
  xmmx = -pnm(2)*axlim(3)/pnm(1); % Compute X limit for medial plane
%
  pxyztll = [0 0 axlim(6); 0 0 axlim(5); xlmx axlim(4) axlim(5); xlmx axlim(4) axlim(6); 0 0 axlim(6)]; % Define lateral trochlear plane vertices
%
  pxyztlm = [0 0 axlim(6); xmmx axlim(3) axlim(6); xmmx axlim(3) axlim(5); 0 0 axlim(5); 0 0 axlim(6)]; % Define medial trochlear plane vertices
%
  patch(pxyztll(:,1),pxyztll(:,2),pxyztll(:,3),pxyztll(:,3),'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33); % Plot lateral trochlear plane
  patch(pxyztlm(:,1),pxyztlm(:,2),pxyztlm(:,3),pxyztlm(:,3),'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','FaceAlpha',0.33); % Plot medial trochlear plane
%
  if iprt % Check if printing is enabled
    print('-dpsc2','-r600','-fillpage',fsav); % Save first view to file
    view(-90,90); % Change to top view
    print('-dpsc2','-r600','-fillpage','-append',fsav); % Append top view to file
  end
%
end
%
return % Exit the function