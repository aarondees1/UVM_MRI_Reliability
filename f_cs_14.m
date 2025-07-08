function [xyzo,xyzax,xyzlnx,r,s,sc,sy,ssep,hf] = f_cs_14(datlb,datmb,datfs,datcp,iplt,hf) % Define function to calculate femur coordinate system
%F_CS_14   Uses MRI femur (F) data cell arrays to calculate a femur
%          based coordinate system (CS).
%
%          [XYZO,XYZAX] = F_CS_14(DATLB,DATMB,DATFS,DATCP)
%          given the MRI femur data cell arrays for the lateral condyle,
%          DATLB, the medial condyle, DATMB, the femoral shaft, DATFS
%          and the center point of the notch, DATCP, calculates a femur
%          coordinate system based on a cylinder fit of the condyle bone
%          surfaces, the centroid of the femur shaft from just above
%          the patella and the center point of the notch.  The function
%          returns the origin of the femur coordinate system, in XYZO,
%          and unit X, Y and Z vectors for the X, Y, and Z axes in the
%          rows of matrix XYZAX (rotation matrix).
%
%          [XYZO,XYZAX,R] = F_CS_14(DATLB,DATMB,DATFS,DATCP) returns
%          the radius of the fitted cylinder, R.
%
%          NOTES:  1.  The Matlab files cyl_fit.m, cyl_plt.m,
%                  dist2cyl.m, plt_datsl.m, pts2lin.m, pt2line.m, and
%                  tri_area.m must be in the current directory or path.
%
%                  2.  Note the femur coordinate system is defined as
%                  the positive X-axis as either anterior (left knees)
%                  or posterior (right knees), positive Y-axis as
%                  lateral and positive Z-axis as superior.
%
%                  3.  Based on the femur coordinate system work of 
%                  Daniel Sturnick.
%
%          29-Jul-2014 * Mack Gardner-Morse
%
%          29-Nov-2022 * Mack Gardner-Morse * Uses the new cyl_fit.m
%          function which uses the Matlab nonlinear least squares
%          function lsqnonlin for fitting the cylinder.
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<4) % Check if fewer than 4 inputs are provided
  error(' *** ERROR in F_CS_14:  Four input variables are required!'); % Throw error for insufficient inputs
end
%
if (nargin<5)||isempty(iplt) % Check if plot flag is missing or empty
  iplt = false; % Set default plot flag to false
end
%
if (nargin<6) % Check if figure handle is provided
  if iplt % Check if plotting is enabled
    hf = figure; % Create new figure
    orient landscape; % Set figure orientation to landscape
    hold on; % Enable hold for multiple plots
    view(3); % Set 3D view
  end
else
  figure(hf); % Activate specified figure
  orient landscape; % Set figure orientation to landscape
  hold on; % Enable hold for multiple plots
  view(3); % Set 3D view
end
%
% Read Femoral Shaft Data and Get Scaling
%
xyzfs = datfs{1}; % Get femoral shaft data from first cell
nptfs = size(xyzfs,1); % Get number of points in femoral shaft data
%
xyzm = mean(xyzfs); % Compute mean of femoral shaft points
xyzp = [xyzfs; xyzm]; % Append mean point to femoral shaft data
%
if iplt % Check if plotting is enabled
  plot3(xyzfs(:,1),xyzfs(:,2),xyzfs(:,3),'k.-','MarkerSize',6,'LineWidth',1); % Plot femoral shaft points
end
%
s = max(xyzfs)-min(xyzfs); % Compute scaling in MRI coordinate system
%
% Get and Plot Centroid
%
tri = [repmat(nptfs+1,nptfs-1,1) (1:nptfs-1)' (2:nptfs)']; % Define triangles for centroid calculation
tri = [tri; [nptfs+1 nptfs 1]]; % Add final triangle
%
[at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri); % Compute triangle areas and centroids
xyzc = sum(repmat(at,1,3).*cgt)./sum(at); % Compute weighted centroid
%
if iplt % Check if plotting is enabled
  plot3(xyzc(1),xyzc(2),xyzc(3),'bs','LineWidth',2,'MarkerSize',6); % Plot centroid
end
%
% Get Notch Coordinates
%
xyz0 = datcp{1}; % Get notch center point
if iplt % Check if plotting is enabled
  plot3(xyz0(1),xyz0(2),xyz0(3),'gs','Color',[0 0.5 0],'LineWidth',2,'MarkerSize',6); % Plot notch point
end
%
% Plot a Line from Notch to Centroid of Slice 1 (Z-Axis Direction)
%
xyzl2 = [xyz0; xyzc]; % Define line from notch to centroid
if iplt % Check if plotting is enabled
  plot3(xyzl2(:,1),xyzl2(:,2),xyzl2(:,3),'g-','Color',[0 0.5 0],'LineWidth',2); % Plot line
end
%
% Add Labels to Figure
%
if iplt % Check if plotting is enabled
  xlabel('X (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
%   title(title_txt,'FontSize',16,'FontWeight','bold','interpreter','none');
  axis equal; % Set equal axis scaling
end   
%
% Get Condyle Sagittal Data
%
xyzlat = cell2mat(datlb); % Convert lateral condyle data to matrix
xyzmed = cell2mat(datmb); % Convert medial condyle data to matrix
%
% [nslice,nsl,isl] = sl_info(datlb);
%
if iplt % Check if plotting is enabled
  plt_datsl(datlb,'b.-',0.5,[0 0 0.7]); % Plot lateral condyle data
  plt_datsl(datmb,'g.-',0.5,[0 0.5 0]); % Plot medial condyle data
end
%
% Fit Cylinder to Slice Data
%
ri = 15; % Set initial cylinder radius
%ri = 25;                % For Knee ID #'s 8,30,62,71
%ri = 32                 % For Knee ID #'s 29,40  
xyz1i = mean(xyzlat); % Compute mean of lateral condyle points
xyz2i = mean(xyzmed); % Compute mean of medial condyle points
[r,xyz1,xyz2,ssep] = cyl_fit(ri,xyz1i,xyz2i,[xyzlat; xyzmed]); % Fit cylinder to condyle data
%
xyzlnx = [xyz1; xyz2]; % Define cylinder axis line
%
if iplt % Check if plotting is enabled
  hc = cyl_plt(r,xyz1,xyz2); % Plot fitted cylinder
  set(hc,'EdgeColor',[0.5 0.5 0.5]); % Set cylinder edge color
%
  plot3(xyzlnx(:,1),xyzlnx(:,2),xyzlnx(:,3),'ko-','LineWidth',2,'Color',[0.5 0.5 0.5]); % Plot cylinder axis
end
%
% Get Coordinate System
%
yax = -diff(xyzlnx); % Compute Y-axis direction
sy = norm(yax); % Compute Y-axis length
yax = yax./sy; % Normalize Y-axis
zax = diff(xyzl2); % Compute Z-axis direction
zax = zax./norm(zax); % Normalize Z-axis
% zax = sign(zax(3))*rvec;
xax = cross(yax,zax); % Compute X-axis via cross product
xax = xax./norm(xax); % Normalize X-axis
zax = cross(xax,yax); % Recompute Z-axis for orthogonality
zax = zax./norm(zax); % Normalize Z-axis
%
% Get Axis Origin on Cylinder Axis Line Nearest Notch Point
% (pt2line.m)
%
xyzo = pt2line(xyz1,yax,xyz0); % Compute origin on cylinder axis
%
xyzax = [xax; yax; zax]'; % Form rotation matrix
%
% Plot Coordinate System Origins
%
if iplt % Check if plotting is enabled
  xc = repmat(xyzo(1),1,3); % Replicate X-coordinate for plotting
  yc = repmat(xyzo(2),1,3); % Replicate Y-coordinate for plotting
  zc = repmat(xyzo(3),1,3); % Replicate Z-coordinate for plotting
%
  plot3(xc(1),yc(1),zc(1),'ro','MarkerSize',9,'LineWidth',2); % Plot origin
  h = quiver3(xc,yc,zc,xyzax(1,:),xyzax(2,:),xyzax(3,:),50,'c'); % Plot coordinate axes
  set(h,'LineWidth',3); % Set line width for axes
%
  axis equal; % Set equal axis scaling
end
%
% Femoral Shaft Scaling
%
xyzfst = xyzfs-repmat(xyzo,nptfs,1); % Transform femoral shaft to femur CS origin
xyzfst = xyzfst*xyzax'; % Rotate to femur CS
%
sc = max(xyzfst)-min(xyzfst); % Compute scaling in femur CS
%
return % Exit the function