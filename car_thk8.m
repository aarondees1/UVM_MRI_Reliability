function [t,ip] = car_thk8(trib,xyzb,tric,xyzc,rad,iplt,fstr,cmprt,hf1); % Define function with input and output arguments
%CAR_THK8 Finds cartilage thicknesses at all the points on the
%         interpolated grid cartilage surface to the nearest bone
%         surface.
%
%         T = CAR_THK8(TRIB,XYZB,TRIC,XYZC) given the subchondral
%         bone triangle connectivity and coordinates matrices TRIB and
%         XYZB, and cartilage triangle connectivity and coordinates
%         matrices TRIC and XYZC, determines the cartilage thicknesses
%         at all the points on the cartilage surface.  The thicknesses
%         are based on the cartilage surface normals.
%
%         T = CAR_THK8(TRIB,XYZB,TRIC,XYZC,RAD) given the radius RAD,
%         checks for nodes and connected triangles within that radius
%         for intersections.  By default, RAD = 12.
%
%         [T,IP] = CAR_THK8(TRIB,...) returns the intersection points
%         with the subchondral bone surface IP in a three column matrix
%         for the X, Y and Z coordinates of the points.
%
%         CAR_THK8(TRIB,XYZB,TRIC,XYZC,RAD,IPLT,FSTR,CMPRT) if IPLT is
%         true and given a title string FSTR and a logical compartment
%         flag (false for lateral, true for medial) CMPRT, plots the
%         cartilage surface color coded by cartilage thicknesses in a
%         figure.
%
%         CAR_THK8(TRIB,XYZB,TRIC,XYZC,RAD,IPLT,FSTR,CMPRT,HF1) plots
%         the cartilage surface color coded by cartilage thicknesses in
%         the figure with figure handle HF1.
%
%         NOTES:  1. Searches for the nearest intersection in the
%                 normal surface direction which may be slower than
%                 just finding an intersection (see CAR_THK6.M).
%
%                 2.  The M-files nod2tri.m, nod_norm.m, tsect4.m and
%                 xprod.m must be in the current path or directory.
%
%                 3.  tsect4.m uses a very strict definition of
%                 intersection.  An alternate is tsect3.m which finds
%                 intersections to within a tolerance.
%
%         04-Jan-2016 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<4) % Check if fewer than 4 inputs are provided
  error([' *** ERROR in CAR_THK8:  CAR_THK8 requires four (4)', ...
         ' inputs!']); % Throw error for insufficient inputs
end
%
if (nargin<5)||isempty(rad) % Check if radius input is missing or empty
  rad = 12; % Set default radius to 12
end
%
if (nargin<6)||isempty(iplt) % Check if plot flag is missing or empty
  iplt = false; % Set default plot flag to false
end
%
if (nargin<7)||isempty(fstr) % Check if title string is missing or empty
  fstr = ''; % Set default title string to empty
end
%
% Compartment Names
%
if (nargin<8)||isempty(cmprt) % Check if compartment flag is missing or empty
  cmprt = false; % Set default compartment to false (lateral)
  txtcmpt = ''; % Set default compartment text to empty
else
  if cmprt % Check if compartment is medial
    txtcmpt = 'Medial'; % Set compartment text to 'Medial'
  else
    txtcmpt = 'Lateral'; % Set compartment text to 'Lateral'
  end
end
%
% Open Figure
%
if iplt % Check if plotting is enabled
  if nargin<9 % Check if figure handle is provided
    hf1 = figure; % Create new figure
  else
    figure(hf1); % Activate specified figure
  end
  orient landscape; % Set figure orientation to landscape
end
%
% Get State of Warnings
%
swrn = warning('query','backtrace'); % Query current warning state for backtrace
warning('off','all'); % Turn off all warnings
%
% Plot Cartilage and Bone Surfaces
%
if iplt % Check if plotting is enabled
  figure(hf1); % Activate figure
  subplot(1,2,1); % Create first subplot
  hc1 = trimesh(tric,xyzc(:,1),xyzc(:,2),xyzc(:,3), ...
                'LineWidth',0.25,'FaceColor','none','EdgeColor','b'); % Plot cartilage mesh
  hold on; % Enable hold for multiple plots
  xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
  ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
  zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
  title({[fstr ' ' txtcmpt]; ...
        'Cartilage and Bone Meshes'},'Interpreter','none', ...
        'FontSize',16,'FontWeight','bold'); % Add title with file and compartment
%
% Plot Bone Surface
%
  hb = trimesh(trib,xyzb(:,1),xyzb(:,2),xyzb(:,3), ...
               'LineWidth',0.25,'FaceColor','none','EdgeColor','k'); % Plot bone mesh
  if cmprt % Check if compartment is medial
    view(-40,60); % Set view angle for medial compartment
  end
  axis equal; % Set equal axis scaling
end
%
% Get Cartilage Surface Normals at the Nodes
%
nnv = nod_norm(tric,xyzc); % Compute cartilage surface normals
%
% Loop through the Cartilage Points
%
nb = size(xyzb,1); % Get number of bone nodes
nc = size(xyzc,1); % Get number of cartilage nodes
xyzi = NaN(nc,3); % Initialize intersection points array
for k = 1:nc % Loop through each cartilage point
   pt = xyzc(k,:); % Get current cartilage point
   dmin = Inf; % Initialize minimum distance to infinity
%
   if ~isnan(pt(:,3)) % Check if z-coordinate is valid
     linv = nnv(k,:); % Get cartilage normal vector
%
% Find Nodes and Connected Triangles Close (<RAD mm) to the Cartilage Points
%
     r = xyzb-repmat(pt,nb,1); % Compute distance vectors from cartilage point to bone nodes
     r = sum(r.*r,2); % Compute squared distances
     idb = find(r<rad*rad); % Identify bone nodes within radius
     idt = nod2tri(idb,trib); % Get triangles connected to these bone nodes
     ntri = size(idt,1); % Get number of triangles
%
% Find the Bone Intersection with Cartilage Normal
%
     for l = 1:ntri % Loop through bone surface triangles
        it = trib(idt(l),:)'; % Get nodes in current triangle
        v1 = xyzb(it(1),:); % Get first triangle vertex
        v2 = xyzb(it(2),:); % Get second triangle vertex
        v3 = xyzb(it(3),:); % Get third triangle vertex
        [xyzit,il,ierr] = tsect4(v1,v2,v3,pt,linv); % Find intersection with normal
        if il % Check if intersection exists
          adis = xyzit-pt'; % Compute distance vector to intersection
          adis = adis'*adis; % Compute squared distance
          if adis<dmin % Check if distance is smallest
            dmin = adis; % Update minimum distance
            xyzi(k,:) = xyzit'; % Store intersection point
          end
        end
        if l==ntri&&dmin==Inf % Check if no intersection found
          warning(sprintf(['  CAR_THK8 is not able to ', ...
                  'find subchondral bone below and normal to ', ...
                  'cartilage point %i in %s compartment in %s.'], ...
                  k,lower(txtcmpt),fstr)); % Issue warning
        end
     end
   end
end
%
% Get Cartilage Thicknesses
%
dp = xyzc-xyzi; % Compute vectors from cartilage to intersection points
cthk = sqrt(sum(dp.*dp,2)); % Compute cartilage thicknesses
%
% Plot the Points for the Cartilage Thicknesses
%
if iplt % Check if plotting is enabled
  figure(hf1); % Activate figure
  subplot(1,2,2); % Create second subplot
  hc2 = trimesh(tric,xyzc(:,1),xyzc(:,2),xyzc(:,3), ...
                'LineWidth',0.25,'FaceColor','none','EdgeColor','b'); % Plot cartilage mesh
  hold on; % Enable hold for multiple plots
  idx = find(~isnan(cthk)); % Identify valid thickness points
  it = nod2tri(idx,tric,2); % Get triangles connected to valid points
  hs = trisurf(tric(it,:),xyzc(:,1),xyzc(:,2), ...
               xyzc(:,3),cthk,'FaceColor','interp', ...
               'EdgeColor','b','LineWidth',0.25); % Plot thickness surface
  if cmprt % Check if compartment is medial
    view(-40,60); % Set view angle for medial compartment
  end
  axis equal; % Set equal axis scaling
  if min(cthk)<0 % Check if minimum thickness is negative
    caxis([min(cthk) max(cthk)]); % Set color axis to thickness range
  else
    caxis([0 max(cthk)]); % Set color axis from 0 to max thickness
  end
  hcb = colorbar; % Add colorbar
  xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
  ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
  zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
  title({[fstr ' ' txtcmpt]; ...
        'Cartilage Thicknesses'},'Interpreter','none', ...
        'FontSize',16,'FontWeight','bold'); % Add title with file and compartment
end
%
% Return Thicknesses and Intersection Points
%
t = cthk; % Assign thicknesses to output
ip = xyzi; % Assign intersection points to output
%
% Restore Warnings
%
warning(swrn); % Restore original warning state
%
return % Exit the function