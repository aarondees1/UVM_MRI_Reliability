function zg = gridproj(trip,xyzp,xg,yg,zdir,dist) % Define function to project grid onto triangular mesh
%GRIDPROJ Projects a two dimensional grid along a vector perpendicular
%         to the grid onto a triangular mesh and finds the most distant
%         intersection along the vector.
%
%         ZG = GRIDPROJ(TRIP,XYZP,XG,YG,ZDIR) given a three (3) columns
%         triangle connectivity matrix, TRIP, a three (3) columns
%         coordinate point data matrix, XYZP, grid X and Y coordinates
%         in vectors, XG and YG, and a Z direction parameter, ZDIR,
%         projects the grid onto the surface defined by the mesh and
%         returns the most distant surface intersection Z coordinate in
%         the positive Z direction if ZDIR is greater than or equal to
%         zero or in the negative Z direction if ZDIR is less than zero.
%         The most distant surface intersection Z coordinates are
%         returned in vector, ZG.
%
%         ZG = GRIDPROJ(TRIP,XYZP,XG,YG,ZDIR,DIST) Finds mesh points (and
%         connecting triangles) within DIST units of the projected lines.
%         Default distance is 2.4 units.
%
%         NOTES:  1.  Function to find the posterior of the bony
%                 patella.
%
%                 2.  Finds mesh points (and connecting triangles)
%                 within 2.4 units of the projected lines.  See line 63.
%
%                 3.  The M-files nod2tri.m, pts2lin.m, tsect4.m and
%                 xprod.m must be in the current path or directory.
%
%         09-Dec-2015 * Mack Gardner-Morse
%
%         15-Mar-2022 * Mack Gardner-Morse * made distance to projected
%                                            line an input variable.
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<4) % Check if fewer than 4 inputs are provided
  error(' *** ERROR in GRIDPROJ:  Not enough inputs!'); % Throw error for insufficient inputs
end
%
if (nargin<6) % Check if distance input is missing
  dist = 2.4*2.4; % Set default squared distance
else
  dist = dist*dist; % Square provided distance
end
%
if (nargin<5) % Check if direction input is missing
  zdir = 1; % Set default positive direction
end
%
% Check Inputs
%
ncolt = size(trip,2); % Get number of columns in triangle connectivity
ncolp = size(xyzp,2); % Get number of columns in point data
%
if (ncolt~=3)||(ncolp~=3) % Check if inputs have 3 columns
  error([' *** ERROR in GRIDPROJ:  Mesh inputs must have', ...
         ' three (3) columns!']); % Throw error for invalid column count
end
%
xg = xg(:); % Convert X grid to column vector
yg = yg(:); % Convert Y grid to column vector
ngp = size(xg,1); % Get number of X grid points
ngpt = size(yg,1); % Get number of Y grid points
%
if ngp~=ngpt % Check if grid vectors have same length
  error([' *** ERROR in GRIDPROJ:  Grid inputs must have', ...
         ' the same lengths!']); % Throw error for mismatched lengths
end
%
% Get Projection Direction
%
if zdir<0 % Check if direction is negative
  zdir = -1; % Set negative direction
else
  zdir = 1; % Set positive direction
end
%
zvec = [0 0 zdir]; % Define projection vector
%
% Loop Through Grid Points
%
zg = zeros(ngp,1); % Initialize output Z coordinates
%
for k = 1:ngp % Loop through grid points
%
% Find Points Within a Distance of the Grid Projection Line
%
   gpt = [xg(k) yg(k) 0]; % Define grid point
   xyzl = pts2lin(gpt,zvec,xyzp); % Project line through grid point
   d = xyzl-xyzp; % Compute differences
   d = sum(d.*d,2); % Calculate squared distances
   idp = find(d<dist); % Find points within distance
%
% Check for Intersections with Mesh Surface
%
   if ~isempty(idp) % Check if points are found
     idt = nod2tri(idp,trip); % Get triangles for nearby points
     nt = size(idt,1); % Get number of triangles
%
% Loop through Triangles Close to Projection Line Looking for
% Intersections
%
     xyzi = NaN(nt,3); % Initialize intersection points
%
     for l = 1:nt % Loop through triangles
        it = trip(idt(l),:)'; % Get triangle node indices
        v1 = xyzp(it(1),:); % Get first vertex
        v2 = xyzp(it(2),:); % Get second vertex
        v3 = xyzp(it(3),:); % Get third vertex
        [xyzit,il] = tsect4(v1,v2,v3,gpt,zvec); % Check for intersection
        if il % Check if intersection exists
          xyzi(l,:) = xyzit'; % Store intersection point
        end
     end
%
% Find Most Distant Intersection Point
%
     idd = find(~isnan(xyzi(:,3))); % Find valid intersections
     if ~isempty(idd) % Check if intersections exist
       zg(k) = zdir*max(zdir*xyzi(idd,3)); % Compute most distant Z
     else
       zg(k) = NaN; % Set Z to NaN
     end
   else
     zg(k) = NaN; % Set Z to NaN
   end
end
%
return % Exit the function