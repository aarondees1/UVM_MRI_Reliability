function hc = cyl_plt(r,xyz1,xyz2,n,fc)
%CYL_PLT  Plots a cylinder given a radius and the centers of the two
%         ends of a cylinder.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2) given the radius of the cylinder R,
%         a row vector of the coordinates of the center of one end of
%         the cylinder XYZ1, and a row vector of the coordinates of the
%         center of the other end of the cylinder XYZ2, plots the
%         cylinder as a lighted surface using SURFL and returns the
%         graphics handle HC.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2,N) sets the number of circumferential
%         facets to N.  The default is 72.
%
%         HC = CYL_PLT(R,XYZ1,XYZ2,N,FC) sets the face color of the
%         cylinder to FC.  See "set(hc,'FaceColor')" for valid color
%         options.  By default there is no face color (only a wire
%         frame cylinder).
%
%         NOTES:  1.  For use with cyl_fit.m to plot the resulting
%                 cylinder.
%
%         17-Jul-2013 * Mack Gardner-Morse
%

% Check if at least three input arguments are provided.
if (nargin<3)
  error(' *** ERROR in CYL_PLT:  Three (3) inputs are required!');
end

% Set default number of circumferential facets if not provided or empty.
if nargin<4||isempty(n)
  n = 72; % Default number of facets
end

% Set default face color if not provided or empty.
if nargin<5||isempty(fc)
  fc = 'none'; % Default to wireframe (no face color)
end

% Validate input dimensions: xyz1 and xyz2 must be 1x3 row vectors.
[nrow1,ncol1] = size(xyz1);
[nrow2,ncol2] = size(xyz2);
if nrow1~=1|nrow2~=1|ncol1~=3|ncol2~=3
  error([' *** ERROR in CYL_PLT:  All input points must have', ...
         ' three columns!']);
end

% Calculate cylinder parameters.
zc = xyz2-xyz1; % Vector from xyz1 to xyz2 (cylinder axis direction)
l = zc*zc'; % Squared length of cylinder axis
l = sqrt(l); % Length of the cylinder
zc = zc./l; % Normalize to unit direction vector
m = (xyz2+xyz1)/2; % Center point of the cylinder

% Generate cylinder surface with specified radius and facets.
[x,y,z] = cylinder([r r r],n); % Create unit cylinder with radius r and n facets
z = l*z-l/2; % Scale and shift z to match cylinder length and center
dim = size(x); % Dimensions of cylinder surface arrays

% Compute rotation matrix to align cylinder with axis.
yc = zeros(1,3); % Initialize arbitrary Y axis
[val,idm] = min(abs(zc)); % Find smallest component of zc
yc(idm) = 1; % Set Y axis along smallest zc component
xc = cross(yc,zc); % X axis: perpendicular to yc and zc
yc = cross(zc,xc); % Y axis: perpendicular to zc and xc
rot = [xc; yc; zc]; % Rotation matrix from local to global coordinates

% Transform cylinder coordinates to align with axis and position.
a = [x(:)'; y(:)'; z(:)']; % Flatten and stack cylinder coordinates
a = rot'*a; % Rotate to align with cylinder axis
a = a+repmat(m',1,prod(dim)); % Translate to cylinder center
x = reshape(a(1,:),dim); % Reshape X coordinates
y = reshape(a(2,:),dim); % Reshape Y coordinates
z = reshape(a(3,:),dim); % Reshape Z coordinates

% Plot cylinder as a lighted surface.
hc = surfl(x,y,z); % Create surface plot with lighting
set(hc,'FaceColor',fc); % Set face color (or none for wireframe)
set(hc,'LineWidth',0.5); % Set line width for edges

return % Exit function