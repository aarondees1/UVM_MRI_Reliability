function dist = dist2cyl(x,xyzp)
%DIST2CYL Finds the distances from a set of points and the surface of a
%         cylinder.
%
%         DIST = DIST2CYL(X,XYZP) given the seven (7) parameter vector
%         X that describes the cylinder and a matrix of points XYZP
%         with the X, Y and and Z coordinates in columns, returns the
%         distances between the points and the surface of a cylinder.
%
%         The seven (7) parameters in X are:
%         1.  radius of the cylinder
%         2.  X coordinate of first point on the axis of the cylinder
%         3.  Y coordinate of first point on the axis of the cylinder
%         4.  Z coordinate of first point on the axis of the cylinder
%         5.  X coordinate of second point on the axis of the cylinder
%         6.  Y coordinate of second point on the axis of the cylinder
%         7.  Z coordinate of second point on the axis of the cylinder
%
%         NOTES:  1.  For use with cyl_fit.m.  See cyl_fit.m for more
%                 details.
%
%                 2.  The M-file pts2lin.m must be in the current path
%                 or directory.
%
%         17-Jul-2013 * Mack Gardner-Morse
%

% Check if two input arguments are provided.
if (nargin<2)
  error(' *** ERROR in DIST2CYL:  Two (2) inputs are required!');
end

% Validate cylinder parameter vector: must have 7 elements.
x = x(:); % Ensure x is a column vector
nparam = size(x,1);
if nparam~=7
  error([' *** ERROR in DIST2CYL:  Seven parameters are required', ...
         ' to define the cylinder!']);
end

% Validate point coordinate matrix: must have 3 columns.
[np,ncol] = size(xyzp);
if ncol~=3
  error([' *** ERROR in DIST2CYL:  Points coordinate matrix', ...
         ' must have three columns!']);
end

% Extract cylinder parameters from input vector.
r = x(1); % Cylinder radius
lpt0 = x(2:4)'; % First point on cylinder axis
lvec = x(5:7)'-lpt0; % Direction vector from first to second point
lvec = lvec./sqrt(lvec*lvec'); % Normalize to unit vector

% Compute closest points on cylinder axis to input points.
xyzc = pts2lin(lpt0,lvec,xyzp); % Coordinates of closest points on axis

% Calculate distances from input points to cylinder surface.
dist = xyzc-xyzp; % Vectors from input points to closest axis points
dist = sqrt(sum(dist.*dist,2)); % Euclidean distances to cylinder axis
dist = dist-r; % Subtract radius to get distances to cylinder surface

return % Exit function