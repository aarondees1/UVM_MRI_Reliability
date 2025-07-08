function [nx,ny,nz,xc,yc,zc] = plnorm(x,y,z) % Define function to compute triangle normal and centroid
%PLNORM Computes the normal and centroid of a 3-D triangle.
%       [Nx,Ny,Nz,Xc,Yc,Zc] = PLNORM(X,Y,Z) returns the normal and
%       centroid of a triangle defined by three (3) points with
%       coordinates X,Y and Z.  X, Y and Z must be vectors of length
%       three (3).  The coordinates should be listed in clockwise
%       direction when seen from above.
%
%       28-Mar-96
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<3) % Check if fewer than 3 inputs are provided
  error('PLNORM requires three input arguments.'); % Throw error for insufficient inputs
end
%
% Check that Inputs are Vectors
%
x = x(:); % Convert X coordinates to column vector
y = y(:); % Convert Y coordinates to column vector
z = z(:); % Convert Z coordinates to column vector
%
[n l] = size(x); % Get dimensions of X coordinates
[n2 l2] = size(y); % Get dimensions of Y coordinates
[n3 l3] = size(z); % Get dimensions of Z coordinates
%
if ((l~=1)|(l2~=1)|(l3~=1)) % Check if inputs are vectors
  error('PLNORM only works with vectors.'); % Throw error for non-vector inputs
end
%
% Check that the Inputs have Three Rows
%
if ((n~=3)|(n2~=3)|(n3~=3)) % Check if inputs have length 3
  error('X, Y and Z must be of length three (3).'); % Throw error for invalid length
end
%
% Normal of Triangle
%
a = [(x(1)-x(2)); (y(1)-y(2)); (z(1)-z(2))]; % Compute vector from point 2 to 1
b = [(x(3)-x(2)); (y(3)-y(2)); (z(3)-z(2))]; % Compute vector from point 2 to 3
n = cross(a,b); % Compute cross product for normal
n = n/norm(n); % Normalize normal vector
%
nx = n(1); % Extract X component of normal
ny = n(2); % Extract Y component of normal
nz = n(3); % Extract Z component of normal
%
% Centroid of Triangle
%
xc = mean(x); % Compute X coordinate of centroid
yc = mean(y); % Compute Y coordinate of centroid
zc = mean(z); % Compute Z coordinate of centroid
%
return % Exit the function