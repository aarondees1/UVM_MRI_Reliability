function [area,cg] = tri_area(x,y,z,tric)
%TRI_AREA Finds the area of 3-D triangles.
%
%         AREA = TRI_AREA(X,Y,Z,TRIC) returns the areas of the
%         triangles defined by the X, Y and Z coordinates and the
%         three (3) column triangle connectivity matrix TRIC.
%
%         [AREA,CG] = TRI_AREA(X,Y,Z,TRIC) returns the area centers of
%         the triangles defined by the X, Y and Z coordinates and the
%         three (3) column triangle connectivity matrix TRIC in the
%         three (3) column matrix CG with the X, Y and Z coordinates of
%         the area centers in the columns
%
%         NOTES:  None.
%
%         13-Dec-04 * Mack Gardner-Morse
%
%         01-Jul-2010 * Mack Gardner-Morse * Include area centers.
%
%         17-Apr-2014 * Mack Gardner-Morse * Made coordinates row
%         vectors so that the function works with just one (1) input
%         triangle.
%

% Check if four input arguments are provided.
if (nargin<4)
  error(' *** ERROR in TRI_AREA:  Not enough input arguments.');
end

% Convert coordinate inputs to row vectors.
x = x(:)'; % X coordinates as row vector
y = y(:)'; % Y coordinates as row vector
z = z(:)'; % Z coordinates as row vector

% Index triangle vertices using connectivity matrix tric.
xe = x(tric); % X coordinates of triangle vertices
ye = y(tric); % Y coordinates of triangle vertices
ze = z(tric); % Z coordinates of triangle vertices

% Calculate triangle edges for cross product.
v1 = [xe(:,2)-xe(:,1) ye(:,2)-ye(:,1) ze(:,2)-ze(:,1)]'; % Edge from vertex 1 to 2
v2 = [xe(:,3)-xe(:,1) ye(:,3)-ye(:,1) ze(:,3)-ze(:,1)]'; % Edge from vertex 1 to 3

% Compute cross product to get area vectors.
c = cross(v1,v2); % Cross product: v1 x v2 for each triangle

% Calculate triangle areas.
area = sqrt(sum(c.*c))'; % Magnitude of cross product (area * 2)
area = area/2; % Divide by 2 to get actual triangle areas

% Compute triangle centroids if requested.
if nargout>1
  cg = [mean(xe,2) mean(ye,2) mean(ze,2)]; % Centroids as mean of vertex coordinates
end

return % Exit function