function [nx,ny,nz,xc,yc,zc] = tri_norm(tri,xyz)
%TRI_NORM Computes the normals and centroids of 3-D triangles.
%
%         [Nx,Ny,Nz,Xc,Yc,Zc] = TRI_NORM(TRI,XYZ) returns the normals,
%         [Nx,Ny,Nz] and centroids [Xc,Yc,Zc] of triangles defined by a
%         three (3) column connectivity matrix, TRI, and X, Y and Z
%         coordinates in a three (3) column matrix, XYZ.
%
%         NOTE:  1.  Must have more than one triangle.
%
%         21-Sep-2010 * Mack Gardner-Morse
%

% Check if both input arguments are provided.
if (nargin<2)
  error(' *** ERROR in TRI_NORM:  Two input arguments are required!');
end

% Validate input dimensions: tri and xyz must have three columns.
[~,ncol1] = size(tri);
ncol2 = size(xyz,2);
if (ncol1~=3)&&(ncol2~=3)
  error([' *** ERROR in TRI_NORM:  Input matrices must have three', ...
        ' (3) columns!']);
end

% Extract X, Y, Z coordinates from xyz matrix as row vectors.
x = xyz(:,1)'; % X coordinates
y = xyz(:,2)'; % Y coordinates
z = xyz(:,3)'; % Z coordinates

% Index triangle vertices using connectivity matrix tri.
xt = x(tri); % X coordinates of triangle vertices
yt = y(tri); % Y coordinates of triangle vertices
zt = z(tri); % Z coordinates of triangle vertices

% Calculate triangle edges for cross product.
a = [xt(:,3)-xt(:,2) yt(:,3)-yt(:,2) zt(:,3)-zt(:,2)]; % Edge from vertex 2 to 3
b = [xt(:,1)-xt(:,2) yt(:,1)-yt(:,2) zt(:,1)-zt(:,2)]; % Edge from vertex 2 to 1

% Compute triangle normals using cross product of edges.
n = cross(a,b,2); % Cross product: a x b for each triangle
n = n./repmat(sqrt(sum(n.^2,2)),1,3); % Normalize to unit vectors

% Extract normal components.
nx = n(:,1); % X components of normals
ny = n(:,2); % Y components of normals
nz = n(:,3); % Z components of normals

% Calculate triangle centroids as mean of vertex coordinates.
xc = mean(xt,2); % X coordinates of centroids
yc = mean(yt,2); % Y coordinates of centroids
zc = mean(zt,2); % Z coordinates of centroids

return % Exit function