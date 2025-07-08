function [nodv,arean,triv,areat,cg] = nod_norm(tri,xyz,iplt)
%NOD_NORM  Calculates corner point (node) normals of a triangular
%          surface mesh. 
%
%          NODV = NOD_NORM(TRI,XYZ) returns a matrix of normal vectors
%          at the corner points (nodes) of a triangular surface mesh in
%          the three (3) column matrix, NODV.  The triangular surface
%          mesh is defined by a three (3) column triangle connectivity
%          matrix, TRI, and the X, Y and Z coordinates of the nodes in
%          a three (3) column matrix, XYZ.
%
%          [NODV,AREAN,TRIV,AREAT,CG] = NOD_NORM(TRI,XYZ) returns the
%          areas of the triangles connected to the node in the vector,
%          AREAN, normal vectors of the triangles in a three (3) column
%          matrix, TRIV, the areas of the triangles in the vector,
%          AREAT, and the X, Y and Z coordinates of the area centers of
%          the triangles in the three (3) column matrix CG.
%
%          NOTES:  None.
%
%          18-April-2011 * Mack Gardner-Morse
%

% Check if at least two input arguments are provided.
if (nargin<2)
  error(' *** ERROR in nod_norm:  Not enough input arguments.');
end

% Set default plotting flag to false if not provided.
if (nargin<3)
  iplt = false;
end

% Validate input dimensions: tri and xyz must have three columns.
ncol1 = size(tri,2);
ncol2 = size(xyz,2);
if (ncol1~=3)&(ncol2~=3)
  error([' *** ERROR in nod_norm:  Input matrices must have three', ...
        ' (3) columns!']);
end

% Get number of triangles and nodes.
nt = size(tri,1); % Number of triangles
nnod = size(xyz,1); % Number of nodes

% Extract X, Y, Z coordinates from xyz matrix.
x = xyz(:,1); % X coordinates
y = xyz(:,2); % Y coordinates
z = xyz(:,3); % Z coordinates

% Index triangle vertices using connectivity matrix tri.
xe = x(tri); % X coordinates of triangle vertices
ye = y(tri); % Y coordinates of triangle vertices
ze = z(tri); % Z coordinates of triangle vertices

% Transpose coordinates for single triangle case to match multi-triangle format.
if nt==1
  xe = xe';
  ye = ye';
  ze = ze';
end

% Calculate triangle edges for cross product.
v1 = [xe(:,2)-xe(:,1) ye(:,2)-ye(:,1) ze(:,2)-ze(:,1)]'; % Edge from vertex 1 to 2
v2 = [xe(:,3)-xe(:,1) ye(:,3)-ye(:,1) ze(:,3)-ze(:,1)]'; % Edge from vertex 1 to 3

% Compute triangle normals using cross product of edges.
triv = cross(v1,v2); % Cross product: v1 x v2 for each triangle

% Calculate triangle areas and normalize normals.
areat = sqrt(sum(triv.*triv))'; % Magnitude of normal vectors (area * 2)
rnorm = repmat(areat,1,3); % Replicate areas for normalization
triv = triv'./rnorm; % Normalize normal vectors to unit length
areat = areat/2; % Compute actual triangle areas
rnorm = repmat(areat,1,3); % Replicate areas for weighted normals

% Initialize node normals and areas.
nodv = zeros(nnod,3); % Node normal vectors
arean = zeros(nnod,1); % Sum of triangle areas per node

% Compute area-weighted triangle normals.
atriv = rnorm.*triv; % Weight normals by triangle areas

% Sum contributions of triangle normals and areas to each node.
for k = 1:nt
   idn = tri(k,:)'; % Node indices for current triangle
   nodv(idn,:) = nodv(idn,:)+repmat(atriv(k,:),3,1); % Add weighted normal
   arean(idn) = arean(idn)+repmat(areat(k),3,1); % Add triangle area
end

% Normalize node normals by area for nodes with non-zero areas.
idg = find(arean); % Indices of nodes with areas to avoid division by zero
nodv(idg,:) = nodv(idg,:)./repmat(arean(idg),1,3); % Divide by total area
rnorm = sqrt(sum(nodv(idg,:).^2,2)); % Compute magnitude of node normals
nodv(idg,:) = nodv(idg,:)./repmat(rnorm,1,3); % Normalize to unit length

% Calculate triangle centroids if requested.
if nargout>4
  cg = [mean(xe,2) mean(ye,2) mean(ze,2)]; % Centroids as mean of vertices
end

return % Exit function