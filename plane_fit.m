function [pxyz,norm_vec,coeff,v,score,pexp,res,sse] = plane_fit(x,y,z)
%PLANE_FIT Fits a plane to X, Y and Z point coordinate data using SVD.
%
%          [PXYZ,NORM_VEC,COEFF] = PLANE_FIT(X,Y,Z) given the X, Y
%          and Z coordinates of a set of points, calculates a least
%          squares fit of a plane. A point on the plane, PXYZ, the
%          normal vector, NORM_VEC and the coefficients of the
%          algebraic equation for a plane, COEFF, are returned.
%          The algebraic form is:
%          coeff(1)*x + coeff(2)*y + coeff(3)*z + coeff(4) = 0.
%
%          [PXYZ,NORM_VEC,COEFF,V,SCORE,PEXP,RES,SSE] = PLANE_FIT(X,Y,Z)
%          returns the rotation matrix, V, PCA scores, SCORE, percent
%          of variance explained by the three orthogonal directions of
%          the plane, PEXP, the residuals (difference between the data
%          and fitted plane), RES, and the sum of squared errors, SSE.
%
%          NOTES:  1.  Must have at least three (3) points.
%
%                  2.  The SVD is used to do a principal component
%                  analysis (PCA) to do an orthogonal regression (total
%                  least squares) fit of the plane.
%
%                  3.  Based on the demonstration of orthogonal
%                  regression using Matlab Statistics Toolbox.  See:
%                  http://www.mathworks.com/products/statistics/
%                  demos.html?file=/products/demos/shipping/stats/
%                  orthoregdemo.html
%
%          22-June-2010 * Mack Gardner-Morse
%

% Check if three input arguments (x, y, z coordinates) are provided.
if (nargin<3)
  error([' *** ERROR in PLANE_FIT:  The X, Y and Z coordinates of ', ...
         'the points to be fit are required as inputs!']);
end

% Combine x, y, z into a matrix and verify sufficient points.
xyz = [x(:) y(:) z(:)]; % Form Nx3 matrix of point coordinates
npts = size(xyz,1); % Number of points
if npts<3
  error(' *** ERROR in PLANE_FIT:  Not enough data points!');
end

% Center the data by subtracting the mean point.
pxyz = mean(xyz); % Mean point (on the fitted plane)
xyz = xyz-repmat(pxyz,npts,1); % Subtract mean to center data

% Perform SVD to fit the plane using PCA.
[u,s,v] = svd(xyz); % Singular value decomposition: u*s*v' = xyz

% Extract normal vector from the third column of V (smallest singular value).
norm_vec = v(:,3); % Normal vector to the fitted plane

% Compute coefficients for the plane equation: coeff(1)*x + coeff(2)*y + coeff(3)*z + coeff(4) = 0.
coeff = [norm_vec; -pxyz*norm_vec]; % [nx, ny, nz, -dot(pxyz, norm_vec)]

% Compute additional outputs if requested.
if nargout>3
  score = u*s; % PCA scores (projections of centered data onto principal components)
  rts = diag(s); % Singular values
  pexp = 100*rts./sum(rts); % Percent variance explained by each principal component
  res = xyz-score(:,1:2)*v(:,1:2)'; % Residuals: difference from plane in 3D space
  err = xyz*norm_vec; % Distances from points to plane along normal
  sse = sum(err.^2); % Sum of squared errors (distances)
end

return % Exit function