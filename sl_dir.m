function [slvec,slx] = sl_dir(dat,iall)
%SL_DIR   Gets the slice direction between MRI slices in the ordered 
%         slice data from the digitized MRI knee cell arrays.
%
%         SLVEC = SL_DIR(DAT) given a cell array containing three (3)
%         columns matrices with slice coordinate point data, DAT,
%         returns the three (3) column slice separation vector, SLVEC,
%         between the second and first slice.
%
%         [SLVEC,SLX] = SL_DIR(DAT) returns the magnitude (separation
%         distance) of the separation vector, SLX.
%
%         [SLVEC,SLX] = SL_DIR(DAT,IALL) if the logical scalar, IALL, is
%         true, returns the separation vector, SLVEC, between each slice
%         and the separation distance, SLX between each slice in matrix
%         with the different slice data in rows.
%
%         NOTES:  1.  The M-file plane_fit.m must be in the current path
%                 or directory.
%
%                 2.  Only checks the direction between the first two
%                 slices.
%
%         03-Jul-2014 * Mack Gardner-Morse
%

% Check if input data is provided.
if (nargin<1)
  error(' *** ERROR in SL_DIR:  No input data!');
end

% Set default for iall to false if not provided.
if (nargin<2)
  iall = false; % Default: process only first two slices
end

% Ensure dat is a column cell array and determine number of slices to process.
dat = dat(:); % Convert dat to column vector
if iall
  nslice = size(dat,1); % Process all slices if iall is true
else
  nslice = 2; % Process only first two slices
end

% Initialize arrays for slice centers and normal vectors.
cntr = zeros(nslice,3); % Centers of fitted planes
nvec = zeros(3,nslice); % Normal vectors of fitted planes

% Fit planes to each slice and compute additional vectors.
for k = 1:nslice
   xyz = dat{k}; % Extract slice data (Nx3 matrix)
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3)); % Fit plane to slice
   vec = xyz(end,:)-xyz(1,:); % Vector from first to last point in slice
   vec1(k,:) = vec./norm(vec); % Normalize vector (not used in output)
end

% Initialize arrays for slice separation vectors and distances.
slx = zeros(nslice-1,1); % Separation distances between consecutive slices
slvec = zeros(nslice-1,3); % Separation vectors between consecutive slices

% Compute separation vectors and distances between consecutive slices.
for k = 2:nslice
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1); % Distance along normal of previous slice
   slvec(k-1,:) = slx(k-1)*nvec(:,k-1)'; % Separation vector using normal of previous slice
end

% Ensure separation distances are positive.
slx = abs(slx); % Convert distances to absolute values

return % Exit function