function [xyz,npts] = fix_pts_AD(xyz,tol,iflag,fstr); % Define function to remove duplicate 3D points with custom message
%FIX_PTS Removes duplicate points from a matrix of 3-D data.
%
%        [XYZ,NPTS] = FIX_PTS(XYZ,TOL) given a three (3) columns matrix
%        with coordinate point data, XYZ, returns the matrix without
%        any duplicate points within a distance, TOL.  The number of
%        returned points, NPTS, may also be returned.
%
%        [XYZ,NPTS] = FIX_PTS(XYZ,TOL,IFLAG) if IFLAG is true (or not
%        equal to zero (0)), prints a message to the screen if any
%        duplicates are found.
%
%        NOTES:  None.
%
%        25-June-2010 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<3) % Check if flag input is missing
  iflag = false; % Set default flag to false
end
%
if (nargin<2) % Check if fewer than 2 inputs are provided
  error(' *** ERROR in FIX_PTS:  Must have two inputs!'); % Throw error for insufficient inputs
end
%
[npts,ncol] = size(xyz); % Get number of points and columns
if ncol~=3 % Check if columns are not 3
  error(' *** ERROR in FIX_PTS:  Data must have three columns!'); % Throw error for invalid column count
end
%
% Find and Remove Duplicates
%
tol2 = tol.*tol; % Calculate squared tolerance
npts = size(xyz,1); % Update number of points
for k = 1:npts-1 % Loop over points
   l = k+1:npts; % Get indices for remaining points
   dxyz = xyz(l,:)-repmat(xyz(k,:),npts-k,1); % Compute differences from current point
   dist2 = sum(dxyz.*dxyz,2); % Calculate squared distances
   idx = find(dist2<tol2); % Find duplicates within tolerance
   if ~isempty(idx) % Check if duplicates are found
     if iflag % Check if flag is true
       disp(fstr); % Display custom message
       fprintf(1,'\n *** Duplicate Points Found!\n'); % Print duplicate message
%        error('\n *** Duplicate Points Found!\n');
       iflag = false; % Disable flag after first message
     end
     idx = [k; l(idx)']; % Combine current and duplicate indices
     ndup = size(idx,1); % Get number of duplicates
     xyz(k,:) = mean(xyz(idx,:)); % Average duplicate points
     xyz(idx(2:ndup),:) = repmat(NaN,ndup-1,3); % Mark duplicates as NaN
   end
end
%
xyz = xyz(~isnan(xyz(:,1)),:); % Remove NaN rows
npts = size(xyz,1); % Update number of points
%
return % Exit the function