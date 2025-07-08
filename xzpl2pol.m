function [varargout] = xzpl2pol(varargin)
%XZPL2POL Transforms X-Z plane coordinates to polar coordinates (Theta-
%         R) in MRI knee data cell arrays.
%
%         DATT = XZPL2POL(DAT) Transforms the X-Z coordinates in MRI
%         knee data cell array, DAT, to polar coordinates (Theta-R).
%         The polar coordinates are returned in a cell array, DATT.
%
%         [DATT1,DATT2,DATT3,...] = XZPL2POL(DAT1,DAT2,DAT3, ...)
%         returns additional cell arrays with the same polar
%         transformation.
%
%         NOTES:  1.  Inputs must be a MRI knee data cell arrays with
%                 the coordinate matrices for each slice in rows in the
%                 cell array.
%
%                 2.  Note the X, Y and Z coordinates are transformed
%                 to the Theta (radians), R and Y coordinates.
%
%                 3.  Note that positive angles greater than 90 degrees
%                 (pi/2 radians and Z-axis) are offset by -2*pi radians
%                 (-360 degrees).
%
%                 4.  See pol2xypl.m for the inverse transformation.
%
%                 5.  M-file sl_info.m must be in the current path or
%                 directory.
%
%         31-Jul-2014 * Mack Gardner-Morse
%

% Check if at least one input cell array is provided.
if nargin<1
  error([' *** ERROR in XZPL2POL:  At least one input MRI knee', ...
         ' data cell array is required!']);
end

% Verify that the number of inputs matches the number of outputs.
if nargout~=nargin
  error([' *** ERROR in XZPL2POL:  Number of inputs do not match', ...
        ' the number of outputs!']);
end

% Process each input cell array in reverse order to assign outputs.
iloop = nargin;
for k = iloop:-1:1
   % Extract current cell array from input arguments.
   dat = varargin{k};

   % Get slice information using sl_info.
   [nslice,~,isd] = sl_info(dat); % nslice: number of slices, isd: slice indices

   % Convert cell array to matrix and transform to polar coordinates.
   xyz = cell2mat(dat); % Combine all slice points into Nx3 matrix
   [theta,r,z] = cart2pol(xyz(:,1),xyz(:,3),xyz(:,2)); % Transform X-Z to Theta-R, keep Y as Z

   % Adjust angles greater than 90 degrees (pi/2 radians).
   id = find(theta>pi/2); % Find indices where theta > pi/2
   if ~isempty(id)
     theta(id) = theta(id)-2*pi; % Subtract 2*pi to adjust angles
   end

   % Combine transformed coordinates into a matrix.
   trz = [theta r z]; % Form Nx3 matrix with Theta, R, Y coordinates

   % Reorganize coordinates back into cell array by slices.
   datt = cell(nslice,1); % Initialize output cell array
   for l = 1:nslice
      datt{l} = trz(isd(l)+1:isd(l+1),:); % Assign points to l-th slice
   end

   % Store transformed cell array in output.
   varargout{k} = datt; % Assign to k-th output
end

return % Exit function