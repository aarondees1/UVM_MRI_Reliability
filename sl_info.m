function [nslice,nsl,isl] = sl_info(dat)
%SL_INFO  Return MRI slice information from MRI cell arrays.
%
%         NSLICE = SL_INFO(DAT) given the MRI data cell array, DAT,
%         returns the number of slices in the data.
%
%         [NSLICE,NSL,ISL] = SL_INFO(DAT) returns a vector with the
%         number of points in each slice, NSL, and an index vector, ISL,
%         such that ISL provides an index into the points in a slice
%         when the data points are one coordinate array (e.g. the index
%         for the points in the k slice is ISL(k)+1:ISL(k+1)).
%
%         NOTES:  1.  Input must be a MRI data cell array with the
%                 coordinate matrices for each slice in rows in the
%                 cell array.
%
%         01-Jul-2014 * Mack Gardner-Morse
%

% Check if input data is provided.
if (nargin<1)
  error(' *** ERROR in SL_INFO:  No input MRI cell array!');
end

% Validate that input is a cell array.
if ~iscell(dat)
  error(' *** ERROR in SL_INFO:  Input must be a MRI cell array!');
end

% Ensure dat is a column cell array.
dat = dat(:); % Convert dat to column vector for consistent row-wise slice access

% Get number of slices.
nslice = size(dat,1); % Total number of slices in the cell array

% Initialize vector for number of points per slice.
nsl = zeros(nslice,1); % Array to store number of points in each slice
for k = 1:nslice
   nsl(k) = size(dat{k},1); % Number of points in k-th slice
end

% Compute cumulative index vector for slice point indexing.
isl = [0; cumsum(nsl)]; % Cumulative sum of points, starting with 0

return % Exit function