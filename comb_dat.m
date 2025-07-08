function dats = comb_dat(datl,datm,datt)
%COMB_DAT Combines the sagittal lateral and medial condyle and trochlea
%         MRI knee data cell arrays into a single data cell array.
%
%         DATS = COMB_DAT(DATL,DATM,DATT) given the data cell arrays
%         containing three (3) columns matrices with slice coordinate
%         point data for the lateral condyle, DATL, medial condyle,
%         DATM, and trochlea, DATT, returns the cell array, DATS, that
%         contains the combined ordered slice data from the three input
%         data cell arrays.
%
%         NOTES:  1.  To ensure the correct ordering of the slices,
%                 the input cell arrays must be input in the correct
%                 order (lateral, medial and trochlea).
%
%                 2.  The function assumes all the slices are in order
%                 within the data cell arrays.
%
%                 3.  The M-files plane_fit.m and sl_dir.m must be in
%                 the current path or directory.
%
%         03-Jul-2014 * Mack Gardner-Morse
%

% Check if three input arguments are provided.
if (nargin<3)
  error(' *** ERROR in COMB_DAT:  No or not enough input data!');
end

% Convert cell arrays to matrices for calculating mean centers.
xyzl = cell2mat(datl); % Combine lateral condyle slices into matrix
xyzm = cell2mat(datm); % Combine medial condyle slices into matrix

% Compute direction vector between mean centers of lateral and medial condyles.
vx = mean(xyzl)-mean(xyzm); % Vector from medial to lateral condyle center

% Determine slice directions for each data set using sl_dir.
slvm = sl_dir(datm); % Slice direction for medial condyle
slvt = sl_dir(datt); % Slice direction for trochlea
slvl = sl_dir(datl); % Slice direction for lateral condyle

% Reverse slice order if direction opposes medial-to-lateral vector.
if slvm*vx'<0
  datm = flipud(datm); % Flip medial condyle slices
end
if slvt*vx'<0
  datt = flipud(datt); % Flip trochlea slices
end
if slvl*vx'<0
  datl = flipud(datl); % Flip lateral condyle slices
end

% Combine slice data cell arrays in order: medial, trochlea, lateral.
dats = [datm; datt; datl]; % Concatenate cell arrays vertically

return % Exit function