function [xyzp,t] = pts2lin(lpt0,lvec,pxyz)
%PTS2LIN  Calculates the coordinates of the points on a line that are
%         closest to the points not on the line.  These points also
%         form perpendicular distances to the line from the points not
%         on the line.
%
%         XYZP = PTS2LIN(LPT0,LVEC,PXYZ) given a line defined by an
%         initial point on the line in a row vector, LPT0, with X, Y
%         and Z coordinates in columns, a direction vector for the line
%         in a row vector, LVEC, with  X, Y and Z components and the
%         coordinates of the points not on the line in a matrix,
%         PXYZ, with X, Y and and Z coordinates in columns, calculates
%         the coordinates of points on the line that are closest to
%         the points not on the line.  These point also form
%         perpendicular distances to the line from the points not
%         on the line.  The X, Y and Z coordinates of the points on the
%         line are returned in the matrix, XYZP.
%
%         [XYZP,T] = PTS2LIN(LPT0,LVEC,PXYZ) returns the parametric
%         distances in column vector, T, along the line to the points on the line that
%         are closest to the points not on the line.
%
%         NOTES:  1.  The line is assumed to have infinite length.
%
%                 2.  All points and the vector must have three columns.
%
%         17-July-2013 * Mack Gardner-Morse
%

% Check if three input arguments are provided.
if (nargin<3)
  error([' *** ERROR in PT2LINE:  At least three (3) inputs', ...
         ' are required!']);
end

% Validate input dimensions: lpt0 and lvec must be 1xN, pxyz must have N columns, where N is 2 or 3.
[nrow1,ncol1] = size(lpt0);
[nrow2,ncol2] = size(lvec);
[nrow3,ncol3] = size(pxyz);
if nrow1~=1|nrow2~=1|ncol1~=ncol2|ncol2~=ncol3|(ncol1~=3&ncol1~=2)
  error([' *** ERROR in PTS2LIN:  All input points and the vector', ...
         ' must have two or three columns!']);
end

% Compute the closest points on the line to the given points.
lpt0mat = repmat(lpt0,nrow3,1); % Replicate initial line point to match number of input points
dp = pxyz-lpt0mat; % Vectors from line's initial point to each input point
t = ((lvec*dp')/(lvec*lvec'))'; % Parametric distances along line to closest points
xyzp = t*lvec+lpt0mat; % Coordinates of the closest points on the line

return % Exit function