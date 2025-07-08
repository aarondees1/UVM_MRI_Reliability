function [ip,t,idx] = plsect2(pp,pn,pwl,tol) % Define function to find plane-piecewise linear line intersection
%PLSECT2 Finds the intersection of a plane and a piecewise linear (PWL)
%        line.
%
%        IP = PLSECT2(PP,PN,PWL) finds the first intersection of a
%        plane defined by a point on the plane, PP, and a vector normal
%        to the plane, PN, and a piecewise linear (PWL) line.  PWL is
%        defined by a series of 3-D points with the X, Y and Z
%        coordinates of the points in columns.  The X, Y, and Z
%        coordinates of the first intersection with the piecewise
%        linear line are returned in IP.  IP is an empty array if there
%        is no intersection.
%
%        [IP,T,IDX] = PLSECT2(PP,PN,PWL) returns the distance, T, along
%        the line with the first intersection.  T is the normalized
%        (0 to 1) distance along the line to the intersection.  The
%        index, IDX, is to the first point of the segment in the
%        piecewise linear line (PWL) that intersects the plane.  If
%        there is no intersection, T and IDX are empty arrays.
%        
%        IP = PLSECT2(PP,PN,PWL,TOL) checks to make make sure the dot
%        product of the line directions and plane normal is not within
%        tolerance TOL of zero indicating the line is (or close to)
%        parallel to the plane.  Default tolerance is 1e-8.
%        
%        NOTES:  1.  The M-file plsect.m must be in the current path or
%                directory.
%
%       02-Dec-2022 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<3) % Check if fewer than 3 inputs are provided
  error(' *** Error in PLSECT2:  Three input arguments are required!'); % Throw error for insufficient inputs
end
%
if (nargin<4)||isempty(tol) % Check if tolerance is missing or empty
  tol = 1e-8; % Set default tolerance
end
%
% Get Column Vectors
%
pp = pp(:); % Convert plane point to column vector
pn = pn(:); % Convert plane normal to column vector
%
% Check Inputs
%
[npts,nc] = size(pwl); % Get dimensions of piecewise linear line
if size(pp,1)~=3||size(pn,1)~=3||nc~=3 % Check if coordinates are 3D
  error([' *** Error in PLSECT2:  Coordinate dimensions must equal', ...
         ' three (3)!']); % Throw error for invalid dimensions
end
if npts<2 % Check if line has at least 2 points
  error([' *** Error in PLSECT2:  Piecewise linear line must have', ...
         ' at least two (2) points!']); % Throw error for insufficient points
end
%
% Initial Values
%
ip = []; % Initialize intersection point as empty
t = []; % Initialize normalized distance as empty
idx = []; % Initialize segment index as empty
%
nl = npts-1; % Compute number of line segments
%
% Loop through Piecewise Linear Line
%
for k = 1:nl % Loop through segments
%
   l = k+1; % Get index of next point
%
   lp = pwl(k,:)' ; % Get current point as column vector
   lv = pwl(l,:)' -lp; % Compute line direction vector
%
% Check for Intersection
%
   [ipp,tp] = plsect(pp,pn,lp,lv,tol); % Find intersection with segment
%
   if ~isempty(ipp) % Check if intersection exists
     ip = ipp; % Store intersection point
     t = tp; % Store normalized distance
     idx = k; % Store segment index
     break; % Exit loop
   end
%
end
%
return % Exit the function