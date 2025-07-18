function [p1,p2,dist] = li_clos(li1,li2) % Define function to find closest points between two 3D lines
%LI_CLOS Finds the closest two points between two nonparallel 3D lines.
%
%        [P1,P2] = LI_CLOS(LI1,LI2) given a two (2) rows by three (3)
%        columns matrix giving 3D coordinates for two points on one
%        line, LI1, and another two (2) rows by three (3) columns
%        matrix giving 3D coordinates for two points on a second line,
%        LI2, returns the three (3) column row vector of the coordinates
%        for the point on the first line that is closest to the second
%        line, P1, and the the three (3) column row vector of the
%        coordinates for the point on the second line that is closest to
%        the first line, P2.  Note that the line between P1 and P2 is
%        also perpendicular to both the input lines.
%
%        [P1,P2,DIST] = LI_CLOS(LI1,LI2) returns the closest distance
%        between the two lines, DIST.
%
%        NOTES:  1.  Based on the algorithm at:
%        http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
%
%                2.  Does not check that the two points on the lines are
%                not the same.
%
%        17-Aug-2010 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<2) % Check if fewer than 2 inputs are provided
  error([' *** ERROR in LI_CLOS:  At least two (2) input lines', ...
         ' are required!']); % Throw error for insufficient inputs
end
%
% Check Inputs
%
[r1,c1] = size(li1); % Get dimensions of first line input
[r2,c2] = size(li2); % Get dimensions of second line input
%
if r1~=2|r2~=2 % Check if each line has 2 points
  error([' *** ERROR in LI_CLOS:  Two (2) points are', ...
         ' required to define each line!']); % Throw error for invalid point count
end
%
if c1~=3|c2~=3 % Check if points have 3 coordinates
  error([' *** ERROR in LI_CLOS:  Three (3) coordinates are', ...
         ' required to define each point!']); % Throw error for invalid coordinate count
end
%
% Get Line Direction Vectors
%
r0 = li1(1,:)-li2(1,:); % Compute vector between line origins
u = diff(li1); % Compute direction vector for first line
v = diff(li2); % Compute direction vector for second line
%
% Get Dot Products
%
a = u*u'; % Compute dot product of first direction vector
b = u*v'; % Compute dot product between direction vectors
c = v*v'; % Compute dot product of second direction vector
d = u*r0'; % Compute dot product of first direction and origin vector
e = v*r0'; % Compute dot product of second direction and origin vector
%
% Solve for Line Parameters
%
den = a*c-b*b; % Compute denominator for line parameters
tol = eps^(2/3); % Set tolerance for parallel check
if abs(den)<tol % Check if lines are nearly parallel
  error([' *** ERROR in LI_CLOS:  The two lines are parallel', ...
         ' or almost parallel!']); % Throw error for parallel lines
else
  sc = (b*e-c*d)./den; % Compute parameter for first line
  tc = (a*e-b*d)./den; % Compute parameter for second line
end
%
% Solve for the Points
%
p1 = li1(1,:)+sc*u; % Compute closest point on first line
p2 = li2(1,:)+tc*v; % Compute closest point on second line
%
% Get Shortest Distance between Lines
%
if nargout>2 % Check if distance output is requested
  dist = p2-p1; % Compute vector between closest points
  dist = sqrt(dist*dist'); % Compute Euclidean distance
end
%
return % Exit the function