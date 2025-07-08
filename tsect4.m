function [ip,il,ierr] = tsect4(v1,v2,v3,lp,lv,tol)
%TSECT4 Finds the intersection of a triangle and a line.
%
%       [IP,IL,IERR] = TSECT4(V1,V2,V3,LP,LV) finds the intersection of
%       a plane triangle defined by the vertices V1, V2 and V3 and a
%       line defined by a point LP and a direction vector LV.
%
%       TSECT4 returns the intersection point IP.  IL is set to true
%       (1) if an intersection was found.  IERR is set to true (1) if
%       the line lies (or is close to lying) in the plane of the
%       triangle.
%
%       [IP,IL,IERR] = TSECT4(V1,V2,V3,LP,LV,TOL) checks the
%       determinate to make sure it is not within tolerance TOL of zero
%       indicating the line is in (or close to) the plane of the
%       triangle.  Default tolerance is 1e-8.
%
%       NOTES:  1.  Based on the algorithm described in:
%               Tomas Moller and Ben Trumbore:  Fast, minimum storage
%               ray/triangle intersection.  Journal of Graphics Tools
%               2(1):21-28, 1997.
%
%               2.  Must have the M-file xprod.m in the current path or
%               directory.
%
%               3.  See also:
% http://www.mathworks.com/matlabcentral/fileexchange/25058-raytriangle-intersection
%
%               4.  Assumes the line is infinitely long.
%
%       10-Jul-2013 * Mack Gardner-Morse

% Check if at least five input arguments are provided.
if (nargin<5)
  error(' *** Error in TSECT4:  Five input arguments are required!');
end

% Set default tolerance if not provided or empty.
if (nargin<6)||isempty(tol)
  tol = 1e-8;
end

% Convert input vectors to column vectors for consistent calculations.
v1 = v1(:); % Vertex 1 of triangle
v2 = v2(:); % Vertex 2 of triangle
v3 = v3(:); % Vertex 3 of triangle
lp = lp(:); % Point on the line
lv = lv(:); % Direction vector of the line

% Initialize output variables.
il = false;   % Intersection flag (true if intersection found)
ierr = false; % Error flag (true if line lies in triangle's plane)
ip = [];      % Intersection point (empty if no intersection)

% Calculate triangle edges from vertices.
e1 = v2-v1; % Edge 1: from vertex 1 to vertex 2
e2 = v3-v1; % Edge 2: from vertex 1 to vertex 3

% Calculate determinant using cross product of line direction and edge 2.
p = xprod(lv,e2); % Cross product: lv x e2
d = p*e1;         % Dot product: (lv x e2) * e1

% Check if determinant is near zero, indicating line is in triangle's plane.
if abs(d)<tol
  ierr = true; % Set error flag
  return       % Exit function
end

% Calculate translation vector from vertex 1 to line point.
tv = lp-v1; % Translation vector: lp - v1

% Calculate barycentric coordinate u and check bounds.
u = p*tv/d; % u = ((lv x e2) * tv) / d
if (u<0)||(u>1) % u must be between 0 and 1 for intersection
  return        % No valid intersection, exit
end

% Calculate barycentric coordinate v and check bounds.
q = xprod(tv,e1); % Cross product: tv x e1
v = q*lv/d;       % v = ((tv x e1) * lv) / d
if (v<0)||(u+v>1) % v >= 0 and u+v <= 1 for intersection
  return          % No valid intersection, exit
end

% Compute intersection point using barycentric coordinates.
ip = (1-u-v)*v1+u*v2+v*v3; % ip = (1-u-v)*v1 + u*v2 + v*v3
il = true;                      % Set intersection flag

return % Exit function
