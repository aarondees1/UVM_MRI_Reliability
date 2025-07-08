function [ip,t,il,ierr] = plsect(pp,pn,lp,lv,tol) % Define function to find plane-line intersection
%PLSECT  Finds the intersection of a plane and a line.
%
%        [IP,T,IL,IERR] = PLSECT(PP,PN,LP,LV) finds the intersection of
%        a plane defined by a point on the plane, PP, and a vector
%        normal to the plane, PN, and a line defined by a
%        point LP and a direction vector LV.  The lengths of all the
%        inputs must be three (3) for the X, Y, and Z components of the
%        coordinates or vectors.  IP is the coordinates of the
%        intersection, T is the parametric distance along the line to
%        the intersection.  IL is set to true (1) if an intersection
%        was found.  IERR is set to true (1) if the line is (or close
%        to) parallel to the plane.  If there is no intersection, IP
%        and T are empty arrays.
%        
%        IP = PLSECT(PP,PN,LP,LV,TOL) checks to make make sure the dot
%        product of the line direction and plane normal is not within
%        tolerance TOL of zero indicating the line is (or close to)
%        parallel to the plane.  Default tolerance is 1e-8.
%        
%        NOTES:  None.
%
%       02-Dec-2022 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<4) % Check if fewer than 4 inputs are provided
  error(' *** Error in PLSECT:  Four input arguments are required!'); % Throw error for insufficient inputs
end
%
if (nargin<5)||isempty(tol) % Check if tolerance is missing or empty
  tol = 1e-8; % Set default tolerance
end
%
% Get Column Vectors
%
pp = pp(:); % Convert plane point to column vector
pn = pn(:); % Convert plane normal to column vector
lp = lp(:); % Convert line point to column vector
lv = lv(:); % Convert line direction to column vector
%
% Initial Values
%
ip = []; % Initialize intersection point as empty
t = []; % Initialize parametric distance as empty
il = false; % Initialize intersection flag as false
ierr = false; % Initialize error flag as false
%
% Check Dot Product of the Line Direction and Plane Normal
%
if abs(lv'*pn)<tol % Check if line is nearly parallel to plane
  ierr = true; % Set error flag
  return; % Exit function
end
%
% Calculate Intersection Point
% Solution of substituting the parametric line equation into the plane
% equation.
%
t = pn'*(pp-lp)./(pn'*lv); % Compute parametric distance
%
% Check Intersection is Within the Line Endpoints (0<=t<=1)
%
if (0<=t)&&(t<=1) % Check if intersection is within segment
  il = true; % Set intersection flag
  ip = lp+t*lv; % Compute intersection point
  ip = ip'; % Transpose to row vector
end
%
return % Exit the function