function [ip,ierr] = psect(pp1,pn1,pp2,pn2,tol)
%PSECT    Finds the intersection of two planes.
%
%         [IP,IERR] = PSECT(PP1,PN1,PP2,PN2) finds the intersection of
%         two planes defined by a point (PP1, PP2) and a normal vector
%         (PN1, PN2).
%
%         PSECT returns two points (IP) from the intersection line.
%         IERR is set to one (1) if there was an error.
%
%         11-Sep-96 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error('PSECT requires four input arguments.');
end
%
if (nargin<5)
  tol = 1e-8;
end
%
% Check Vectors
%
pp1 = pp1(:);
pn1 = pn1(:);
pp2 = pp2(:);
pn2 = pn2(:);
%
[n1 l1] = size(pp1);
[n2 l2] = size(pn1);
[n3 l3] = size(pp2);
[n4 l4] = size(pn2);
%
if ((l1~=1)|(l2~=1)|(l3~=1)|(l4~=1))
  error('PSECT only works with vectors.')
end
%
% Check that the Inputs have Three Rows
%
if ((n1~=3)|(n2~=3)|(n3~=3)|(n4~=3))
  error('Point and vector inputs must be of length three (3).')
end
%
% Initialize IERR
%
ierr = 0;
%
% Solve for Point in X-Y Plane
%
A = [pn1(1:2)'; pn2(1:2)'];
b = [pn1'*pp1; pn2'*pp2];
ip = [(A\b)' 0];
%
% Solve for Point in Y-Z Plane
%
A = [pn1(2:3)'; pn2(2:3)'];
ip = [ip; 0 (A\b)'];
%
% Check for NaNs
%
inan = isnan(ip);
inan = find(any(inan,2));
if ~isempty(inan)
%
% Solve for Point in X-Z Plane
%
  A = [pn1([1 3])'; pn2([1 3])'];
  b = [pn1'*pp1; pn2'*pp2];
  c = (A\b)';
  ip(inan,:) = [c(1) 0 c(2)];
end
%
% Check Intersection
%
if (abs(pn1'*(ip(1,:)'-pp1))>tol)|(abs(pn1'*(ip(2,:)'-pp1))>tol)| ...
   (abs(pn2'*(ip(1,:)'-pp2))>tol)|(abs(pn2'*(ip(2,:)'-pp2))>tol)
%
  ierr = 1;
  disp(' *** PSECT:  Intersection points not in planes.')
end
%
return