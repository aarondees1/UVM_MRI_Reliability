function [sides,ns,out] = meshbndq(quad)
%MESHBNDQ Finds the boundary edges of a quadrilateral connectivity
%         matrix.
%
%         SIDES = MESHBNDQ(QUAD) given a four (4) column quadrilateral
%         connectivity matrix, QUAD, returns a two (2) column list of
%         nodes forming boundary edges, SIDES.
%
%         [SIDES,NS] = MESHBNDQ(QUAD) returns the number of sides, NS.
%
%         [SIDES,NS,OUT] = MESHBNDQ(QUAD) returns a ordered list of
%         boundary nodes, OUT.
%
%         NOTES:  1.  See MESHBND4.M and MESHBND5.M for similar methods
%                 for triangular connectivity matrices.
%
%         28-Mar-2014 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<1
  error([' *** ERROR in MESHBNDQ:  A quadrilateral connectivity' ...
         ' matrix is required as an input!']);
end
%
[~,nc] = size(quad);
if nc~=4
  error([' *** ERROR in MESHBNDQ:  quadrilateral connectivity' ...
         ' matrix must have four (4) columns!']);
end
%
% Get Sides of the Quadrilaterals
%
sides = [quad(:,1:2); quad(:,2:3); quad(:,3:4); quad(:,[4 1])];
sides = sort(sides,2);
ns = size(sides,1);
%
% Find Duplicate Edges (Duplicate Edge == Interior Edge)
%
sides = sortrows(sides);
ds = diff(sides);
ds = sum(abs(ds),2);
idup = find(ds==0);
%
% Delete Duplicate Edges
%
if ~isempty(idup)
  idup = [idup; idup+1];
  idx = true(ns,1);
  idx(idup) = false;
  sides = sides(idx,:);
  ns = size(sides,1);
end
%
% Get Unique and Sorted Node IDs
%
if nargout>2
  out = unique(sides(:));
end
%
return