function sides = meshbnd5(tri)
%MESHBND5 Finds the boundary edges of a triangle connectivity matrix.
%
%         SIDES = MESHBND5(TRI) given a three (3) column triangle
%         connectivity matrix, TRI, returns a two (2) column list of
%         nodes forming boundary edges, SIDES.
%
%         NOTES:  1.  See MESHBND3.M, MESHBND4.M and MESHBND2.M for
%                 similar methods where only the boundary nodes are
%                 returned.
%
%         06-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<1
  error([' *** ERROR in MESHBND5:  A triangle connectivity matrix' ...
         ' is required as an input!']);
end
%
[~,nc] = size(tri);
if nc~=3
  error([' *** ERROR in MESHBND5:  Triangle connectivity matrix' ...
         ' must have three columns!']);
end
%
% Get Sides of the Triangles
%
sides = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
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
end
%
return