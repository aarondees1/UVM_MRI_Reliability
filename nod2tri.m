function [ir,nt] = nod2tri(nodlst,tric,ncn)
%NOD2TRI  Finds the triangles connected to a vector of nodes.
%
%         IR = NOD2TRI(NODLST,TRIC) returns the row index to triangles
%         connected to the nodes in the vector, NODLST, that are in the
%         three (3) column triangle connectivity matrix, TRIC.
%
%         [IR,NT] = NOD2TRI(NODLST,TRIC) returns the number of
%         triangles, NT, connected to the nodes in NODLST.
%
%         IR = NOD2TRI(NODLST,TRIC,NCN) returns the row index to
%         triangles connected to the nodes in NODLST with greater
%         than NCN nodes from NODLST in the triangle.  NCN must be
%         between zero (0) and two (2) since there are only three
%         nodes in a triangle.  The default is zero (0) (at least
%         one node from the triangle is in the input list).
%
%         NOTES:  1.  Note that the number of connected nodes from
%                 NODLST must be greater than the number of connected
%                 nodes, NCN, in a triangle.
%
%         27-Sep-2010 * Mack Gardner-Morse
%

% Check if at least two input arguments are provided.
if (nargin<2)
  error(' *** ERROR in NOD2TRI:  Not enough input arguments!');
end

% Set default number of connected nodes if not provided.
if (nargin<3)
  ncn = 0; % Default: at least one node in triangle must match nodlst
end

% Validate triangle connectivity matrix: must have three columns.
if size(tric,2)~=3
  error([' *** ERROR in NOD2TRI:  Triangle connectivity matrix,', ...
         ' TRIC, must have three (3) columns!']);
end

% Validate number of connected nodes: must be between 0 and 2.
if ncn<0|ncn>2
  error([' *** ERROR in NOD2TRI:  Number of connected nodes, NCN,' ...
         ' must be between zero (0) and two (2)!']);
end

% Ensure nodlst is a column vector and get its size.
nodlst = nodlst(:); % Convert node list to column vector
nn = size(nodlst,1); % Number of nodes in nodlst

% Get number of triangles.
nt = size(tric,1); % Number of triangles in tric

% Initialize logical array to track node connections.
nnc = false(nt,3); % Logical array: true where tric matches nodlst nodes
for k = 1:nn
   nnc = nnc|tric==nodlst(k); % Mark true where triangle nodes match current node
end

% Find triangles with more than ncn matching nodes.
ir = find(sum(nnc,2)>ncn); % Row indices of triangles with > ncn matches
nt = size(ir,1); % Number of matching triangles

return % Exit function