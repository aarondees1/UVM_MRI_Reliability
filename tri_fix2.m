function tri = tri_fix2(tri,xyz) % Define function to improve triangle aspect ratios
%TRI_FIX2 Trys to improve the aspect ratios of triangles in a
%         three-dimensional (3-D) triangular mesh.
%
%         TRI = TRI_FIX2(TRI,XYZ) given a three column triangular
%         connectivity matrix TRI and the X, Y and Z coordinates in the
%         columns of matrix XYZ, checks the sum of the opposite angles
%         of the triangles sharing a common edge and flips the edge if
%         the sum of the opposite angles is greater than pi.
%
%         NOTES:  1.  Only operates locally (pairs of triangles) to
%                 improve the aspect ratio.  Other more global
%                 functions or functions that change the number of
%                 triangles may work better to improve the mesh.
%
%                 2.  The coordinates of the vertices of the triangles
%                 must be a matrix with the vertices in the rows and
%                 the X, Y and Z coordinates in the columns.
%
%                 3.  The M-file tri_norm.m must be in the current path
%                 or directory.
%
%         27-Aug-2013 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<2) % Check if fewer than 2 inputs are provided
  error(' *** ERROR in TRI_FIX2:  TRI_FIX2 requires two inputs!'); % Throw error for insufficient inputs
end
%
% Check Inputs
%
[nt,nc1] = size(tri); % Get number of triangles and columns
[~,nc2] = size(xyz); % Get number of columns in coordinates
if nc1~=3||nc2~=3 % Check if matrices have 3 columns
  error([' *** ERROR in TRI_FIX2:  Triangle connectivity matrix' ...
         ' and coordinate matrix must have three columns!']); % Throw error for invalid column count
end
%
% Get Edges
%
edges = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])]; % Extract all triangle edges
edges = sort(edges,2); % Sort vertices within each edge
%
% ne = size(edges,1);                    % Number of edges
idt = repmat((1:nt)',3,1); % Create triangle indices for edges
%
% Find Duplicate or Common Edges
%
[edges,is] = sortrows(edges); % Sort edges by rows
idt = idt(is); % Reorder triangle indices
ids = repmat(1:3,nt,1); % Create side indices
ids = ids(:); % Flatten side indices
ids = ids(is); % Reorder side indices
%
ds = diff(edges); % Compute differences between consecutive edges
ds = sum(abs(ds),2); % Sum absolute differences
idup = find(ds==0); % Find indices of duplicate edges
nd = size(idup,1); % Get number of common edges
%
idc = [idup idup+1]'; % Create indices for paired edges
idc = idc(:); % Flatten paired indices
edges = edges(idc,:); % Select common edges
idt = idt(idc); % Select corresponding triangle indices
ids = ids(idc); % Select corresponding side indices
%
% Get Scores for Opposite Angles
%
sc = zeros(nd,1); % Initialize score array
%
for k = 1:nd % Loop over common edges
%
   l = 2*k; % Get index for second edge in pair
   n23 = edges(l,:); % Get common edge vertices
%
   n1 = setdiff(tri(idt(l-1),:),n23); % Get opposite vertex of first triangle
   p1 = xyz(n1,:); % Get coordinates of opposite vertex
   v1 = xyz(n23(1),:)-p1; % Compute vector to first vertex
   v2 = xyz(n23(2),:)-p1; % Compute vector to second vertex
   ang1 = v1*v2'; % Compute dot product
   ang1 = ang1./sqrt((v1*v1')*(v2*v2')); % Normalize dot product
   ang1 = acos(ang1); % Compute angle
%
   n2 = setdiff(tri(idt(l),:),n23); % Get opposite vertex of second triangle
   p2 = xyz(n2,:); % Get coordinates of opposite vertex
   v3 = xyz(n23(1),:)-p2; % Compute vector to first vertex
   v4 = xyz(n23(2),:)-p2; % Compute vector to second vertex
   ang2 = v3*v4'; % Compute dot product
   ang2 = ang2./sqrt((v3*v3')*(v4*v4')); % Normalize dot product
   ang2 = acos(ang2); % Compute angle
   sc(k) = ang1+ang2-pi; % Compute score for flipping
end
%
% Sort Scores
%
[ssc,iss] = sort(sc,1,'descend'); % Sort scores in descending order
%
% Flip Edges and Update Score
%
iter = 0; % Initialize iteration counter
while ssc(1)>0 % Loop while top score is positive
     iter = iter+1; % Increment iteration counter
%
% Form New Triangles by Flipping Edge
%
     l = 2*iss(1); % Get index for top-scoring edge pair
     n23 = edges(l,:); % Get common edge vertices
     idt1 = idt(l-1); % Get first triangle index
     n1 = setdiff(tri(idt1,:),n23); % Get opposite vertex of first triangle
     idt2 = idt(l); % Get second triangle index
     n2 = setdiff(tri(idt2,:),n23); % Get opposite vertex of second triangle
     nid1 = [n1 n23(1) n2]; % Form new first triangle
     nid2 = [n2 n23(2) n1]; % Form new second triangle
%
% Check Node Ordering of the Triangles to Preserve Triangles Normal
% Directions
%
     [xn,yn,zn] = tri_norm(tri([idt1;idt2],:),xyz); % Compute normals for original triangles
     nv = [mean(xn) mean(yn) mean(zn)]; % Compute mean normal vector
     nv = nv/norm(nv); % Normalize mean normal vector
%
     [xn,yn,zn] = tri_norm(nid1,xyz); % Compute normal for new first triangle
     nvn = [xn,yn,zn]; % Store normal vector
     ang = acos(nvn*nv'); % Compute angle between normals
     if ang>pi/2 % Check if normal is flipped
       nid1 = nid1([1 3 2]); % Reverse vertex order
     end
%
     [xn,yn,zn] = tri_norm(nid2,xyz); % Compute normal for new second triangle
     nvn = [xn,yn,zn]; % Store normal vector
     ang = acos(nvn*nv'); % Compute angle between normals
     if ang>pi/2 % Check if normal is flipped
       nid2 = nid2([1 3 2]); % Reverse vertex order
     end
%
% Update Triangle Connectivity Matrix, Edges, Triangle Index and Side Index
%
     tri(idt1,:) = nid1; % Update first triangle
     tri(idt2,:) = nid2; % Update second triangle
%
     ne1 = [tri(idt1,1:2); tri(idt1,2:3); tri(idt1,[3 1])]; % Get edges of first new triangle
     ne2 = [tri(idt2,1:2); tri(idt2,2:3); tri(idt2,[3 1])]; % Get edges of second new triangle
%
     ne1 = sort(ne1,2); % Sort vertices within edges
     ne2 = sort(ne2,2); % Sort vertices within edges
%
     [s,i1,i2] = intersect(ne1,ne2,'rows'); % Find common edges
     edges(l-1:l,:) = repmat(s,2,1); % Update common edges
     idt(l-1:l) = [idt1;idt2]; % Update triangle indices
     ids(l-1:l) = [i1;i2]; % Update side indices
%
     [s1,i1] = setdiff(ne1,ne2,'rows'); % Find unique edges in first triangle
     [s2,i2] = setdiff(ne2,ne1,'rows'); % Find unique edges in second triangle
     si = [i1; i2]; % Combine side indices
     [~,i3,i4] = intersect(edges,[s1;s2],'rows'); % Find matching edges
     i3 = i3+rem(i3,2); % Adjust edge indices
     idx = [i3-1 i3]'; % Create paired edge indices
     idx = idx(:); % Flatten indices
     idx2 = idt(idx)==idt1|idt(idx)==idt2; % Identify affected triangles
     idx2 = idx(idx2); % Select relevant indices
     idts = [idt1; idt1; idt2; idt2]; % Create triangle index array
     idts = idts(i4); % Select corresponding triangle indices
%     iter
% if iter>=120; keyboard; end
     idt(idx2) = idts; % Update triangle indices
     si = si(i4); % Reorder side indices
     ids(idx2) = si; % Update side indices
%
% Update Scores of Adjoining Edges
%
     sc(iss(1)) = -sc(iss(1)); % Negate score of flipped edge
     ids1 = ids(l-1); % Get side index of first triangle
     ids2 = ids(l); % Get side index of second triangle
     ide1 = find(idt==idt1&ids~=ids1); % Find other edges of first triangle
     ide2 = find(idt==idt2&ids~=ids2); % Find other edges of second triangle
     ide = unique([ide1; ide2]); % Combine unique edge indices
     ide = ide+rem(ide,2); % Adjust edge indices
     ns = size(ide,1); % Get number of affected edges
     if ns>0 % Check if there are affected edges
       for k = 1:ns % Loop over affected edges
          l = ide(k); % Get edge index
          n23 = edges(l,:); % Get common edge vertices
%
          n1 = setdiff(tri(idt(l-1),:),n23); % Get opposite vertex of first triangle
          p1 = xyz(n1,:); % Get coordinates of opposite vertex
          v1 = xyz(n23(1),:)-p1; % Compute vector to first vertex
          v2 = xyz(n23(2),:)-p1; % Compute vector to second vertex
          ang1 = v1*v2'; % Compute dot product
          ang1 = ang1./sqrt((v1*v1')*(v2*v2')); % Normalize dot product
          ang1 = acos(ang1); % Compute angle
%
          n2 = setdiff(tri(idt(l),:),n23); % Get opposite vertex of second triangle
          p2 = xyz(n2,:); % Get coordinates of opposite vertex
          v3 = xyz(n23(1),:)-p2; % Compute vector to first vertex
          v4 = xyz(n23(2),:)-p2; % Compute vector to second vertex
          ang2 = v3*v4'; % Compute dot product
          ang2 = ang2./sqrt((v3*v3')*(v4*v4')); % Normalize dot product
          ang2 = acos(ang2); % Compute angle
          m = l/2; % Compute score index
          sc(m) = ang1+ang2-pi; % Update score for edge
       end
     end
%
% Get New Sort Order for Scores
%
     [ssc,iss] = sort(sc,1,'descend'); % Re-sort scores
end
%
return % Exit the function