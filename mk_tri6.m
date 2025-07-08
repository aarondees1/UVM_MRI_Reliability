function [tri,nt,slx] = mk_tri6(dat,tol,iplt); % Define function to create triangular mesh from MRI slice data
%MK_TRI6 Makes a triangular mesh by using the ordered slice data from
%        the digitized MRI.  The normals of the triangles are oriented
%        so the Z component is in the positive Z direction.
%
%        [TRI,NT] = MK_TRI6(DAT) given a cell array containing three (3)
%        columns matrices with slice coordinate point data, DAT, returns
%        the three (3) column triangle connectivity matrix, TRI.  The
%        number of returned triangles, NT, may also be returned.
%
%        NOTES:  1.  Each slice coordinate data matrix must correspond
%                to one index into the cell array DAT.
%
%                2.  The coordinates should be ordered in the same
%                direction in every slice.  The dot product of the
%                directions of adjacent slices are used to check the
%                ordering direction and the ordering direction is
%                reversed if the dot product is negative.
%
%                3.  The arclength along each slice is used to determine
%                the triangulation.  See mk_tris.m for a triangulation
%                based on each slice starting at the same in slice
%                point.  See mk_triq.m for a quicker, but different
%                triangulation based on the number of points in each
%                slice.
%
%                4.  The M-file plane_fit.m and tri_norm.m must be in
%                the current path or directory.
%
%        27-Aug-2013 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<3) % Check if plotting flag is missing
  iplt = false; % Set default plotting flag to false
end
%
if iplt % Check if plotting is enabled
  hf = figure; % Create new figure
  orient tall; % Set figure orientation to tall
end
%
if (nargin<2) % Check if tolerance is missing
  tol = 0.1; % Set default tolerance for dot product
end
%
if (nargin<1) % Check if data input is missing
  error(' *** ERROR in MK_TRI6:  No input data!'); % Throw error for missing input
end
%
% Get Arc Lengths
%
dat = dat(:); % Ensure data is a column cell array
nslice = size(dat,1); % Get number of slices
slen = cell(nslice,1); % Initialize cell array for arc lengths
npts = zeros(nslice,1); % Initialize array for point counts
rpt1 = zeros(nslice,3); % Initialize array for first points
vec1 = zeros(nslice,3); % Initialize array for direction vectors
cntr = zeros(nslice,3); % Initialize array for plane centers
nvec = zeros(3,nslice); % Initialize array for normal vectors
irev = false; % Initialize reverse flag
for k = 1:nslice % Loop through slices
   xyz = dat{k}; % Get slice data
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3)); % Fit plane to slice
   [mxv,idmx] = max(abs(nvec(:,k))); % Find maximum normal component
   if nvec(idmx,k)<0 % Check if normal is negative
     nvec(:,k) = -nvec(:,k); % Flip normal direction
   end
   vec = xyz(end,:)-xyz(1,:); % Compute direction vector
   vec1(k,:) = vec./norm(vec); % Normalize direction vector
%
% Check for Slices with a Reverse Digitization
%
   if k>1 % Check if not first slice
     dotp = vec1(k-1,:)*vec1(k,:)' ; % Compute dot product with previous slice
     if dotp<tol % Check if directions differ significantly
       irev = true; % Set reverse flag
       xyz = flipud(xyz); % Reverse point order
       vec = xyz(end,:)-xyz(1,:); % Recompute direction vector
       vec1(k,:) = vec./norm(vec); % Normalize new direction vector
       dotp2 = vec1(k-1,:)*vec1(k,:)' ; % Compute new dot product
       if dotp2<dotp % Check if reversal worsened alignment
         warning([' *** WARNING in mk_tri6:  Ordering of points', ...
                  ' in the slices may not be in the same direction!']); % Issue warning
         irev = false; % Revert reverse flag
         xyz = flipud(xyz); % Restore original order
         vec = xyz(end,:)-xyz(1,:); % Recompute direction vector
         vec1(k,:) = vec./norm(vec); % Normalize direction vector
       end
     else
       irev = false; % Clear reverse flag
     end
   end
   rpt1(k,:) = xyz(1,:); % Store first point
   npts(k) = size(xyz,1); % Store number of points
   dd = diff(xyz); % Compute point differences
   dlen = sqrt(sum(dd.*dd,2)); % Compute segment lengths
   if irev % Check if slice was reversed
     slen{k} = flipud([0; cumsum(dlen)]); % Store reversed arc lengths
   else
     slen{k} = [0; cumsum(dlen)]; % Store arc lengths
   end
end
%
n = [0; cumsum(npts)]; % Compute cumulative point indices
tri = []; % Initialize triangle connectivity matrix
slx = zeros(nslice-1,1); % Initialize slice separation array
for k = 2:nslice % Loop through slice pairs
%
% Slice Separations and Offsets
%
   ds = rpt1(k,:)-rpt1(k-1,:); % Compute difference between first points
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1); % Compute slice separation
   offst = ds*vec1(k-1,:)' ; % Compute offset along direction
%
% Delaunay Triangulation
%
   xt = [zeros(npts(k-1),1); slx(k-1)*ones(npts(k),1)]; % Create X coordinates for triangulation
   yt = [slen{k-1}-offst; slen{k}]; % Create Y coordinates with offset
 %  xt = [zeros(npts(k-1),1); ones(npts(k),1)]; % Commented-out alternative X coordinates
 %  yt = [(0:1/(npts(k-1)-1):1)'; (0:1/(npts(k)-1):1)']; % Commented-out alternative Y coordinates
   tril = delaunay(xt,yt); % Perform Delaunay triangulation
   [nx,ny,nz] = tri_norm(tril,[xt yt zeros(npts(k-1)+npts(k),1)]); % Compute triangle normals
   idn = find(nz<0); % Find triangles with negative Z normals
   if ~isempty(idn) % Check if negative normals exist
     tril(idn,:) = tril(idn,[1 3 2]); % Reorient triangles
   end
   nid = n(k-1)+1:n(k+1); % Get point indices for slices
   tri = [tri; nid(tril)]; % Append triangles to connectivity matrix
   if iplt % Check if plotting is enabled
     ntril = size(tril,1); % Get number of triangles
     cla; % Clear axes
     plot(xt,yt,'k.'); % Plot points
     hold on; % Enable hold
     trimesh(tril,xt,yt); % Plot triangulation
     text(xt,yt,int2str((1:length(xt))'),'Color','b','FontSize',12); % Label points
     text(mean(xt(tril),2),mean(yt(tril),2),int2str((1:ntril)'), ...
          'Color','r','FontSize',12); % Label triangles
     pause; % Pause for inspection
   end
%
end
%
nt = size(tri,1); % Compute number of triangles
%
return % Exit the function