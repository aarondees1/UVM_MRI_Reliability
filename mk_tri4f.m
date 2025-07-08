function [tri,nt,slx] = mk_tri4f(dat,ang,tol,iplt) % Define function to create triangular mesh from MRI femur slice data
%MK_TRI4F Makes a triangular mesh by using the ordered slice data from
%         the digitized MRI femur data.
%
%         [TRI,NT] = MK_TRI4F(DAT) given a cell array containing three
%         (3) columns matrices with slice coordinate point data, DAT,
%         returns the three (3) column triangle connectivity matrix,
%         TRI.  The number of returned triangles, NT, may also be
%         returned.
%
%         TRI = MK_TRI4F(DAT,ANG) The data is converted to cylindrical
%         coordinates with the Y-axis forming the Z-axis of the cylinder
%         (X-Z plane transformed to polar coordinates). The mesh is only
%         defined between the angles of the smaller of the two slices
%         plus an additional angle of ANG degrees.  By default, ANG is
%         three (3) degrees.
%
%         NOTES:  1.  Each slice coordinate data matrix must correspond
%                 to one index into the cell array DAT.
%
%                 2.  The coordinates should be ordered in the same
%                 direction in every slice.  The dot product of the
%                 directions of adjacent slices are used to check the
%                 ordering direction and the ordering direction is
%                 reversed if the dot product is negative.
%
%                 3.  The arclength along each slice is used to
%                 determine the triangulation.  The data is converted
%                 to cylindrical coordinates with the Y-axis forming
%                 the Z-axis of the cylinder (X-Z plane transformed to
%                 polar coordinates).
% 
%                 4.  The M-files plane_fit.m and xzpl2pol.m must be in
%                 the current path or directory.
%
%         31-Jul-2014 * Mack Gardner-Morse
%
%#######################################################################
%
% Check for Inputs
%
if (nargin<4) % Check if plotting flag is missing
  iplt = false; % Set default plotting flag to false
end
%
if iplt % Check if plotting is enabled
  figure; % Create new figure
  orient landscape; % Set figure orientation to landscape
end
%
if (nargin<3)||isempty(tol) % Check if tolerance is missing or empty
  tol = 0.1; % Set default tolerance for dot product
end
%
if (nargin<2)||isempty(ang) % Check if angle is missing or empty
  ang = 3; % Set default angle extension
end
%
if (nargin<1) % Check if data input is missing
  error(' *** ERROR in MK_TRI4F:  No input data!'); % Throw error for missing input
end
%
% Convert to Polar Coordinates
%
dat = dat(:); % Ensure data is a column cell array
dat = xzpl2pol(dat); % Convert to cylindrical coordinates
%
% Get Arc Lengths (Angles)
%

nslice = size(dat,1); % Get number of slices
slen = cell(nslice,1); % Initialize cell array for arc lengths
npts = zeros(nslice,1); % Initialize array for point counts
vec1 = zeros(nslice,3); % Initialize array for direction vectors
cntr = zeros(nslice,3); % Initialize array for plane centers
nvec = zeros(3,nslice); % Initialize array for normal vectors
angs = zeros(nslice,2); % Initialize array for angle ranges
irev = false; % Initialize reverse flag
for k = 1:nslice % Loop through slices
   xyz = dat{k}; % Get slice data
   xyz = [xyz(:,1)*180/pi xyz(:,3), xyz(:,2)]; % Convert to degrees and reorder
   angs(k,:) = [min(xyz(:,1)) max(xyz(:,1))]; % Store min and max angles
%
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3)); % Fit plane to slice
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
         warning([' *** WARNING in mk_tri4:  Ordering of points', ...
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
   npts(k) = size(xyz,1); % Store number of points
   if irev % Check if slice was reversed
     slen{k} = flipud(xyz(:,1)); % Store reversed angles
   else
     slen{k} = xyz(:,1); % Store angles
   end
end
%
n = [0; cumsum(npts)]; % Compute cumulative point indices
tri = []; % Initialize triangle connectivity matrix
slx = zeros(nslice-1,1); % Initialize slice separation array
%
% Get Triangulation between Slices
%
for k = 2:nslice % Loop through slice pairs
%
% Limit Angles
%
   angl = [max(angs(k-1:k,1))-ang min(angs(k-1:k,2))+ang]; % Compute angle limits
   idx1 = find(slen{k-1}>angl(1)&slen{k-1}<angl(2)); % Find points within angle limits for first slice
   idx2 = find(slen{k}>angl(1)&slen{k}<angl(2)); % Find points within angle limits for second slice
%
   npt1 = size(idx1,1); % Get number of points in first slice
   npt2 = size(idx2,1); % Get number of points in second slice
%
% Slice Separation
%
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1); % Compute slice separation
%
% Delaunay Triangulation
%
   xt = [zeros(npt1,1); slx(k-1)*ones(npt2,1)]; % Create X coordinates for triangulation
   yt = [slen{k-1}(idx1); slen{k}(idx2)]; % Create Y coordinates for triangulation
   tril = delaunay(xt,yt); % Perform Delaunay triangulation
   nid1 = n(k-1)+1:n(k); % Get point indices for first slice
   nid1 = nid1(idx1); % Filter indices by angle limits
   nid2 = n(k)+1:n(k+1); % Get point indices for second slice
   nid2 = nid2(idx2); % Filter indices by angle limits
   nid = [nid1 nid2]; % Combine indices
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