function [tri,nt,slx] = mk_tri4f(dat,ang,tol,iplt)
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
if (nargin<4)
  iplt = false;
end
%
if iplt
  figure;
  orient landscape;
end
%
if (nargin<3)||isempty(tol)
  tol = 0.1;            % Tolerance on dot product (projection of directions)
end
%
if (nargin<2)||isempty(ang)
  ang = 3;
end
%
if (nargin<1)
  error(' *** ERROR in MK_TRI4F:  No input data!');
end
%
% Convert to Polar Coordinates
%
dat = dat(:);
dat = xzpl2pol(dat);
%
% Get Arc Lengths (Angles)
%
nslice = size(dat,1);
slen = cell(nslice,1);
npts = zeros(nslice,1);
% rpt1 = zeros(nslice,3);
vec1 = zeros(nslice,3);
cntr = zeros(nslice,3);
nvec = zeros(3,nslice);
angs = zeros(nslice,2);
irev = false;
for k = 1:nslice
   xyz = dat{k};
   xyz = [xyz(:,1)*180/pi xyz(:,3), xyz(:,2)];
   angs(k,:) = [min(xyz(:,1)) max(xyz(:,1))];
%
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
   vec = xyz(end,:)-xyz(1,:);
   vec1(k,:) = vec./norm(vec);
%
% Check for Slices with a Reverse Digitization
%
   if k>1
     dotp = vec1(k-1,:)*vec1(k,:)';
     if dotp<tol
       irev = true;
       xyz = flipud(xyz);
       vec = xyz(end,:)-xyz(1,:);
       vec1(k,:) = vec./norm(vec);
       dotp2 = vec1(k-1,:)*vec1(k,:)';
       if dotp2<dotp    % Revert back to original ordering
         warning([' *** WARNING in mk_tri4:  Ordering of points', ...
                  ' in the slices may not be in the same direction!']);
         irev = false;
         xyz = flipud(xyz);
         vec = xyz(end,:)-xyz(1,:);
         vec1(k,:) = vec./norm(vec);
       end
     else
       irev = false;
     end
   end
   npts(k) = size(xyz,1);
   if irev
     slen{k} = flipud(xyz(:,1));
   else
     slen{k} = xyz(:,1);
   end
end
%
n = [0; cumsum(npts)];
tri = [];
slx = zeros(nslice-1,1);               % Slice separation
%
% Get Triangulation between Slices
%
for k = 2:nslice
%
% Limit Angles
%
   angl = [max(angs(k-1:k,1))-ang min(angs(k-1:k,2))+ang];
   idx1 = find(slen{k-1}>angl(1)&slen{k-1}<angl(2));
   idx2 = find(slen{k}>angl(1)&slen{k}<angl(2));
%
   npt1 = size(idx1,1);
   npt2 = size(idx2,1);
%
% Slice Separation
%
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1);
%
% Delaunay Triangulation
%
   xt = [zeros(npt1,1); slx(k-1)*ones(npt2,1)];
   yt = [slen{k-1}(idx1); slen{k}(idx2)];
   tril = delaunay(xt,yt);
   nid1 = n(k-1)+1:n(k);
   nid1 = nid1(idx1);
   nid2 = n(k)+1:n(k+1);
   nid2 = nid2(idx2);
   nid = [nid1 nid2];
   tri = [tri; nid(tril)];
   if iplt
     ntril = size(tril,1);
     cla;
     plot(xt,yt,'k.');
     hold on;
     trimesh(tril,xt,yt);
     text(xt,yt,int2str((1:length(xt))'),'Color','b','FontSize',12);
     text(mean(xt(tril),2),mean(yt(tril),2),int2str((1:ntril)'), ...
          'Color','r','FontSize',12);
     pause;
   end
%
end
%
nt = size(tri,1);
%
return