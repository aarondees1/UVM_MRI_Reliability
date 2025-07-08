function [r,xyz1,xyz2,sse,res,exitflag] = cyl_fit(ri,xyz1i,xyz2i,xyzp,tol)
%CYL_FIT  Fits a cylinder to a set of points using nonlinear least
%         squares.
%
%         [R,XYZ1,XYZ2] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP) given an initial
%         radius RI, a row vector of the coordinates of an initial
%         point on the cylinder axis XYZ1I, an initial second point on
%         the cylinder axis XYZ2I, and a matrix of points XYZP with the
%         X, Y, and and Z coordinates in columns, returns the nonlinear
%         least squares fit of a cylinder with radius R, a point on the
%         axis of the cylinder XYZ1, and a second point on the axis of
%         the cylinder XYZ2.
%
%         [R,XYZ1,XYZ2] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP,TOL) set the
%         parameter and function tolerances in the optimization to TOL.
%
%         [R,XYZ1,XYZ2,SSE,RES] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP) returns
%         the sum of squared errors SSE and the residuals RES.
%
%         [R,XYZ1,XYZ2,SSE,RES,EXITFLAG] = CYL_FIT(RI,XYZ1I,XYZ2I,XYZP)
%         returns the EXITFLAG from MATLAB function lsqnonlin or
%         fminsearch.  Use "help lsqnonlin" or "help fminsearch" to see
%         the exit conditions for the different values of EXITFLAG.
%
%         NOTES:  1.  See cyl_plt.m to plot the resulting cylinder.
%
%                 2.  The M-file dist2cyl.m and pts2lin.m must be in
%                 the current path or directory.
%
%         17-Jul-2013 * Mack Gardner-Morse
%
%         21-Nov-2014 * Mack Gardner-Morse * Added the METHOD option.
%
%         25-Nov-2014 * Mack Gardner-Morse * Check cylinder axis matches
%                                            input axis direction.
%
%         10-Aug-2022 * Mack Gardner-Morse * Removed the METHOD option.
%         Only uses the Matlab nonlinear least squares function
%         lsqnonlin.
%

% Check if at least four input arguments are provided.
if (nargin<4)
  error(' *** ERROR in CYL_FIT:  Four (4) inputs are required!');
end

% Set default tolerance if not provided or empty.
if (nargin<5)||isempty(tol)
  tol = eps^(2/3); % Use machine epsilon-based tolerance
end

% Validate input dimensions: xyz1i and xyz2i must be 1x3, xyzp must have 3 columns.
[nrow1,ncol1] = size(xyz1i);
[nrow2,ncol2] = size(xyz2i);
[~,ncol3] = size(xyzp);
if nrow1~=1||nrow2~=1||ncol1~=3||ncol2~=3||ncol3~=3
  error([' *** ERROR in CYL_FIT:  All input points must have', ...
         ' three columns!']);
end

% Initialize parameter vector for optimization.
xi = zeros(7,1); % 7 parameters: radius (1), xyz1 (2:4), xyz2 (5:7)
if ri>0
  xi(1) = ri; % Set initial radius
else
  xi(1) = 1; % Default radius if invalid
end
xi(2:4) = xyz1i'; % Initial point 1 on cylinder axis
xi(5:7) = xyz2i'; % Initial point 2 on cylinder axis
zci = xyz2i-xyz1i; % Initial cylinder axis direction
zci = zci./sqrt(zci*zci'); % Normalize to unit direction vector

% Configure optimization settings for lsqnonlin.
opt = optimset('lsqnonlin');
opt = optimset(opt,'Display','off','Algorithm','levenberg-marquardt');
opt.MaxFunEvals = 5e+7; % Maximum function evaluations
opt.MaxIter = 2e+7;     % Maximum iterations
opt.TolFun = tol;       % Function tolerance
opt.TolX = tol;         % Parameter tolerance

% Perform nonlinear least squares fit to minimize distances to cylinder.
lb = -Inf(7,1); % Lower bounds for parameters
lb(1) = 0;      % Radius must be positive
[x,sse,res,exitflag] = lsqnonlin(@(x) dist2cyl(x,xyzp),xi,lb,[],opt);
% x: optimized parameters, sse: sum of squared errors, res: residuals

% Extract output parameters from optimization result.
r = x(1);       % Optimized radius
xyz1 = x(2:4)'; % Optimized point 1 on cylinder axis
xyz2 = x(5:7)'; % Optimized point 2 on cylinder axis

% Ensure cylinder axis direction matches input direction.
zc = xyz2-xyz1; % Optimized cylinder axis direction
zc = zc./sqrt(zc*zc'); % Normalize to unit direction vector
cdir = zci*zc'; % Dot product to check direction alignment
if cdir<0
  xyzt = xyz1; % Swap xyz1 and xyz2 if direction is opposite
  xyz1 = xyz2;
  xyz2 = xyzt;
  clear xyzt; % Clear temporary variable
  zc = -zc;   % Reverse direction vector
end

% Adjust cylinder ends to fit the range of input data points.
[xyzl,t] = pts2lin(xyz1,zc,xyzp); % Project points onto axis, get distances
[~,idmn] = min(t); % Index of minimum distance along axis
[~,idmx] = max(t); % Index of maximum distance along axis
xyz1 = xyzl(idmn,:); % Set first end of cylinder
xyz2 = xyzl(idmx,:); % Set second end of cylinder

return % Exit function