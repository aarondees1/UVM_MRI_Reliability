%#######################################################################
%
%               * Femoral Cartilage Thickness Calculation Program *
%                         (Script Name: fcart08_sag_AD)
%
%          This M-File processes femoral bone and cartilage data for
%     selected subjects, establishing a 1 mm by 2 degree grid. It scales
%     and projects this grid onto the bone surface mesh, then calculates
%     cartilage thicknesses by finding the intersection of bone normals
%     with the cartilage surface mesh. 
%
%     The program requires pre-existing bone and cartilage MAT files,
%     which are outputs from 'femurs08b_AD' and 'femurs08c_AD'
%     respectively.
%
%     NOTES:  1.  The M-files car_thk8.m, gridproj.m, line_fit.m,
%             nod2tri.m, nod_norm.m, pts2lin.m, tri_fix2.m,
%             tri_norm.m, tsect4.m and xprod.m must be in the current
%             path or directory.
%
%             2.  Bone and cartilage MAT files must already exist in
%             the selected subject directories.
%
%             3.  This program produces a MAT file containing the grid
%             scaling for femurs in the format ***_fgrid08_.mat
%             for a given scan type (e.g., FFE_fgrid08_.mat).
%             Final thicknesses are saved in a MAT file in the format
%             ***_fcart08_thk.mat for a given subject, knee (left or
%             right) and scan type (e.g., 001_L_FFE_fcart08_thk.mat).
%            
%             4.  This M-file outputs PDF files with plots of the
%             femoral bony and cartilage surfaces, and femoral cartilage
%             thicknesses.
%
%             5. The code processes all three scan types sequentially.
%
%     7-July-2025 * Mack Gardner-Morse & Aaron Dees
%
%#######################################################################

clear; % Clear all variables from workspace
close all; % Close all open figures
clc; % Clear command window
iprt =true; % Enable printing flag
csnam = '_femurCS.mat'; % Define suffix for femur coordinate system MAT file
cnam = '_femurCart.mat'; % Define suffix for femur cartilage MAT file
%
%
% Get All Subject Subdirectories
%
div = uigetdir; % Open dialog to select base directory
ddir = fullfile(div); % Construct full path to selected directory
rdir = fullfile(div,'Results'); % Construct path to results directory
bdir = fullfile(rdir,'Bone'); % Construct path to bone data directory
cdir = fullfile(rdir, 'Cartilage'); % Construct path to cartilage data directory
tdir = fullfile(rdir,'Thickness'); % Construct path to thickness data directory
gdir = fullfile(rdir,'Grids'); % Construct path to grids directory
load(fullfile(rdir, 'All Subjects Full Tibia and Femur Bone and Cartilage Data.mat')); % Load subject data from MAT file
ns=size(sd,1); % Get number of subjects
%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1); % Initialize array for leg indicators (false for left)
kids = repmat(blanks(ns)',1,5); % Initialize array for knee IDs
tribfs = cell(ns,1); % Initialize cell array for bone triangular mesh connectivity
xyzbfs = cell(ns,1); % Initialize cell array for bone point coordinates
xyzmnl = zeros(ns,3); % Initialize array for minimum lateral femur coordinates
xyzmxl = zeros(ns,3); % Initialize array for maximum lateral femur coordinates
xyzmnm = zeros(ns,3); % Initialize array for minimum medial femur coordinates
xyzmxm = zeros(ns,3); % Initialize array for maximum medial femur coordinates
xyzs = zeros(ns,3); % Initialize array for proximal femur outline range
rf = zeros(ns,3); % Initialize array for femur condyle radii
scan(1,:)='FFE'; % Set first scan type to FFE
scan(2,:)='RHO'; % Set second scan type to RHO
scan(3,:)='T2S'; % Set third scan type to T2S
%
% Initialize Overall Minimums and Maximums
%
hpi = pi/2; % Define half pi constant
dpi = 2*pi; % Define double pi constant
rad2deg = 180/pi; % Define radians to degrees conversion factor
%
% Loop through MAT Files
%
%% Load Bone Data and Generate min/max and range of Grid
tpw=zeros(ns,3); % Initialize array for tibial plateau widths
zsc=zeros(ns,3); % Initialize array for Z-scaling factors
for j=1:3 % Loop through scan types (FFE, RHO, T2S)
    tzrmn = zeros(1,3); % Initialize minimums for cylindrical coordinates
    tzrmx = zeros(1,3); % Initialize maximums for cylindrical coordinates
    for i = 1:ns % Loop through each subject
        i_str=int2str(i); % Convert subject index to string
        if j==1 % Check for FFE scan type
            fstr=sd(i).FFE.femur.bfnam; % Get FFE bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        elseif j==2 % Check for RHO scan type
            fstr=sd(i).RHO.femur.bfnam; % Get RHO bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        elseif j==3 % Check for T2S scan type
            fstr=sd(i).T2S.femur.bfnam; % Get T2S bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        end
        bone = load(fullfile(rdir,'Bone',[fstr csnam])); % Load bone data from MAT file
        cart = load(fullfile(rdir,'Cartilage',[fstr cnam])); % Load cartilage data from MAT file
        trifb = mk_tri4f(bone.datfb); % Create triangular mesh for bone data
        xyzfb = cell2mat(bone.datfb); % Convert bone data to matrix
        trifb = tri_fix2(trifb,xyzfb); % Fix triangulation errors in bone mesh
        datfc = comb_dat(cart.datlct,cart.datmct,cart.dattct); % Combine cartilage data from all regions
        trifc = mk_tri4f(datfc); % Create triangular mesh for cartilage data
        xyzfc = cell2mat(datfc); % Convert cartilage data to matrix
        trifc = tri_fix2(trifc,xyzfc); % Fix triangulation errors in cartilage mesh
        eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.trimesh = trifb;']); % Store bone trimesh in subject data
        eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.points = xyzfb;']); % Store bone points in subject data
        eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.trimesh = trifc;']); % Store cartilage trimesh in subject data
        eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.points = xyzfc;']); % Store cartilage points in subject data
        eval(['xyzs(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline)-min(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline);']); % Calculate tibia outline range
        % Transform to Cylindrical Coordinate System
        %
        [th,r,z] = cart2pol(xyzfb(:,1),xyzfb(:,3),xyzfb(:,2)); % Convert bone points to cylindrical coordinates
        id = find(th>hpi); % Find angles greater than pi/2
        th(id) = th(id)-dpi; % Adjust angles by subtracting 2*pi
        tzr = [th,z,r]; % Combine cylindrical coordinates
        %
        % Get Minimums and Maximums
        %
        tzrn = min(tzr); % Calculate minimum cylindrical coordinates
        tzrx = max(tzr); % Calculate maximum cylindrical coordinates
        itst = tzrn<tzrmn; % Identify where subject minimums are less than global minimums
        tzrmn(itst) = tzrn(itst); % Update global minimums
        itst = tzrx>tzrmx; % Identify where subject maximums exceed global maximums
        tzrmx(itst) = tzrx(itst); % Update global maximums
    end
    tpw(:,j) = xyzs(:,2); % Store Y-range of tibia outline for scan type
    tpw_mn = mean(tpw(:,j)); % Calculate mean tibial plateau width
    zsc(:,j) = tpw(:,j)./tpw_mn; % Compute Z-scaling factors
    %
    % Get Range and Make a Rectangular Grid
    %
    tzrmn(1) = tzrmn(1)*180/pi; % Convert minimum theta to degrees
    tzrmx(1) = tzrmx(1)*180/pi; % Convert maximum theta to degrees
    tzmn = floor(tzrmn(1:2)); % Get minimum theta and Z rounded down
    tzmx = ceil(tzrmx(1:2)); % Get maximum theta and Z rounded up
    tzmn(1) = floor(tzmn(1)/2)*2; % Adjust minimum theta to even degree
    tzmx(1) = ceil(tzmx(1)/2)*2; % Adjust maximum theta to even degree
    t = (tzmn(1):2:tzmx(1))'; % Create theta range in 2-degree increments
    nt = size(t,1); % Get number of theta points
    z = (tzmn(2):tzmx(2)+2)'; % Create Z range in 1-mm increments
    nz = size(z,1); % Get number of Z points
    [T,Z] = meshgrid(t,z); % Create rectangular grid of theta and Z
    tq = T(:); % Flatten theta grid to vector
    zq = Z(:); % Flatten Z grid to vector
    nq = size(tq,1); % Get number of grid points
    %
    % FROM QUAD CONNECT
    quad = (1:nz-1)'; % Initialize quadrilateral connectivity for first row
    quad = [quad quad+nz quad+nz+1 quad+1]; % Define quadrilateral connectivity
    quad = repmat(quad,nt-1,1); % Replicate quadrilaterals for all theta
    addcol = repmat(nz*(0:nt-2),(nz-1)*4,1); % Create offset for quadrilateral indices
    addcol = reshape(addcol,4,(nt-1)*(nz-1))'; % Reshape offset to match quadrilaterals
    quad = quad+addcol; % Adjust quadrilateral indices
    %
    % Triangular Grid
    %
    ne = size(quad,1); % Get number of quadrilateral elements
    trig = [quad(:,1:3) quad(:,1) quad(:,3:4)]'; % Convert quadrilaterals to triangles
    trig = reshape(trig,3,2*ne)'; % Reshape triangle connectivity
    for i=1:ns % Loop through subjects to save grid file names
        i_str=int2str(i); % Convert subject index to string
        gnam = [scan(j,:) '_fgrid08_.mat']; % Construct grid MAT file name
        eval(['sd(' i_str ').' scan(j,:) '.femur.grid = gnam;']); % Store grid file name in subject data
        gnam = fullfile(gdir,gnam); % Construct full path to grid MAT file
    end
    %
    % Save Data
    %
    save(gnam, 'nq', 'nt', 'nz', 'ne', 'trig','quad', 'tq', 'zq', 'zsc', 'tzrmn', 'tzrmx'); % Save grid data to MAT file
    clear quad addcol nq tq nz T Z nt t tzrmn tzrmx tzrn tzmn tzmx tzrx xyzs; % Clear temporary variables
end
for j=1:3 % Loop through scan types (FFE, RHO, T2S)
    for i = 1:ns % Loop through each subject
        i_str=int2str(i); % Convert subject index to string
        if j==1 % Check for FFE scan type
            fstr=sd(i).FFE.femur.bfnam; % Get FFE bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        elseif j==2 % Check for RHO scan type
            fstr=sd(i).RHO.femur.bfnam; % Get RHO bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        elseif j==3 % Check for T2S scan type
            fstr=sd(i).T2S.femur.bfnam; % Get T2S bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
        end
        gnam = [scan(j,:) '_fgrid08_.mat']; % Construct grid MAT file name
        load(fullfile(gdir,gnam)); % Load grid data from MAT file
        cthk = zeros(nq,1); % Initialize array for cartilage thicknesses
        xyzic = zeros(nq,3,1); % Initialize array for cartilage intersection points
        rq = zeros(nq,1); % Initialize array for radial coordinates
        zqs = zeros(nq,1); % Initialize array for scaled Z coordinates
        eval(['trifb = sd(' i_str ').' scan(j,:) '.femur.bone.full.trimesh;']); % Extract bone trimesh from subject data
        eval(['xyzfb = sd(' i_str ').' scan(j,:) '.femur.bone.full.points;']); % Extract bone points from subject data
        eval(['trifc = sd(' i_str ').' scan(j,:) '.femur.cart.full.trimesh;']); % Extract cartilage trimesh from subject data
        eval(['xyzfc = sd(' i_str ').' scan(j,:) '.femur.cart.full.points;']); % Extract cartilage points from subject data
        %
        % Get Cylindrical Coordinates
        %
        [t,r,z] = cart2pol(xyzfb(:,1),xyzfb(:,3),xyzfb(:,2)); % Convert bone points to cylindrical coordinates
        id = find(t>pi/2); % Find angles greater than pi/2
        t(id) = t(id)-2*pi; % Adjust angles by subtracting 2*pi
        td = t*rad2deg; % Convert angles to degrees
        %
        % Scale Polar "Z" Coordinates
        %
        zs = zsc(i,j)*zq; % Scale Z coordinates for grid
        %
        % Interpolate to Scaled "Standard" Grid
        %
        dist=12; % Set distance parameter for projection
        rg = gridproj(trifb,[td z r],tq,zs,1,dist); % Project bone mesh onto grid
        %
        % Back to Cartesian Coordinates
        %
        [xg,zg,yg] = pol2cart(tq/rad2deg,rg,zs); % Convert grid back to Cartesian coordinates
        xyzg = [xg yg zg]; % Combine Cartesian coordinates
        %
        % Improve the Mesh
        %
        trig = tri_fix2(trig,xyzg); % Fix triangulation errors in grid mesh
        %
        % Calculate Cartilage Thicknesses
        %
        [ct,bi] = car_thk8(trifc,xyzfc,trig,xyzg); % Calculate cartilage thicknesses
        %
        % Save Thicknesses and Print Plot
        %
        cthk(:,1) = ct'; % Store cartilage thicknesses
        xyzic(:,:,1) = bi; % Store cartilage intersection points
        % Plot the Cartilage Thicknesses
        %
        hf2 = figure; % Create new figure for cartilage and bone meshes
        hb2 = trimesh(trig,xyzg(:,1),xyzg(:,2),xyzg(:,3),'LineWidth',0.25,'FaceColor','none','EdgeColor','k'); % Plot grid mesh
        hold on; % Enable hold for multiple plots
        hc2 = trimesh(trifc,xyzfc(:,1),xyzfc(:,2),xyzfc(:,3),'LineWidth',0.25,'FaceColor','none','EdgeColor','b'); % Plot cartilage mesh
        xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
        ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
        zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
        title({fstr; 'Sagittal Cartilage and Bone Meshes'},'Interpreter','none','FontSize',16,'FontWeight','bold'); % Add title with subject ID
        view(-60,12); % Set view angle
        axis equal; % Set equal scaling for axes
        psnam = fullfile(tdir,[fstr '_fcart08_sagittal_thk.pdf']); % Construct PDF file path
        orient landscape; % Set figure orientation to landscape
        set(hf2, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
        exportgraphics(hf2, psnam, "Resolution", 300); % Export mesh plot to PDF
        close(hf2); % Close mesh figure
        % Plot Sagittal Cartilage Thicknesses
        %
        hf4 = figure; % Create new figure for cartilage thicknesses
        hb4 = trimesh(trig,xyzg(:,1),xyzg(:,2),xyzg(:,3),'LineWidth',0.25,'FaceColor','none','EdgeColor','b'); % Plot grid mesh
        hold on; % Enable hold for multiple plots
        ctf = squeeze(cthk(:,1)); % Extract cartilage thickness values
        idx = find(~isnan(ctf)); % Find non-NaN thickness indices
        it = nod2tri(idx,trig,2); % Convert nodes to triangle indices
        hs4l = trisurf(trig(it,:),xyzg(:,1),xyzg(:,2),xyzg(:,3),ctf,'FaceColor','interp','EdgeColor','b','LineWidth',0.25); % Plot thickness surface
        colormap jet; % Set colormap to jet
        view(-60,12); % Set view angle
        axis equal; % Set equal scaling for axes
        if min(ctf)<0 % Check if minimum thickness is negative
            caxis([min(ctf) max(ctf)]); % Set color axis to full range
        else
            caxis([0 max(ctf)]); % Set color axis from zero to maximum
        end
        hcb4 = colorbar; % Add colorbar
        xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
        ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
        zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
        title({fstr; 'Sagittal Cartilage Thicknesses'},'Interpreter','none','FontSize',16,'FontWeight','bold'); % Add title with subject ID
        psnam = fullfile(tdir,[fstr '_fcart08_sagittal_thk.pdf']); % Construct PDF file path
        orient landscape; % Set figure orientation to landscape
        set(hf4, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
        exportgraphics(hf4, psnam, "Resolution", 300, 'Append',true); % Append thickness plot to PDF
        close(hf4);
        %
        % Save Tibia Data
        %
        kid=fstr(1:5); % Extract knee identifier
        ileg=false; % Set leg indicator to false (left knee)
        rnam = fullfile(tdir,[fstr '_fcart08_thk.mat']); % Construct MAT file path
        save(rnam,'cthk', 'trig', 'trifc', 'xyzg', 'xyzfc'); % Save thickness data to MAT file
        %
    end
end
main_file=fullfile(rdir, 'All Subjects Tibia and Femur Cartilage Thickness Data.mat'); % Construct path to main MAT file
save (main_file, 'sd'); % Save subject data to main MAT file
return % Exit the script