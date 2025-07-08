%#######################################################################
%
%           * Tibial Cartilage Thickness Calculation, Gridding, and Visualization *
%                         (Script Name: tcart08_sag_AD)
%
%          This program guides the user to select a data directory,
%     then processes digitized tibia bone and cartilage data for
%     selected subjects. It is designed to run after 'tibias08c_AD',
%     utilizing the processed bone and cartilage data generated
%     by that script.
%     The program defines a 1mm by 1mm grid, scales and projects this
%     grid onto the bone surface mesh, and calculates cartilage
%     thicknesses by intersecting bone normals with the cartilage surface.
%
%     Results include calculated thicknesses, generated grids, and plots
%     of bone/cartilage surfaces and thickness maps. All plots are saved
%     as PDF files into the results directory.
%     The calculated cartilage thicknesses are saved in a MAT file
%     named '***_tcart08_thk.mat', and the generated grid data is saved
%     as 'scanType_tgrid08_.mat'. Finally, the 'sd' structure, containing 
%     updated subject information including file paths, is also saved for 
%     comprehensive data tracking and subsequent scripts.
%
%     NOTES:  1. The M-files car_thk8.m, gridproj.m, line_fit.m,
%             nod2tri.m, nod_norm.m, pts2lin.m, tri_fix2.m,
%             tri_norm.m, tsect4.m and xprod.m must be in the current
%             path or directory.
%
%             2. **This script requires the MAT files output by 'tibias08b_AD'
%             (for bone data and coordinate systems) and 'tibias08c_AD'
%             (for transformed cartilage data and meshes) to be present
%             in the selected subject directories.**
%
%             3. This program stores generated grid data in a MAT file
%             named 'scanType_tgrid08_.mat' and calculated cartilage
%             thicknesses in '***_tcart08_thk.mat' for a given subject,
%             knee (left or right) and scan type (e.g. 001_L_FFE_tcart08_thk.mat).
%             These are saved within the results directory of the knee data.
%
%             4. This M-file outputs various PDF files including
%             plots of the tibial bony and cartilage surfaces, as well
%             as cartilage thickness maps, all saved into the subjects'
%             results directory.
%
%     7-July-2025 * Mack Gardner-Morse & Aaron Dees
%
%#######################################################################
clear; % Clear all variables from the workspace
close all; % Close all open figures
clc; % Clear the command window
%
% Plot Parameter
%
iplt = true; % Enable plots
%
iprt = true; % Enable printing of plots

csnam = '_tibiaCS.mat'; % Suffix for tibia coordinate system MAT file
cnam = '_tibiaCart.mat'; % Suffix for tibia cartilage MAT file
%
% Get All Subject Subdirectories
%
%
%
div = uigetdir; % Prompt user to select a directory containing subject data
ddir = fullfile(div); % Construct full path to the selected directory
rdir = fullfile(div,'Results'); % Construct path to results directory
tdir = fullfile(rdir,'Thickness'); % Construct path to thickness results directory
gdir = fullfile(rdir,'Grids'); % Construct path to grids directory
if ~isfolder(tdir)
    mkdir(tdir); % Create thickness directory if it doesn't exist
end
if ~isfolder(gdir)
    mkdir(gdir); % Create grids directory if it doesn't exist
end
load(fullfile(rdir,'All Subjects Tibia Bone and Cartilage Data.mat')); % Load subject data from MAT file
ns=size(sd,1); % Number of subjects

%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1); % Array to store leg information (false for left, true for right)
kids = repmat(blanks(ns)',1,5); % Array to store knee IDs
tribl = cell(ns,1); % Cell array for sagittal bone lateral triangular mesh connectivity
tribm = cell(ns,1); % Cell array for sagittal bone medial triangular mesh connectivity
xyzbl = cell(ns,1); % Cell array for sagittal lateral bone point coordinates
xyzbm = cell(ns,1); % Cell array for sagittal medial bone point coordinates
xyz_min_l = zeros(ns,3); % Array for minimum individual lateral tibia coordinates
xyz_max_l = zeros(ns,3); % Array for maximum individual lateral tibia coordinates
xyz_min_m = zeros(ns,3); % Array for minimum individual medial tibia coordinates
xyz_max_m = zeros(ns,3); % Array for maximum individual medial tibia coordinates
xyzs = zeros(ns,3); % Array for range of proximal tibia outline
%

scan(1,:)='FFE'; % Define scan type: FFE
scan(2,:)='RHO'; % Define scan type: RHO
scan(3,:)='T2S'; % Define scan type: T2S
%% Load Bone Data and Generate min/max and range of Grid
for i = 1:ns
    i_str=int2str(i); % Convert subject index to string

    for j=1:3 % Loop through FFE, RHO, and T2S scans
        if j==1
            fstr=sd(i).FFE.tibia.bfnam; % FFE bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name

        elseif j==2
            fstr=sd(i).RHO.tibia.bfnam; % RHO bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name

        elseif j==3
            fstr=sd(i).T2S.tibia.bfnam; % T2S bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name
        end

        bone = load(fullfile(rdir,'Bone',[fstr csnam])); % Load bone data from MAT file

        if isfield(bone,'xyzpto') % Check if proximal tibia outline exists
            % Store lateral bone triangular mesh
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.trimesh = bone.tril;']);
            % Store medial bone triangular mesh
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.medial.trimesh = bone.trim;']);
            % Store lateral bone point coordinates
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points = bone.xyzl;']);
            % Store medial bone point coordinates
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points = bone.xyzm;']);
            % Store proximal tibia outline
            eval(['sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline = bone.xyzpto;']);

        else
            % Throw error if proximal tibia outline is missing
            error([' *** ERROR in tcart08:  Not able to find proximal ', ...
                'tibial outline.  Please run tibias08b again.']);
        end

        cart = load(fullfile(rdir,'Cartilage',[fstr cnam])); % Load cartilage data from MAT file

        % Store lateral cartilage triangular mesh
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.trimesh = cart.tril;']);
        % Store medial cartilage triangular mesh
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.medial.trimesh = cart.trim;']);
        % Store lateral cartilage point coordinates
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.points = cart.xyzl;']);
        % Store medial cartilage point coordinates
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.medial.points = cart.xyzm;']);
    end
end
for j=1:3 % Loop through scan types
    for i=1:ns % Loop through subjects
        i_str=int2str(i); % Convert subject index to string
        % Get Minimum and Maximum Coordinates
        % Lateral minimum coordinates
        eval(['xyz_min_l(' i_str ',:) = min(sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points);']);
        % Lateral maximum coordinates
        eval(['xyz_max_l(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points);']);
        % Medial minimum coordinates
        eval(['xyz_min_m(' i_str ',:) = min(sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points);']);
        % Medial maximum coordinates
        eval(['xyz_max_m(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points);']);

        % Get Range of Proximal Tibia Outline
        %
        % Calculate range of proximal tibia outline
        eval(['xyzs(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline)-min(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline);']);
        %
    end
    % Plot and Calculate Grid
    xyz_min_l=min(xyz_min_l); % Minimum lateral coordinates across all subjects
    xyz_max_l=max(xyz_max_l); % Maximum lateral coordinates across all subjects
    xyz_min_m=min(xyz_min_m); % Minimum medial coordinates across all subjects
    xyz_max_m=max(xyz_max_m); % Maximum medial coordinates across all subjects
    %
    % Generate Scaling (Based on Proximal Tibia Outline) and Uniform Grid
    %

    xyzs_avg = mean(xyzs); % Average range of proximal tibia outlines
    sc = xyzs./repmat(xyzs_avg,ns,1); % Scaling factors for each subject
    scx = sc(:,1); % X-axis scaling factors
    scy = sc(:,2); % Y-axis scaling factors
   
    %
    %
    % Get Range and Make a Uniform Rectangular Grid in X and Y for the
    % Lateral Compartment
    %

    x = (floor(xyz_min_l(1)):ceil(xyz_max_l(1)))'; % X range in 1 mm increments for lateral
    nxl = size(x,1); % Number of X grid points
    y = (floor(xyz_min_l(2)):ceil(xyz_max_l(2)))'; % Y range in 1 mm increments for lateral
    nyl = size(y,1); % Number of Y grid points
    [X,Y] = meshgrid(x,y); % Create rectangular grid of square quadrilaterals
    xql = X(:); % Grid X points as a vector
    yql = Y(:); % Grid Y points as a vector
    nl = size(xql,1); % Number of grid points in lateral compartment
    %
    quadl = (1:nyl-1)'; % Initialize quadrilateral connectivity for lateral
    quadl = [quadl quadl+nyl quadl+nyl+1 quadl+1]; % Define quadrilateral elements
    quadl = repmat(quadl,nxl-1,1); % Replicate for all X grid points
    addcol = repmat(nyl*(0:nxl-2),(nyl-1)*4,1); % Offset for quadrilateral indices
    addcol = reshape(addcol,4,(nxl-1)*(nyl-1))'; % Reshape offset matrix
    quadl = quadl+addcol; % Final quadrilateral connectivity
    %
    nel = size(quadl,1); % Number of lateral quadrilateral elements

    %
    % Get Range and Make a Uniform Rectangular Grid in X and Y for the
    % Medial Compartment
    %

    x = (floor(xyz_min_m(1)):ceil(xyz_max_m(1)))'; % X range in 1 mm increments for medial
    nxm = size(x,1); % Number of X grid points
    y = (floor(xyz_min_m(2)):ceil(xyz_max_m(2)))'; % Y range in 1 mm increments for medial
    nym = size(y,1); % Number of Y grid points
    [X,Y] = meshgrid(x,y); % Create rectangular grid of square quadrilaterals
    xqm = X(:); % Grid X points as a vector
    yqm = Y(:); % Grid Y points as a vector
    nm = size(xqm,1); % Number of grid points in medial compartment
    %
    quadm = (1:nym-1)'; % Initialize quadrilateral connectivity for medial
    quadm = [quadm quadm+nym quadm+nym+1 quadm+1]; % Define quadrilateral elements
    quadm = repmat(quadm,nxm-1,1); % Replicate for all X grid points
    addcol = repmat(nym*(0:nxm-2),(nym-1)*4,1); % Offset for quadrilateral indices
    addcol = reshape(addcol,4,(nxm-1)*(nym-1))'; % Reshape offset matrix
    quadm = quadm+addcol; % Final quadrilateral connectivity
    %
    nem = size(quadm,1); % Number of medial quadrilateral elements
    %
    % Triangular Grids
    %
    trigl = [quadl(:,1:3) quadl(:,1) quadl(:,3:4)]'; % Define triangular connectivity for lateral
    trigl = reshape(trigl,3,2*nel)'; % Reshape to triangular mesh format
    %
    trigm = [quadm(:,1:3) quadm(:,1) quadm(:,3:4)]'; % Define triangular connectivity for medial
    trigm = reshape(trigm,3,2*nem)'; % Reshape to triangular mesh format
    %
    % Save Analysis Grids
    %
    for i=1:ns
        i_str=int2str(i); % Convert subject index to string
        gnam = [scan(j,:) '_tgrid08_.mat']; % Construct grid MAT file name
        eval(['sd(' i_str ').' scan(j,:) '.tibia.grid = gnam;']); % Store grid file name in structure
        gnam = fullfile(gdir,gnam); % Construct full path for grid MAT file
    end
    %
    % if exist(rdir,'dir') % Commented out: check if results directory exists
    %   ierr = menu('Overwrite Previous Results','Yes','No')-1;
    %   if ierr
    %     return;
    %   end
    % else
    %   mkdir(rdir);
    % end
    %
    kid=fstr(1:5); % Extract knee identifier
    % Save grid data to MAT file
    save(gnam,'nel','nem','nl','nm','nxl','nxm','nyl','nym','kid', ...
        'quadl','quadm','sc','scx','scy','trigl','trigm','xql','xqm', ...
        'yql','yqm','xyz_min_l','xyz_max_l','xyz_min_m','xyz_max_m');

    %% Plot And Calculate Cartilage Thicknesses

    %
    % Initialize Cartilage Arrays
    %

    % Initialize Cartilage Thickness Arrays
    cthkl = zeros(nl,1); % Array for lateral sagittal cartilage thicknesses
    xyzil = zeros(nl,3,1); % Array for lateral sagittal cartilage intersection points
    %
    cthkm = zeros(nm,1); % Array for medial sagittal cartilage thicknesses
    xyzim = zeros(nm,3,1); % Array for medial sagittal cartilage intersection points
    %
    zgl = zeros(nl,1); % Array for most superior Z coordinates on lateral bony grid surface
    zgm = zeros(nm,1); % Array for most superior Z coordinates on medial bony grid surface

    for i=1:ns % Loop through subjects
        i_str=int2str(i); % Convert subject index to string
        % Load lateral cartilage triangular mesh
        eval(['trils = sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.trimesh;']);
        % Load medial cartilage triangular mesh
        eval(['trims = sd(' i_str ').' scan(j,:) '.tibia.cart.medial.trimesh;']);
        % Load lateral cartilage point coordinates
        eval(['xyzls = sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.points;']);
        % Load medial cartilage point coordinates
        eval(['xyzms = sd(' i_str ').' scan(j,:) '.tibia.cart.medial.points;']);
        %
        if j==1
            fstr=sd(i).FFE.tibia.bfnam; % FFE bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier

        elseif j==2
            fstr=sd(i).RHO.tibia.bfnam; % RHO bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier

        elseif j==3
            fstr=sd(i).T2S.tibia.bfnam; % T2S bone file name
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier
        end

        % Loop through Subjects (Tibias)
        %
        % Scale Grid
        %
        xgl = xql*scx(i); % Scale lateral X grid points
        ygl = yql*scy(i); % Scale lateral Y grid points
        xgm = xqm*scx(i); % Scale medial X grid points
        ygm = yqm*scy(i); % Scale medial Y grid points

        %
        % Get Bony Tibia Data
        %
        % Load lateral bone coordinates
        eval(['xyzl = sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points;']);
        % Load medial bone coordinates
        eval(['xyzm = sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points;']);
        % Load lateral bone mesh
        eval(['tril = sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.trimesh;']);
        % Load medial bone mesh
        eval(['trim = sd(' i_str ').' scan(j,:) '.tibia.bone.medial.trimesh;']);
        %
        % Project Grid onto Bony Surface Mesh
        %
        if j==1
            dist = 2.4; % Twice slice thickness for FFE images
        elseif j==2
            dist = 6; % Twice slice thickness for T1rho images
        elseif j==3
            dist = 4; % Twice slice thickness for T2S images
        end

        zgl(:,1) = gridproj(tril,xyzl,xgl,ygl,1,dist); % Project lateral grid onto bone surface
        zgm(:,1) = gridproj(trim,xyzm,xgm,ygm,1,dist); % Project medial grid onto bone surface
        %
        xyzgl = [xgl ygl zgl(:,1)]; % Combine lateral grid coordinates
        xyzgm = [xgm ygm zgm(:,1)]; % Combine medial grid coordinates
        %
        % Improve the Mesh
        %
        trigli = tri_fix2(trigl,xyzgl); % Fix lateral triangular mesh
        trigmi = tri_fix2(trigm,xyzgm); % Fix medial triangular mesh
        %%

        %
        % Plot Cartilage and Bone Surfaces
        %
        if iplt % Check if plotting is enabled

            % Plot Sagittal Cartilage and Gridded Bone Surfaces
            %
            hf2 = figure; % Create new figure for sagittal cartilage and bone surfaces
      
            % Plot lateral gridded bone surface
            hb2l = trimesh(trigli,xgl,ygl,zgl(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','k');
            hold on;
            % Plot medial gridded bone surface
            hb2m = trimesh(trigmi,xgm,ygm,zgm(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','k');
            %
            % Plot lateral cartilage surface
            hc2l = trimesh(trils,xyzls(:,1),xyzls(:,2), ...
                xyzls(:,3),'LineWidth',0.25, ...
                'FaceColor','none','EdgeColor','b');
            % Plot medial cartilage surface
            hc2m = trimesh(trims,xyzms(:,1),xyzms(:,2), ...
                xyzms(:,3),'LineWidth',0.25, ...
                'FaceColor','none','EdgeColor','b');
            xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
            ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
            zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
            % Add title to the plot
            title({fstr; 'Sagittal Cartilage and Bone Meshes'}, ...
                'Interpreter','none','FontSize',16,'FontWeight','bold');
            view(-60,12); % Set view angle
            axis equal; % Set equal scaling for axes

            psnam = fullfile(tdir,[fstr '_tcart08_sagittal_thk.pdf']); % Construct PDF file name
            orient landscape; % Set figure orientation to landscape
            set(hf2, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
            exportgraphics(hf2, psnam, "Resolution", 300); % Export figure to PDF
            close(hf2); % Close the figure

           
        end
        % %

        %
        % Calculate and Save Sagittal Cartilage Thicknesses
        %
        [ct,bi] = car_thk8(trils,xyzls,trigli,xyzgl); % Calculate lateral cartilage thicknesses
        cthkl(:,1) = ct; % Save lateral thicknesses
        xyzil(:,:,1) = bi; % Save lateral cartilage intersection points
        %
        [ct,bi] = car_thk8(trims,xyzms,trigmi,xyzgm); % Calculate medial cartilage thicknesses
        cthkm(:,1) = ct; % Save medial thicknesses
        xyzim(:,:,1) = bi; % Save medial cartilage intersection points
        %
        % Plot the Cartilage Thicknesses
        %
        if iplt % Check if plotting is enabled
            % %

            % Plot Sagittal Cartilage Thicknesses
            %

            hf4 = figure; % Create new figure for cartilage thicknesses
        
            % Plot lateral gridded bone surface
            hb4l = trimesh(trigli,xgl,ygl,zgl(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','b');
            hold on;
            % Plot medial gridded bone surface
            hb4m = trimesh(trigmi,xgm,ygm,zgm(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','b');
            %
            ctl = squeeze(cthkl(:,1)); % Extract lateral cartilage thicknesses
            idx = find(~isnan(ctl)); % Find non-NaN thickness indices
            it = nod2tri(idx,trigli,2); % Convert node indices to triangle indices
            % Plot lateral cartilage thickness surface
            hs4l = trisurf(trigli(it,:),xgl,ygl,zgl(:,1),ctl, ...
                'FaceColor','interp', ...
                'EdgeColor','b','LineWidth',0.25);
            ctm = squeeze(cthkm(:,1)); % Extract medial cartilage thicknesses
            idx = find(~isnan(ctm)); % Find non-NaN thickness indices
            it = nod2tri(idx,trigmi,2); % Convert node indices to triangle indices
            % Plot medial cartilage thickness surface
            hs4m = trisurf(trigmi(it,:),xgm,ygm,zgm(:,1),ctm, ...
                'FaceColor','interp', ...
                'EdgeColor','b','LineWidth',0.25);
            colormap jet; % Set color map to jet
            view(-60,12); % Set view angle
            axis equal; % Set equal scaling for axes
            %
            ct = [ctl; ctm]; % Combine lateral and medial thicknesses
            if min(ct)<0
                caxis([min(ct) max(ct)]); % Set color axis based on thickness range
            else
                caxis([0 max(ct)]); % Set color axis starting from 0
            end
            hcb4 = colorbar; % Add colorbar
            xlabel('X','FontSize',12,'FontWeight','bold'); % Label X-axis
            ylabel('Y','FontSize',12,'FontWeight','bold'); % Label Y-axis
            zlabel('Z','FontSize',12,'FontWeight','bold'); % Label Z-axis
            % Add title to the plot
            title({fstr; 'Sagittal Cartilage Thicknesses'}, ...
                'Interpreter','none','FontSize',16,'FontWeight','bold');
            orient landscape; % Set figure orientation to landscape
            set(hf4, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
            exportgraphics(hf4, psnam, "Resolution", 300, 'Append',true); % Append figure to PDF
            close(hf4); % Close the figure
        end % End if iplt
        %

        % %% Save all Figures to Combined PDF
        % h = findobj('type','figure'); % Commented out: find all open figures
        % n = length(h);
        % for i=1:n
        %     exportgraphics(figure(i), 'output.pdf', 'Append', true);
        % end

        %% Save Data and Close Plots

        %
        % Save Tibia Data
        %
        kid=fstr(1:5); % Extract knee identifier
        ileg=false; % Flag for leg (not used in this context)
        rnam = fullfile(tdir,[fstr '_tcart08_thk.mat']); % Construct results MAT file name
        % Save cartilage thickness data to MAT file
        save(rnam,'cthkl','cthkm','ileg','kid', ...
            'tribl','tribm','trils','trims', ...
            'xyzbl','xyzbm','xyzil','xyzim', ...
            'xyzls','xyzms','zgl','zgm');
        %
        %close all; % Commented out: close all figures
        %
        % if i ~= ns % Commented out: clear variables except for specified ones
        %     clearvars -except k ns bpnams bfnams cpnams cfnams scans conds ilegs kids...
        %         iplt iprt
        % end
    end
end
% Save subject data to a MAT file
main_file=fullfile(rdir, 'All Subjects Tibia Cartilage Thickness Data.mat');
save (main_file, 'sd');
return % End of script