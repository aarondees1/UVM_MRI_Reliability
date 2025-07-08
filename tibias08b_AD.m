%#######################################################################
%
%           * Tibia Bone Data Processing and Tibial Coordinate System Analysis *
%                         (Script Name: tibias08b_AD)
%
%          This program guides the user to select a data directory, then
%     reads and plots digitized points from a tibia in the tibia
%     coordinate system. It plots medial and lateral compartment
%     subchondral bone data for user-selected knees. The tibial
%     coordinate system for each knee is saved to a Matlab MAT file.
%     Additionally, a structure 'sd' containing selected subject data
%     is saved for subsequent processing by other scripts.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  See tibia_cs8.m for more information on the tibia
%             coordinate system.
%
%             3. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi6.m, rotxyz.m, tibia_cs8.m,
%             tri_area.m, tri_fix2.m and tri_norm.m must be in the
%             current directory or path.
%
%             4. This file outputs PDF files
%             with plots of the tibia coordinate system and the
%             subchondral bone into the directory of the knee data. 
%             ***_tibias08b.pdf,for a given subject, knee (left or right) 
%             and scan type (e.g. 001_L_FFE_tibias08b.pdf).
%
%             5. The subchondral bone data in the tibia coordinate
%             system and the coordinate transformation matrix and
%             origin are saved in the directory of the knee data in a
%             Matlab MAT file of similar format, ***_tibiaCS.mat.
%
%             6. This is an updated version of tibias07b.m.
%
%     7-July-2025 * Mack Gardner-Morse & Aaron Dees
%
%#######################################################################
%%
clear; % Clear all variables from the workspace
close all; % Close all open figures
clc; % Clear the command window
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2; % Minimum distance between distinct points to identify duplicates
iflag = true; % Flag to print a message if duplicate points are found

% Output Pdf and MAT File Names
%
psnam_suffix = '_tibias08b.pdf'; % Suffix for the output PDF file name
mnam_suffix = '_tibiaCS.mat'; % Suffix for the output MAT file name
%
% Compartment Colors and Labels
%
tclrs = [0 0 0.7; 0 0.5 0]; % Color definitions: deep blue for lateral, dark green for medial
%
tcmpt = ['lateral'
    'medial ']; % Labels for the tibial compartments

% Prompt user to select a directory containing subject data
div = uigetdir;
ddir = fullfile(div); % Construct full path to the selected directory
sd = dir(ddir); % Get all subdirectories in the selected directory
snams = {sd.name}'; % Extract names of subdirectories
id = startsWith(snams,{'.'; '..'}); % Identify current and parent directories
snams = snams(~id); % Remove current and parent directories from the list
sd=sd(~id);
sd=rmfield(sd,{'bytes','date','datenum','isdir'}); % Remove unnecessary fields from structure
% Display a dialog for the user to select subdirectories
idx = listdlg('ListString',snams,'ListSize',[200 600],'Name', ...
    'Tibias','PromptString', ...
    {'    Select all subdirectories for'; ...
    'calculating cartilage thicknesses.'});
%
% Check if user cancelled the selection
if isempty(idx)
    return; % Exit the script if no directories are selected
end
sd=sd(idx); % Keep only selected subdirectories
snams = string(snams(idx)); % Convert selected subject names to string array

%
ns = size(snams,1); % Number of selected subjects
nss = int2str(ns); % Convert number of subjects to string
% Initialize arrays to store paths for each subject
snam_paths=strings(ns,1);
ffe_paths=strings(ns,1);
rho_paths=strings(ns,1);
t2s_paths=strings(ns,1);

% Loop through each subject to collect file paths
for i=1:ns
    snam_paths(i)=fullfile(ddir,snams(i)); % Path to subject directory
    ffe_paths(i)=fullfile(snam_paths(i),'Visit 2','FFE'); % Path to FFE data
    rho_paths(i)=fullfile(snam_paths(i),'Visit 2','RHO'); % Path to RHO data
    t2s_paths(i)=fullfile(snam_paths(i),'Visit 2','T2S'); % Path to T2S data
    sd(i).folder=char(snam_paths(i)); % Store subject folder path

    % Get FFE Tibia Bone Names
    bnams=dir(fullfile(ffe_paths(i),'*_L_SAG_TIB*.csv')); % Find sagittal tibia CSV files
    bnams={bnams.name};
    % Filter out files with 'dup', 'MGG', or 'LD' in the name for any
    % duplicates, Mack tests, and loaded conditions
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).FFE.tibia.bfnam=char(bnams); % Store filtered FFE tibia file name

    % Get RHO Tibia Bone Names
    bnams=dir(fullfile(rho_paths(i),'Tibia','*_L_SAG_TIB*.csv')); % Find RHO sagittal tibia CSV files
    bnams={bnams.name};
    % Filter out files with 'dup', 'MGG', or 'LD' in the name for any
    % duplicates, Mack tests, and loaded conditions
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).RHO.tibia.bfnam=char(bnams); % Store filtered RHO tibia file name

    % Get T2S Tibia Bone Names
    bnams=dir(fullfile(t2s_paths(i),'Tibia','*_L_SAG_TIB*.csv')); % Find T2S sagittal tibia CSV files
    bnams={bnams.name};
    % Filter out files with 'dup', 'MGG', or 'LD' in the name for any
    % duplicates, Mack tests, and loaded conditions
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).T2S.tibia.bfnam=char(bnams); % Store filtered T2S tibia file name

    % Get AX Tibia Files for FFE
    ax_list=dir(fullfile(ffe_paths(i),'*_L_AX_TIB*.csv')); % Find axial tibia CSV files
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out files with 'LD'
    sd(i).FFE.tibia.axfnam=char(ax_list); % Store filtered FFE axial file name

    % Get AX Tibia Files for RHO
    ax_list=dir(fullfile(rho_paths(i),'Tibia','*_L_AX_TIB*.csv')); % Find RHO axial tibia CSV files
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out files with 'LD'
    sd(i).RHO.tibia.axfnam=char(ax_list); % Store filtered RHO axial file name

    % Get AX Tibia Files for T2S
    ax_list=dir(fullfile(t2s_paths(i),'Tibia','*_L_AX_TIB*.csv')); % Find T2S axial tibia CSV files
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out files with 'LD'
    sd(i).T2S.tibia.axfnam=char(ax_list); % Store filtered T2S axial file name
end
ffe_paths=ffe_paths'; % Transpose FFE paths
rho_paths=rho_paths'; % Transpose RHO paths
t2s_paths=t2s_paths'; % Transpose T2S paths

% Check if no subjects were selected
if isequal(snams,0)
    return; % Exit the script if no subjects are selected
end
%%

% Setup Tibial Coordinate System Figure
%
% Make Results Folder
%
rdir = fullfile(ddir,'Results','Bone'); % Define path for results directory
if ~isfolder(rdir)
    mkdir(rdir); % Create results directory if it doesn't exist
end
%
% Loop through each subject
for i=1:ns
    % Loop through FFE, RHO, and T2S data
    for j=1:3
        % Select appropriate file paths based on data type
        if j==1
            ax_path=fullfile(ffe_paths(i),sd(i).FFE.tibia.axfnam); % FFE axial file path
            fstr=sd(i).FFE.tibia.bfnam; % FFE bone file name
            bone_path=fullfile(ffe_paths(i),sd(i).FFE.tibia.bfnam); % FFE bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name

        elseif j==2
            ax_path=fullfile(rho_paths(i),'Tibia',sd(i).RHO.tibia.axfnam); % RHO axial file path
            fstr=sd(i).RHO.tibia.bfnam; % RHO bone file name
            bone_path=fullfile(rho_paths(i),'Tibia',fstr); % RHO bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name

        elseif j==3
            ax_path=fullfile(t2s_paths(i),'Tibia',sd(i).T2S.tibia.axfnam); % T2S axial file path
            fstr=sd(i).T2S.tibia.bfnam; % T2S bone file name
            bone_path=fullfile(t2s_paths(i),'Tibia',fstr); % T2S bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract identifier from file name
        end

        % Create a new figure for the tibial coordinate system
        hf1 = figure;
        orient tall; % Set figure orientation to tall
        view(3); % Set 3D view
        hold on; % Allow multiple plots on the same figure
        %
        % Get Tibial Coordinate System, Proximal Tibia Outline (PTO)
        %
        % NOTE: PTO is in the tibia coordinate system.
        %

        % Call function to compute tibia coordinate system and related data
        [xyzc,xyzr,aspect,widt,height,xyzpto,prox_tib_ml,prox_tib_ap] = tibia_cs8(ax_path, ...
            false,true);
        %
        % Finish Plot
        %
        % Add title to the plot with file identifier and description
        title({fstr; 'Tibia Coordinate System'},'FontSize',16, ...
            'FontWeight','bold', 'Interpreter','none');
        axis equal; % Set equal scaling for axes
        %
        psnam = [fstr psnam_suffix]; % Construct PDF file name
        psnam = fullfile(rdir,psnam); % Construct full path for PDF
        %
        % Maximize figure window
        set(hf1, 'units','normalized','outerposition',[0 0 1 1]);
        % Export figure to PDF
        exportgraphics(hf1, psnam, "Resolution", 300);
        close(hf1); % Close the figure
       
        % Get Raw Tibia Bone Data
        %
        roi = rd_roi6(bone_path); % Read region of interest data from CSV
        %
        roinams = upper(char(roi.name)); % Convert ROI names to uppercase
        idc(2) = find(strcmp('MTCS',cellstr(roinams))); % Find index for medial compartment
        idc(1) = find(strcmp('LTCS',cellstr(roinams))); % Find index for lateral compartment
        %
        % Raw Data Figure
        %
        % Create a new figure for raw data
        hf2 = figure;
        orient landscape; % Set figure orientation to landscape
        view(3); % Set 3D view
        hold on; % Allow multiple plots
        xlabel('X (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis
        ylabel('Y (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
        zlabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
        % Add title to raw data plot
        title({[fstr ' - MRI CS']; 'Blue - Lateral, Green - Medial'}, ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        % Transformed Data Figure
        %
        % Create a new figure for transformed data

        hf3 = figure;
        orient landscape; % Set figure orientation to landscape
        view(3); % Set 3D view
        hold on; % Allow multiple plots
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis (anterior-posterior)
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
        % Add title to transformed data plot
        title([fstr ' - Tibial CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %%
        %
        % Loop through Compartments
        %

        % Process lateral and medial compartments
        for l = 1:2

            %
            % Get Data
            %
            dat1 = roi(idc(l)).data'; % Extract data for the current compartment
            nsl1 = size(dat1,1); % Number of slices in the data
            %
            % Check for Duplicates, Direction of Digitization and Transform Data to
            % Tibia Coordinate System
            %
            datt = cell(nsl1,1); % Initialize cell array for transformed data
            %
            % Process each slice
            for n = 1:nsl1
                xyz = dat1{n}; % Get coordinates for the current slice
                % Check if slice data is empty
                if isempty(xyz)
                    cmprt = deblank(tcmpt(l,:)); % Get compartment name
                    % Throw error if no coordinates are found
                    error([' *** ERROR in tibias08b:  No coordinates for ', ...
                        'slice ' int2str(n) ' in ' cmprt ' compartment!']);
                end
                % Remove duplicate points
                xyz = fix_pts_AD(xyz,tol,iflag);
                npts = size(xyz,1); % Number of points in the slice
                [~,imx] = max(xyz(:,2)); % Find index of maximum Y-coordinate
                [~,imn] = min(xyz(:,2)); % Find index of minimum Y-coordinate
                % Ensure anterior-to-posterior ordering
                if imn<imx
                    xyz = xyz(1:npts,:); % Keep order if already anterior to posterior
                else
                    xyz = xyz(npts:-1:1,:); % Reverse order to anterior to posterior
                end
                xyzt = xyz-repmat(xyzc,npts,1); % Center data by subtracting origin
                xyzt = xyzt*xyzr; % Transform data to tibia coordinate system
                datt{n} = xyzt; % Store transformed data
                %
                %
                % Plot Raw Data
                %
                figure(hf2); % Switch to raw data figure
                %
                % Plot raw data points with compartment-specific color
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','Color',tclrs(l,:), ...
                    'MarkerSize',8,'LineWidth',1);
                % Add slice number label
                if xyz(1,2)<0
                    text(xyz(1,1),xyz(1,2)-0.5,xyz(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                else
                    text(xyz(1,1),xyz(1,2)+0.5,xyz(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                end
                %
                % Plot Transformed Data
                %
                figure(hf3); % Switch to transformed data figure
                %
                % Plot transformed data points with compartment-specific color
                plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color',tclrs(l,:), ...
                    'MarkerSize',8,'LineWidth',1);
                % Add slice number label
                if xyzt(1,1)<0
                    text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                else
                    text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                end
            end

            %
            % Save Data into Compartment Specific Variables
            %

            %
            % Get Surface Triangulation
            %
            trit = mk_tri6(datt); % Create triangulation for the compartment
            xyzt = cell2mat(datt); % Convert cell array to matrix
            trit = tri_fix2(trit, xyzt); % Fix triangulation issues

            % Store data in compartment-specific variables
            eval(['dat' tcmpt(l,1) ' = datt;']);
            eval(['tri' tcmpt(l,1) ' = trit;']);
            eval(['xyz' tcmpt(l,1) ' = xyzt;']);

            %
            %   clear datt trit xyzt; % Clear temporary variables (commented out)
            %
        end

        %%
        %
        % Save Plots
        %
        % Maximize figure windows
        set(hf2, 'units','normalized','outerposition',[0 0 1 1]);
        set(hf3, 'units','normalized','outerposition',[0 0 1 1]);
        % Append raw and transformed data plots to PDF
        exportgraphics(hf2, psnam, "Resolution", 300, 'Append', true);
        exportgraphics(hf3, psnam, "Resolution", 300, 'Append', true);
        close(hf2); % Close raw data figure
        close(hf3); % Close transformed data figure
        %
        % Plot Transformed Bone Surface Data
        %
        % Create a new figure for surface data
        hf4 = figure;
        orient landscape; % Set figure orientation to landscape
        %
        % Plot lateral compartment surface
        trisurf(tril,xyzl(:,1),xyzl(:,2),xyzl(:,3),'FaceColor','interp', ...
            'EdgeColor',tclrs(1,:),'LineWidth',1);
        hold on;
        % Plot medial compartment surface
        trisurf(trim,xyzm(:,1),xyzm(:,2),xyzm(:,3),'FaceColor','interp', ...
            'EdgeColor',tclrs(2,:),'LineWidth',1);
        axis equal; % Set equal scaling for axes
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
        % Add title to surface plot
        title([fstr ' - Tibial CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        % Maximize figure window
        set(hf4, 'units','normalized','outerposition',[0 0 1 1]);
        % Append surface plot to PDF
        exportgraphics(hf4, psnam, "Resolution", 300, 'Append', true);
        close(hf4); % Close surface figure
        %
        % Save Data into a Matlab MAT File for Further Processing
        %
        mnam = [fstr mnam_suffix]; % Construct MAT file name
        mnam = fullfile(rdir,mnam); % Construct full path for MAT file
        kid=fstr(1:5); % Extract knee identifier
        ileg=false; % Flag for leg (not used in this context)
        %
        % Save data to MAT file
        save(mnam, 'datl','datm','tril','trim','xyzc','xyzr', ...
            'xyzl','xyzm','xyzpto','prox_tib_ml','prox_tib_ap');
        %
    end
end
% Save subject data to a MAT file
main_file=fullfile(ddir,'Results', 'All Subjects Tibia Bone Data.mat');
save (main_file, 'sd');
return % End of script