%#######################################################################
%
%           * Femoral Bone Data Processing and Femoral Coordinate System Analysis *
%                         (Script Name: femurs08b_AD)
%
%          This program guides the user to select a data directory,
%     then reads and plots digitized points from a femur. It computes
%     and transforms the raw data into the femoral coordinate system.
%     Specifically, it processes and plots subchondral bone data for
%     the lateral and medial condyles, and the trochlear region of
%     user-selected knees.
%     This script utilizes the 'sd' structure from 'All Subjects Tibia Cartilage Thickness Data.mat',
%     which is an output of the 'tcart08_AD' script, indicating its position
%     in the overall data processing pipeline.
%
%     The femoral coordinate system and the transformed subchondral bone
%     data are saved to a Matlab MAT file. Various plots, including the
%     femoral coordinate system, raw MRI data, transformed data, and the
%     transformed bone surface, are generated and saved as PDF files.
%     The 'sd' structure is updated with new femoral bone data
%     and saved for subsequent processing by other scripts.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  See f_cs_14.m for more information on the femoral
%             coordinate system.
%
%             3.  The angles between slice ends in the triangulation is
%             controlled by the angle input (in degrees) into
%             mk_tri4f.m. Typical angles are 5-15 degrees.
%
%             4. The Matlab M-files comb_dat.m, cyl_fit.m, cyl_plt.m,
%             dist2cyl.m, f_cs_14.m, fix_pts.m, mk_tri4f.m, plane_fit.m,
%             plt_datsl.m, pt2line.m, pts2lin.m, rd_roi6.m, sl_dir.m,
%             tri_area.m, tri_fix2.m, tri_norm.m and xzpl2pol.m must be
%             in the current directory or path.
%
%             5.  This file outputs PDF file*, e.g., '***_femur08b.pdf',
%             containing plots of the femoral coordinate system, raw and
%             transformed subchondral bone data, and the bone surface. 
%             The format follows ***_femur08b.pdf for a given subject, knee 
%             (left or right) and scan type (e.g. 001_L_FFE_femur08b.pdf).
%
%             6.  The subchondral bone data in the femoral coordinate
%             system, along with the coordinate transformation matrix and
%             origin, are saved in a Matlab MAT file named
%             '***_femurCS.mat'.
%
%     7-July-2025 * Mack Gardner-Morse & Aaron Dees
%
%#######################################################################
clear; % Clear all variables from workspace
close all; % Close all open figures
clc; % Clear command window
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2; % Set minimum distance threshold for distinct points
iflag = true; % Enable flag to print duplicate point messages
%
% Output PS and MAT File Names
%
psnam_suffix = '_femur08b.pdf'; % Define suffix for PDF output file
mnam_suffix = '_femurCS.mat'; % Define suffix for MAT output file
%
% Regional Colors and Labels
%
fclrs = [0 0 0.7; 0 0.5 0; 0.7 0 0]; % Define colors: deep blue, dark green, red
%
fregs = ['lateral '
         'medial  '
         'trochlea']; % Define labels for femoral regions
%
% Get Sagittal Bone CSV File Name
%
div = uigetdir; % Open dialog to select base directory
ddir = fullfile(div); % Construct full path to selected directory
rdir = fullfile(div,'Results'); % Construct path to results directory
bdir = fullfile(rdir,'Bone'); % Construct path to bone data directory
tdir = fullfile(rdir,'Thickness'); % Construct path to thickness data directory
gdir = fullfile(rdir,'Grids'); % Construct path to grids directory
load(fullfile(rdir,'All Subjects Tibia Cartilage Thickness Data.mat')); % Load subject data from MAT file
ns=size(sd,1); % Get number of subjects
nss = int2str(ns); % Convert number of subjects to string
ffe_paths=strings(ns,1); % Initialize array for FFE scan paths
rho_paths=strings(ns,1); % Initialize array for RHO scan paths
t2s_paths=strings(ns,1); % Initialize array for T2S scan paths
%
%
% Check for Axial Plane File
%
for i=1:ns % Loop through each subject
    ffe_paths(i)=fullfile(sd(i).folder,'Visit 2','FFE'); % Set path for FFE scan directory
    rho_paths(i)=fullfile(sd(i).folder,'Visit 2','RHO'); % Set path for RHO scan directory
    t2s_paths(i)=fullfile(sd(i).folder,'Visit 2','T2S'); % Set path for T2S scan directory
%
    %Get FFE Femur Bone Names
    bnams=dir(fullfile(ffe_paths(i),'*_L_SAG_FEM*.csv')); % Find FFE sagittal femur CSV files
    bnams={bnams.name}; % Extract file names
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)&~contains(bnams,'MGG')&~contains(bnams,'LD')); % Filter out unwanted files
    sd(i).FFE.femur.bfnam=char(bnams); % Store FFE femur file name in subject data
%
    %Get RHO Femur Bone Names
    bnams=dir(fullfile(rho_paths(i),'Femur','*_L_SAG_FEM*.csv')); % Find RHO sagittal femur CSV files
    bnams={bnams.name}; % Extract file names
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)&~contains(bnams,'MGG')&~contains(bnams,'LD')); % Filter out unwanted files
    sd(i).RHO.femur.bfnam=char(bnams); % Store RHO femur file name in subject data
%
    %Get T2S Femur Bone Names
    bnams=dir(fullfile(t2s_paths(i),'Femur','*_L_SAG_FEM*.csv')); % Find T2S sagittal femur CSV files
    bnams={bnams.name}; % Extract file names
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)&~contains(bnams,'MGG')&~contains(bnams,'LD')); % Filter out unwanted files
    sd(i).T2S.femur.bfnam=char(bnams); % Store T2S femur file name in subject data
%
    %Get AX Femur Files
    ax_list=dir(fullfile(ffe_paths(i),'*_L_AX_FEM*.csv')); % Find FFE axial femur CSV files
    ax_list={ax_list.name}; % Extract file names
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out unwanted files
    sd(i).FFE.femur.axfnam=char(ax_list); % Store FFE axial femur file name in subject data
%
    ax_list=dir(fullfile(rho_paths(i),'Femur','*_L_AX_FEM*.csv')); % Find RHO axial femur CSV files
    ax_list={ax_list.name}; % Extract file names
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out unwanted files
    sd(i).RHO.femur.axfnam=char(ax_list); % Store RHO axial femur file name in subject data
%
    ax_list=dir(fullfile(t2s_paths(i),'Femur','*_L_AX_FEM*.csv')); % Find T2S axial femur CSV files
    ax_list={ax_list.name}; % Extract file names
    ax_list=ax_list(~contains(ax_list,'LD')); % Filter out unwanted files
    sd(i).T2S.femur.axfnam=char(ax_list); % Store T2S axial femur file name in subject data
%
end
ffe_paths=ffe_paths'; % Transpose FFE paths array
rho_paths=rho_paths'; % Transpose RHO paths array
t2s_paths=t2s_paths'; % Transpose T2S paths array
%
% Setup Femoral Coordinate System Figure
%
if ~isfolder(bdir) % Check if bone directory exists
    mkdir(rdir); % Create results directory if it doesn't exist
end
for i=1:ns % Loop through each subject
    for j=1:3 % Loop through scan types (FFE, RHO, T2S)
        if j==1 % Check for FFE scan type
            ax_path=fullfile(ffe_paths(i),sd(i).FFE.femur.axfnam); % Construct FFE axial file path
            fstr=sd(i).FFE.femur.bfnam; % Get FFE femur file name
            bone_path=fullfile(ffe_paths(i),sd(i).FFE.femur.bfnam); % Construct FFE bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
%
        elseif j==2 % Check for RHO scan type
            ax_path=fullfile(rho_paths(i),'Femur',sd(i).RHO.femur.axfnam); % Construct RHO axial file path
            fstr=sd(i).RHO.femur.bfnam; % Get RHO femur file name
            bone_path=fullfile(rho_paths(i),'Femur',fstr); % Construct RHO bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
%
        elseif j==3 % Check for T2S scan type
            ax_path=fullfile(t2s_paths(i),'Femur',sd(i).T2S.femur.axfnam); % Construct T2S axial file path
            fstr=sd(i).T2S.femur.bfnam; % Get T2S femur file name
            bone_path=fullfile(t2s_paths(i),'Femur',fstr); % Construct T2S bone file path
            fstr=[fstr(1:6) fstr(15:17)]; % Extract subject ID and scan suffix
%
        end
hf1 = figure; % Create new figure for femoral coordinate system
orient landscape; % Set figure orientation to landscape
view(3); % Set 3D view
hold on; % Enable hold for multiple plots
%
% Get MRI Sagittal Femur Bone Data
%
roi = rd_roi6(bone_path); % Read sagittal femur ROI data from CSV
%
roinams = {roi.name}'; % Extract ROI names into cell array
%roinams = upper(char(roi.name));
ids(3) = find(startsWith(roinams,'TRO')); % Find index of trochlea ROI
ids(2) = find(startsWith(roinams,'MED')); % Find index of medial condyle ROI
ids(1) = find(startsWith(roinams,'LAT')); % Find index of lateral condyle ROI
%
datlb = roi(ids(1)).data'; % Extract lateral condyle data
datmb = roi(ids(2)).data'; % Extract medial condyle data
dattb = roi(ids(3)).data'; % Extract trochlea data
%
% Get MRI Axial Femoral Shaft Outline (FS) and Center Point (CP) Data
%
roia = rd_roi6(ax_path); % Read axial femur ROI data from CSV
%
roinama = {roia.name}'; % Extract axial ROI names into cell array
ida(2) = find(startsWith(roinama,'PROX')); % Find index of proximal femoral shaft ROI
ida(1) = find(startsWith(roinama,'CENT')); % Find index of center point ROI
%
datcp = roia(ida(1)).data; % Extract center point data
datfs = roia(ida(2)).data'; % Extract proximal femoral shaft data
%
% Get Femoral Coordinate System
%
[xyzc,xyzr,xyzl,r] = f_cs_14(datlb,datmb,datfs,datcp,true,hf1); % Compute femoral coordinate system
%
% Finish Plot
%
title({fstr; 'Femoral Coordinate System'},'FontSize',16,'FontWeight','bold', 'Interpreter','none'); % Add title with subject ID and description
axis equal; % Set equal scaling for axes
%
psnam = [fstr psnam_suffix]; % Construct PDF file name with suffix
psnam = fullfile(bdir,psnam); % Construct full path to PDF file
%
set(hf1, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window to full screen
exportgraphics(hf1, psnam, "Resolution", 300); % Export figure to PDF with 300 DPI
close(hf1); % Close femoral coordinate system figure
%
% Raw Data Figure
%
hf2 = figure; % Create new figure for raw MRI data
orient landscape; % Set figure orientation to landscape
view(3); % Set 3D view
hold on; % Enable hold for multiple plots
xlabel('X (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis
ylabel('Y (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
zlabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
title({[fstr ' - MRI CS']; ['Blue - Lateral, Green - Medial, ','Red - Trochlea']},'FontSize',16,'FontWeight','bold','Interpreter','none'); % Add title with color legend
%
% Transformed Data Figure
%
hf3 = figure; % Create new figure for transformed data
orient landscape; % Set figure orientation to landscape
view(3); % Set 3D view
hold on; % Enable hold for multiple plots
xlabel('AP (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis (anterior-posterior)
ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
zlabel('Superior (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
title({[fstr ' - Femoral CS']; ['Blue - Lateral, Green - Medial, ','Red - Trochlea']},'FontSize',16,'FontWeight','bold','Interpreter','none'); % Add title with color legend
%
% Loop through Regions
%
for l = 1:3 % Loop through lateral, medial, and trochlea regions
%
% Get Data
%
   eval(['dat1 = dat' fregs(l,1) 'b;']); % Dynamically extract region data
   nsl1 = size(dat1,1); % Get number of slices in region
%
% Check for Duplicates, Direction of Digitization and Transform Data to
% Femoral Coordinate System
%
   datt = cell(nsl1,1); % Initialize cell array for transformed data
%
   for n = 1:nsl1 % Loop through each slice
      xyz = dat1{n}; % Extract coordinates for current slice
      if isempty(xyz) % Check if slice data is empty
        reg = deblank(fregs(l,:)); % Get region name without trailing spaces
        error([' *** ERROR in femurs08b:  No coordinates for ','slice ' int2str(n) ' in ' reg ' region!']); % Throw error for empty slice
      end
      xyz = fix_pts_AD(xyz,tol,iflag,fstr); % Remove duplicate points using tolerance
      npts = size(xyz,1); % Get number of points in slice
      [~,imx] = max(xyz(:,2)); % Find index of maximum Y-coordinate
      [~,imn] = min(xyz(:,2)); % Find index of minimum Y-coordinate
      if imn<imx % Check digitization direction
        xyz = xyz(1:npts,:); % Keep anterior-to-posterior order
      else
         xyz = xyz(npts:-1:1,:); % Reverse to anterior-to-posterior order
      end
      xyzt = xyz-repmat(xyzc,npts,1); % Center data relative to coordinate origin
      xyzt = xyzt*xyzr; % Transform data to femoral coordinate system
      datt{n} = xyzt; % Store transformed data for slice
% 
% Plot Raw Data
%
      figure(hf2); % Select raw data figure
%
      plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','Color',fclrs(l,:),'MarkerSize',8,'LineWidth',1); % Plot raw slice data in region color
      if xyz(1,2)<0 % Check Y-position for label placement
        text(xyz(1,1),xyz(1,2)-0.5,xyz(1,3),int2str(n),'Color','k','FontSize',10); % Add slice number below point
      else
        text(xyz(1,1),xyz(1,2)+0.5,xyz(1,3),int2str(n),'Color','k','FontSize',10); % Add slice number above point
      end
% 
% Plot Transformed Data
%
      figure(hf3); % Select transformed data figure
%
      plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color',fclrs(l,:),'MarkerSize',8,'LineWidth',1); % Plot transformed slice data in region color
      if xyzt(1,1)<0 % Check X-position for label placement
        text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n),'Color','k','FontSize',10); % Add slice number left of point
      else
        text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n),'Color','k','FontSize',10); % Add slice number right of point
      end
   end
%
% Save Data into Region Specific Variables
%
   eval(['dat' fregs(l,1) 'bt = datt;']); % Store transformed data in region-specific variable
%
   clear datt; % Clear temporary transformed data
%
end
%
% Save Plots
%
set(hf2, 'units','normalized','outerposition',[0 0 1 1]); % Maximize raw data figure window
set(hf3, 'units','normalized','outerposition',[0 0 1 1]); % Maximize transformed data figure window
exportgraphics(hf2, psnam,"Resolution", 300, 'Append', true); % Append raw data plot to PDF
exportgraphics(hf3, psnam,"Resolution", 300, 'Append', true); % Append transformed data plot to PDF
close(hf2); % Close raw data figure
close(hf3); % Close transformed data figure
%
% Combine Regions
%
datfb = comb_dat(datlbt,datmbt,dattbt); % Combine transformed data from all regions
%
% Get Surface Triangulation
%
trifb = mk_tri4f(datfb,15); % Create triangular mesh with 15-degree angle smoothing
xyzfb = cell2mat(datfb); % Convert combined data to matrix
trifb = tri_fix2(trifb,xyzfb); % Fix triangulation errors
%
% Plot Transformed Bone Surface Data
% 
hf4 = figure; % Create new figure for bone surface
orient landscape; % Set figure orientation to landscape
%
trisurf(trifb,xyzfb(:,1),xyzfb(:,2),xyzfb(:,3),'FaceColor','interp','EdgeColor',fclrs(1,:),'LineWidth',1); % Plot triangulated bone surface
axis equal; % Set equal scaling for axes
xlabel('AP (mm)','FontSize',12,'FontWeight','bold'); % Label X-axis (anterior-posterior)
ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold'); % Label Y-axis
zlabel('Superior (mm)','FontSize',12,'FontWeight','bold'); % Label Z-axis
title([fstr ' - Femoral CS'],'FontSize',16,'FontWeight','bold','Interpreter','none'); % Add title with subject ID
%
set(hf4, 'units','normalized','outerposition',[0 0 1 1]); % Maximize bone surface figure window
exportgraphics(hf4, psnam,"Resolution", 300, 'Append', true); % Append bone surface plot to PDF
close(hf4); % Close bone surface figure
%
%Save Data into a Matlab MAT File for Further Processing
%
        mnam = [fstr mnam_suffix]; % Construct MAT file name with suffix
        mnam = fullfile(bdir,mnam); % Construct full path to MAT file
        kid=fstr(1:5); % Extract knee identifier
        ileg=false; % Set leg indicator to false (left knee)
%
save(mnam,'datlb','datlbt','datmb','datmbt','dattb','dattbt','datfb','kid','ileg','r','trifb','xyzc','xyzfb','xyzl','xyzr'); % Save data to MAT file
%
    end
end
main_file=fullfile(ddir,'Results', 'All Subjects Full Tibia and Femur Bone Data.mat'); % Construct path to main MAT file
save (main_file, 'sd'); % Save subject data to main MAT file
return % Exit the script