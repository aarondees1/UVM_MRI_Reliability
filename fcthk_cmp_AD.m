%#######################################################################
%
%           * Femoral Cartilage Thickness Comparison Program *
%                      (Script Name: fcthk_cmp_AD)
%
%          This M-File processes and compares femoral cartilage
%     thicknesses obtained from T1FFE, T1rho, and T2S MRI scans. It reads
%     the previously calculated thickness data, establishes a combined
%     1 mm by 2 degree grid that covers all data sets, and then calculates
%     the differences in cartilage thicknesses between the scan types.
%
%          The script generates plots showing the individual cartilage
%     thickness maps for each scan type, as well as the thickness
%     differences (T1rho - T1FFE and T2S - T1FFE). It also produces
%     statistical plots of these differences (means, standard deviations,
%     and subject numbers contributing to each grid point). All generated
%     plots are saved as a single PDF file, and regional average
%     thicknesses are exported to Excel spreadsheets for detailed analysis.
%
%     NOTES:  1.  Both grids must have integer coordinates.
%
%             2.  Bone and cartilage MAT files must already exist in
%             the selected subject directories. These are outputs from
%             'femurs08b_AD' and 'femurs08c_AD' respectively.
%
%             3.  This program requires the grid scaling MAT files in the
%             format ***_fgrid08_.mat for a given scan type (e.g.,
%             FFE_fgrid08_.mat). It also requires the previously calculated
%             thickness MAT files in the format ***_fcart08_thk.mat for a
%             given subject, knee (left or right) and scan type (e.g.,
%             001_L_FFE_fcart08_thk.mat).
%
%             4.  This M-file outputs PDF files containing plots of
%             the femoral cartilage thicknesses and their differences
%             (Femoral Thickness Differences.pdf).
%             ROI data for each scan type is also exported to
%             Excel files in the format fem_lat_pos.xlsx,
%             fem_med_pos.xlsx, etc.
%
%     7-July-2025 * Mack Gardner-Morse & Aaron Dees
%
%#######################################################################
%%
clear; % Clear all variables from workspace
close all; % Close all open figures
clc; % Clear command window
grayColor = [.7 .7 .7]; % Define gray color for potential use
div = uigetdir; % Open dialog to select base directory
rdir = fullfile(div,'Results'); % Construct path to results directory
bdir = fullfile(rdir,'Bone'); % Construct path to bone data directory
tdir = fullfile(rdir,'Thickness'); % Construct path to thickness data directory
gdir = fullfile(rdir,'Grids'); % Construct path to grids directory
load(fullfile(rdir, 'All Subjects Tibia and Femur Cartilage Thickness Data.mat')); % Load subject data from MAT file
ns=size(sd,1); % Get number of subjects
% Get Analysis Grids
%
ffeg = load(fullfile(gdir,'FFE_fgrid08_.mat')); % Load FFE grid data
rhog = load(fullfile(gdir,'RHO_fgrid08_.mat')); % Load RHO grid data
t2sg = load(fullfile(gdir,'T2S_fgrid08_.mat')); % Load T2S grid data
% Get Coordinates and Ranges for Lateral Compartment
%
tq_ffe = ffeg.tq; % Extract theta coordinates for FFE grid
zq_ffe = ffeg.zq; % Extract Z coordinates for FFE grid
%
tq_rho = rhog.tq; % Extract theta coordinates for RHO grid
zq_rho = rhog.zq; % Extract Z coordinates for RHO grid
%
tq_t2s = t2sg.tq; % Extract theta coordinates for T2S grid
zq_t2s = t2sg.zq; % Extract Z coordinates for T2S grid
%
tmin_ffe = min(tq_ffe); % Get minimum theta for FFE
tmax_ffe = max(tq_ffe); % Get maximum theta for FFE
nc_ffe = ffeg.nt; % Get number of theta columns for FFE
%
tmin_rho = min(tq_rho); % Get minimum theta for RHO
tmax_rho = max(tq_rho); % Get maximum theta for RHO
nc_rho = rhog.nt; % Get number of theta columns for RHO
%
tmin_t2s = min(tq_t2s); % Get minimum theta for T2S
tmax_t2s = max(tq_t2s); % Get maximum theta for T2S
nc_t2s = t2sg.nt; % Get number of theta columns for T2S
%
zmin_ffe = min(zq_ffe); % Get minimum Z for FFE
zmax_ffe = max(zq_ffe); % Get maximum Z for FFE
nr_ffe = ffeg.nz; % Get number of Z rows for FFE
zmin_rho = min(zq_rho); % Get minimum Z for RHO
zmax_rho = max(zq_rho); % Get maximum Z for RHO
nr_rho = rhog.nz; % Get number of Z rows for RHO
zmin_t2s = min(zq_t2s); % Get minimum Z for T2S
zmax_t2s = max(zq_t2s); % Get maximum Z for T2S
nr_t2s = t2sg.nz; % Get number of Z rows for T2S
%%
%Get Combined Grid for Both Data Sets
tmin=[tmin_ffe; tmin_rho; tmin_t2s]; % Combine minimum theta values
tmin=min(tmin); % Get overall minimum theta
tmax = [tmax_ffe; tmax_rho; tmax_t2s]; % Combine maximum theta values
tmax=max(tmax); % Get overall maximum theta
zmin = [zmin_ffe; zmin_rho; zmin_t2s]; % Combine minimum Z values
zmin = min(zmin); % Get overall minimum Z
zmax = [zmax_ffe; zmax_rho; zmax_t2s]; % Combine maximum Z values
zmax = max(zmax); % Get overall maximum Z
[tg,zg] = meshgrid(tmin:2:tmax,zmin:zmax); % Create combined grid
[nr,nc] = size(tg); % Get dimensions of combined grid
%%
%
% Get Indexes into Combined Grid
%
n = nr*nc; % Calculate total number of grid points
idx = reshape(1:n,nr,nc); % Create index matrix for grid points
%
offstcf = (tmin_ffe-tmin)/2+1; % Calculate column offset for FFE
offstrf = (zmin_ffe-zmin+1); % Calculate row offset for FFE
idxf = idx(offstrf:offstrf+nr_ffe-1,offstcf:offstcf+nc_ffe-1); % Get FFE indices in combined grid
%
offstcr = (tmin_rho-tmin)/2+1; % Calculate column offset for RHO
offstrr = (zmin_rho-zmin+1); % Calculate row offset for RHO
idxr = idx(offstrr:offstrr+nr_rho-1,offstcr:offstcr+nc_rho-1); % Get RHO indices in combined grid
%
offstct = (tmin_t2s-tmin)/2+1; % Calculate column offset for T2S
offstrt = (zmin_t2s-zmin+1); % Calculate row offset for T2S
idxt = idx(offstrt:offstrt+nr_t2s-1,offstct:offstct+nc_t2s-1); % Get T2S indices in combined grid
%%
%
% Read Cartilage Thicknesses
%
cthkd_rf = zeros(nr,nc,ns); % Initialize array for RHO-FFE thickness differences
cthkd_tf = zeros(nr,nc,ns); % Initialize array for T2S-FFE thickness differences
%%
cthkf = []; % Initialize FFE thickness array
cthkr = []; % Initialize RHO thickness array
cthkt = []; % Initialize T2S thickness array
%%
for i=1:ns % Loop through each subject
    fstr=sd(i).FFE.femur.bfnam; % Get FFE bone file name
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat']; % Construct FFE thickness file name
    ffe = load(fullfile(tdir,fstr)); % Load FFE thickness data
    fstr=sd(i).RHO.femur.bfnam; % Get RHO bone file name
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat']; % Construct RHO thickness file name
    rho = load(fullfile(tdir,fstr)); % Load RHO thickness data
    fstr=sd(i).T2S.femur.bfnam; % Get T2S bone file name
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat']; % Construct T2S thickness file name
    t2s = load(fullfile(tdir,fstr)); % Load T2S thickness data
    cthkf1 = ffe.cthk; % Extract FFE thickness values
    cthkr1 = rho.cthk; % Extract RHO thickness values
    cthkt1 = t2s.cthk; % Extract T2S thickness values
    %%
    %
    % Get Cartilage Thicknesses at Similar Coordinates
    %
    cthkf2 = NaN(nr,nc); % Initialize FFE thickness grid with NaNs
    cthkr2 = NaN(nr,nc); % Initialize RHO thickness grid with NaNs
    cthkt2 = NaN(nr,nc); % Initialize T2S thickness grid with NaNs
    %
    cthkf2(idxf) = cthkf1; % Map FFE thicknesses to combined grid
    cthkr2(idxr) = cthkr1; % Map RHO thicknesses to combined grid
    cthkt2(idxt) = cthkt1; % Map T2S thicknesses to combined grid
    cthkf = cat(3,cthkf,cthkf2); % Append FFE thickness grid to 3D array
    cthkr = cat(3,cthkr,cthkr2); % Append RHO thickness grid to 3D array
    cthkt = cat(3,cthkt,cthkt2); % Append T2S thickness grid to 3D array
    %
    % Calculate Differences
    %
    cthkd_rf(:,:,i) = cthkr(:,:,i) - cthkf(:,:,i); % Compute RHO-FFE thickness differences
    cthkd_tf(:,:,i) = cthkt(:,:,i) - cthkf(:,:,i); % Compute T2S-FFE thickness differences
end
idx_cthkf = ~isnan(cthkf); % Identify non-NaN FFE thickness indices
idx_cthkr = ~isnan(cthkr); % Identify non-NaN RHO thickness indices
idx_cthkt = ~isnan(cthkt); % Identify non-NaN T2S thickness indices
cols = sum(max(idx_cthkf),3) + sum(max(idx_cthkr),3) + sum(max(idx_cthkt),3); % Sum valid columns across scans
rows = sum(max(idx_cthkf,[],2),3) + sum(max(idx_cthkr,[],2),3) + sum(max(idx_cthkt,[],2),3); % Sum valid rows across scans
cols = logical(cols); % Convert column sums to logical
rows = logical(rows); % Convert row sums to logical
col1 = find(cols,1,'first'); % Find first valid column
col2 = find(cols,1,'last'); % Find last valid column
row1 = find(rows,1,'first'); % Find first valid row
row2 = find(rows,1,'last'); % Find last valid row
cthkf = cthkf(row1:row2,col1:col2,:); % Crop FFE thickness array
cthkr = cthkr(row1:row2,col1:col2,:); % Crop RHO thickness array
cthkt = cthkt(row1:row2,col1:col2,:); % Crop T2S thickness array
cthkd_rf = cthkd_rf(row1:row2,col1:col2,:); % Crop RHO-FFE difference array
cthkd_tf = cthkd_tf(row1:row2,col1:col2,:); % Crop T2S-FFE difference array
tg = tg(row1:row2,col1:col2); % Crop theta grid
zg = zg(row1:row2,col1:col2); % Crop Z grid
[nr,nc] = size(tg); % Get dimensions of cropped grid
tg = tg(:); % Flatten theta grid to vector
zg = zg(:); % Flatten Z grid to vector
quad = quadconn(nr,nc); % Generate quadrilateral connectivity
zsc = t2sg.zsc; % Extract Z-scaling factors from T2S grid
t0 = -145; % Set theta cutoff angle
[idm_ant,idl_ant,idm_ctr,idl_ctr,idm_pos,idl_pos] = freg_axpf2_AD(tg,zg,zsc,t0,sd,div); % Get ROI indices
A = reshape(1:size(tg),nr,nc); % Create index matrix for grid points
A(~(idm_ant|idl_ant)) = NaN; % Set non-anterior points to NaN
A = min(A,[],2); % Get minimum valid indices per row
roi(1) = "fem_lat_pos.xlsx"; % Define lateral posterior ROI file
roi(2) = "fem_med_pos.xlsx"; % Define medial posterior ROI file
roi(3) = "fem_lat_ctr.xlsx"; % Define lateral central ROI file
roi(4) = "fem_med_ctr.xlsx"; % Define medial central ROI file
roi(5) = "fem_lat_ant.xlsx"; % Define lateral anterior ROI file
roi(6) = "fem_med_ant.xlsx"; % Define medial anterior ROI file
roi_idx(:,1) = idl_pos; % Assign lateral posterior ROI indices
roi_idx(:,2) = idm_pos; % Assign medial posterior ROI indices
roi_idx(:,3) = idl_ctr; % Assign lateral central ROI indices
roi_idx(:,4) = idm_ctr; % Assign medial central ROI indices
roi_idx(:,5) = idl_ant; % Assign lateral anterior ROI indices
roi_idx(:,6) = idm_ant; % Assign medial anterior ROI indices
coord = 1:size(tg); % Create coordinate indices
coord = coord(:); % Convert coordinates to column vector
ii=0; % Initialize row offset for Excel output
for i=1:ns % Loop through each subject
    %%
    % Turning the thicknesses into vectors
    cthkf_v = cthkf(:,:,i); % Extract FFE thickness for subject
    cthkr_v = cthkr(:,:,i); % Extract RHO thickness for subject
    cthkt_v = cthkt(:,:,i); % Extract T2S thickness for subject
    cthkf_v = cthkf_v(:); % Convert FFE thickness to vector
    cthkr_v = cthkr_v(:); % Convert RHO thickness to vector
    cthkt_v = cthkt_v(:); % Convert T2S thickness to vector
    % Get averages for specific ROIs
    avg_lats(1,1) = mean(cthkf_v(idl_pos),"omitnan"); % Compute mean FFE thickness for lateral posterior ROI
    avg_lats(2,1) = mean(cthkr_v(idl_pos),"omitnan"); % Compute mean RHO thickness for lateral posterior ROI
    avg_lats(3,1) = mean(cthkt_v(idl_pos),"omitnan"); % Compute mean T2S thickness for lateral posterior ROI
    avg_lats(1,2) = mean(cthkf_v(idl_ctr),"omitnan"); % Compute mean FFE thickness for lateral central ROI
    avg_lats(2,2) = mean(cthkr_v(idl_ctr),"omitnan"); % Compute mean RHO thickness for lateral central ROI
    avg_lats(3,2) = mean(cthkt_v(idl_ctr),"omitnan"); % Compute mean T2S thickness for lateral central ROI
    avg_lats(1,3) = mean(cthkf_v(idl_ant),"omitnan"); % Compute mean FFE thickness for lateral anterior ROI
    avg_lats(2,3) = mean(cthkr_v(idl_ant),"omitnan"); % Compute mean RHO thickness for lateral anterior ROI
    avg_lats(3,3) = mean(cthkt_v(idl_ant),"omitnan"); % Compute mean T2S thickness for lateral anterior ROI
    avg_meds(1,1) = mean(cthkf_v(idm_pos),"omitnan"); % Compute mean FFE thickness for medial posterior ROI
    avg_meds(2,1) = mean(cthkr_v(idm_pos),"omitnan"); % Compute mean RHO thickness for medial posterior ROI
    avg_meds(3,1) = mean(cthkt_v(idm_pos),"omitnan"); % Compute mean T2S thickness for medial posterior ROI
    avg_meds(1,2) = mean(cthkf_v(idm_ctr),"omitnan"); % Compute mean FFE thickness for medial central ROI
    avg_meds(2,2) = mean(cthkr_v(idm_ctr),"omitnan"); % Compute mean RHO thickness for medial central ROI
    avg_meds(3,2) = mean(cthkt_v(idm_ctr),"omitnan"); % Compute mean T2S thickness for medial central ROI
    avg_meds(1,3) = mean(cthkf_v(idm_ant),"omitnan"); % Compute mean FFE thickness for medial anterior ROI
    avg_meds(2,3) = mean(cthkr_v(idm_ant),"omitnan"); % Compute mean RHO thickness for medial anterior ROI
    avg_meds(3,3) = mean(cthkt_v(idm_ant),"omitnan"); % Compute mean T2S thickness for medial anterior ROI
    fstr=sd(i).FFE.femur.bfnam; % Get FFE bone file name
    f = figure; % Create new figure
    tiledlayout(2,6); % Set up 2x6 tiled layout
    sgtitle({[fstr(1:5) ' - Femur'];},'FontSize',16,'FontWeight','bold','Interpreter','none'); % Add supertitle with subject ID
    hf1=nexttile(1,[1 2]); % Create first tile spanning 2 columns
    cthk = cthkr_v; % Assign RHO thickness vector
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot RHO thickness patch
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    view(-90,90); % Set view to top-down
    title('T1\rho Cartilage Thicknesses','FontSize',16,'FontWeight','bold'); % Add title
    colorbar(hf1, 'Ticks',0:6); % Add colorbar with ticks
    clim([0 6]); % Set color limits
    xlim(hf1,[tmin tmax]); % Set X-axis limits
    ylim(hf1,[zmin zmax]); % Set Y-axis limits
    plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
    plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
    plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
    hf2=nexttile(3,[1 2]); % Create second tile spanning 2 columns
    cthk = cthkf_v; % Assign FFE thickness vector
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot FFE thickness patch
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    view(-90,90); % Set view to top-down
    title('T1FFE Cartilage Thicknesses','FontSize',16,'FontWeight','bold'); % Add title
    colorbar(hf2, 'Ticks',0:6); % Add colorbar with ticks
    clim([0 6]); % Set color limits
    xlim(hf2,[tmin tmax]); % Set X-axis limits
    ylim(hf2,[zmin zmax]); % Set Y-axis limits
    plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
    plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
    plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
    hf3=nexttile(5,[1 2]); % Create third tile spanning 2 columns
    cthk = cthkt_v; % Assign T2S thickness vector
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot T2S thickness patch
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    view(-90,90); % Set view to top-down
    title('T2S Cartilage Thicknesses','FontSize',16,'FontWeight','bold'); % Add title
    colorbar(hf3, 'Ticks',0:6); % Add colorbar with ticks
    clim([0 6]); % Set color limits
    xlim(hf3,[tmin tmax]); % Set X-axis limits
    ylim(hf3,[zmin zmax]); % Set Y-axis limits
    plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
    plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
    plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
    colormap(hf1, 'jet'); % Set colormap for RHO plot
    colormap(hf2, 'jet'); % Set colormap for FFE plot
    colormap(hf3, 'jet'); % Set colormap for T2S plot
    hf4=nexttile(8,[1 2]); % Create fourth tile spanning 2 columns
    dthkrf = cthkd_rf(:,:,i); % Extract RHO-FFE difference for subject
    patch(tg(quad'),zg(quad'),dthkrf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot RHO-FFE difference patch
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    view(-90,90); % Set view to top-down
    title('T1\rho - T1FFE Cartilage Thickness Differences','FontSize',16,'FontWeight','bold'); % Add title
    colorbar; % Add colorbar
    clim([-1.5 1.5]); % Set color limits
    xlim(hf4,[tmin tmax]); % Set X-axis limits
    ylim(hf4,[zmin zmax]); % Set Y-axis limits
    plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
    plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
    plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
    rflat1 = mean(dthkrf(idl_pos),"omitnan"); % Compute RHO-FFE lateral posterior mean difference
    rflat2 = mean(dthkrf(idl_ctr),"omitnan"); % Compute RHO-FFE lateral central mean difference
    rflat3 = mean(dthkrf(idl_ant),"omitnan"); % Compute RHO-FFE lateral anterior mean difference
    rfmed1 = mean(dthkrf(idm_pos),"omitnan"); % Compute RHO-FFE medial posterior mean difference
    rfmed2 = mean(dthkrf(idm_ctr),"omitnan"); % Compute RHO-FFE medial central mean difference
    rfmed3 = mean(dthkrf(idm_ant),"omitnan"); % Compute RHO-FFE medial anterior mean difference
    hf5=nexttile(10,[1 2]); % Create fifth tile spanning 2 columns
    dthktf = cthkd_tf(:,:,i); % Extract T2S-FFE difference for subject
    patch(tg(quad'),zg(quad'),dthktf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot T2S-FFE difference patch
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    view(-90,90); % Set view to top-down
    title('T2S - T1FFE Cartilage Thickness Differences','FontSize',16,'FontWeight','bold'); % Add title
    colorbar; % Add colorbar
    clim([-1.5 1.5]); % Set color limits
    xlim(hf5,[tmin tmax]); % Set X-axis limits
    ylim(hf5,[zmin zmax]); % Set Y-axis limits
    colormap(hf4, parula); % Set colormap for RHO-FFE difference
    colormap(hf5, parula); % Set colormap for T2S-FFE difference
    plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
    plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
    plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
    tflat1 = mean(dthktf(idl_pos),"omitnan"); % Compute T2S-FFE lateral posterior mean difference
    tflat2 = mean(dthktf(idl_ctr),"omitnan"); % Compute T2S-FFE lateral central mean difference
    tflat3 = mean(dthktf(idl_ant),"omitnan"); % Compute T2S-FFE lateral anterior mean difference
    tfmed1 = mean(dthktf(idm_pos),"omitnan"); % Compute T2S-FFE medial posterior mean difference
    tfmed2 = mean(dthktf(idm_ctr),"omitnan"); % Compute T2S-FFE medial central mean difference
    tfmed3 = mean(dthktf(idm_ant),"omitnan"); % Compute T2S-FFE medial anterior mean difference
    pic_nam=fullfile(rdir, 'Femoral Thickness Differences.pdf'); % Construct PDF file path
    if i==1 && exist(pic_nam,'file')==2 % Check if PDF exists for first subject
        delete(pic_nam); % Delete existing PDF
    end
    set(f,'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true); % Append plot to PDF
    close(f); % Close figure
    %%
    %
    % Statistics and Plot of Differences
    %
    dthkrf = dthkrf(:); % Convert RHO-FFE difference to vector
    dthktf = dthktf(:); % Convert T2S-FFE difference to vector
    dmeanrf = mean(dthkrf,"omitnan"); % Compute mean RHO-FFE difference
    dmeantf = mean(dthktf,"omitnan"); % Compute mean T2S-FFE difference
    dstdrf = std(dthkrf,"omitnan"); % Compute standard deviation of RHO-FFE differences
    dstdtf = std(dthktf,"omitnan"); % Compute standard deviation of T2S-FFE differences
    dmaxrf = max(dthkrf); % Get maximum RHO-FFE difference
    dminrf = min(dthkrf); % Get minimum RHO-FFE difference
    dmaxtf = max(dthktf); % Get maximum T2S-FFE difference
    dmintf = min(dthktf); % Get minimum T2S-FFE difference
    %
    f = figure; % Create new figure
    tiledlayout(2,1); % Set up 2x1 tiled layout
    sgtitle({[fstr(1:5) ' - Femur'];},'FontSize',16,'FontWeight','bold','Interpreter','none'); % Add supertitle with subject ID
    %
    nexttile; % Create first tile
    plot(dthkrf,'k.'); % Plot RHO-FFE differences
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    plot(axlim(1:2),[dmeanrf dmeanrf],'b--','LineWidth',1); % Plot mean RHO-FFE difference
    plot(axlim(1:2),[dmeanrf+3*dstdrf dmeanrf+3*dstdrf],'r--','LineWidth',1); % Plot upper 3 SD line
    plot(axlim(1:2),[dmeanrf-3*dstdrf dmeanrf-3*dstdrf],'r--','LineWidth',1); % Plot lower 3 SD line
    title('Rho/FFE Differences','FontSize',16,'FontWeight','bold'); % Add title
    nexttile; % Create second tile
    plot(dthktf,'k.'); % Plot T2S-FFE differences
    hold on; % Enable hold for multiple plots
    axlim = axis; % Get current axis limits
    plot(axlim(1:2),[dmeantf dmeantf],'b--','LineWidth',1); % Plot mean T2S-FFE difference
    plot(axlim(1:2),[dmeantf+3*dstdtf dmeantf+3*dstdtf],'r--','LineWidth',1); % Plot upper 3 SD line
    plot(axlim(1:2),[dmeantf-3*dstdtf dmeantf-3*dstdtf],'r--','LineWidth',1); % Plot lower 3 SD line
    title('T2S/FFE Differences','FontSize',16,'FontWeight','bold'); % Add title
    jj = 0; % Initialize counter for ROI divisions
    % Write Thicknesses to CSV Spreadsheet
    for k = 1:6 % Loop through ROIs
        output = fullfile(rdir,roi(k)); % Construct Excel file path
        if i==1 && exist(output,'file')==2 % Check if Excel file exists for first subject
            delete(output); % Delete existing Excel file
        end
        col_header1 = {'Subject'}; % Define subject column header
        col_header2 = {'ScanType'}; % Define scan type column header
        col_header3 = {'Compartment'}; % Define compartment column header
        col_header4 = {'Division'}; % Define division column header
        if rem(k,2) == 1 % Check if ROI is lateral
            jj = jj+1; % Increment division counter
            cmprt = 0; % Set compartment to lateral
            if rem(jj,3) == 1 % Check for posterior division
                division = 0; % Set division to posterior
            elseif rem(jj,3) == 2 % Check for central division
                division = 1; % Set division to central
            elseif rem(jj,3) == 0 % Check for anterior division
                division = 2; % Set division to anterior
            end
        elseif rem(k,2) == 0 % Check if ROI is medial
            cmprt = 1; % Set compartment to medial
            if rem(jj,3) == 1 % Check for posterior division
                division = 0; % Set division to posterior
            elseif rem(jj,3) == 2 % Check for central division
                division = 1; % Set division to central
            elseif rem(jj,3) == 0 % Check for anterior division
                division = 2; % Set division to anterior
            end
        end
        if i == 1 % Check if first subject
            if sum(roi_idx(:,k)) < 16385 % Verify ROI indices fit within Excel limits
                col_header5 = repmat('pt_', size(coord(roi_idx(:,k)),1),1); % Create point column headers
                col_header5 = [col_header5 int2str(coord(roi_idx(:,k)))]; % Append coordinate indices
                col_header5 = string(col_header5); % Convert to string array
                col_header5 = col_header5'; % Transpose header array
            else
                error([' *** ERROR: DATA EXCEEDS TOTAL NUMBER OF COLUMNS',' AVAILABLE IN EXCEL SHEET']); % Throw error if too many columns
            end
            writecell(col_header1,output,'Range','A1'); % Write subject header
            writecell(col_header2,output,'Range','B1'); % Write scan type header
            writecell(col_header3,output,'Range','C1'); % Write compartment header
            writecell(col_header4,output,'Range','D1'); % Write division header
            writematrix(col_header5,output, 'Range','E1'); % Write point headers
        end
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+1+ii)]); % Write subject ID for FFE
        writematrix(0,output,'Range',['B' int2str(i+1+ii)]); % Write FFE scan type
        writematrix(cmprt,output,'Range',['C' int2str(i+1+ii)]); % Write compartment
        writematrix(division,output,'Range',['D' int2str(i+1+ii)]); % Write division
        writematrix(cthkf_v(roi_idx(:,k))',output,'Range',['E' int2str(i+1+ii)]); % Write FFE thicknesses
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+2+ii)]); % Write subject ID for RHO
        writematrix(1,output,'Range',['B' int2str(i+2+ii)]); % Write RHO scan type
        writematrix(cmprt,output,'Range',['C' int2str(i+2+ii)]); % Write compartment
        writematrix(division,output,'Range',['D' int2str(i+2+ii)]); % Write division
        writematrix(cthkr_v(roi_idx(:,k))',output,'Range',['E' int2str(i+2+ii)]); % Write RHO thicknesses
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+3+ii)]); % Write subject ID for T2S
        writematrix(2,output,'Range',['B' int2str(i+3+ii)]); % Write T2S scan type
        writematrix(cmprt,output,'Range',['C' int2str(i+3+ii)]); % Write compartment
        writematrix(division,output,'Range',['D' int2str(i+3+ii)]); % Write division
        writematrix(cthkt_v(roi_idx(:,k))',output,'Range',['E' int2str(i+3+ii)]); % Write T2S thicknesses
        %
    end
    set(f,'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true); % Append plot to PDF
    close(f); % Close figure
    %
    ii = ii+2; % Increment Excel row offset
end
%%
n_rf = sum(~isnan(cthkd_rf),3); % Count valid RHO-FFE differences
idxn = find(n_rf<3); % Find indices with fewer than 3 subjects
n_rf(idxn) = NaN; % Set low-count indices to NaN
avg_rf = mean(cthkd_rf,3,"omitnan"); % Compute mean RHO-FFE differences
avg_rf(idxn) = NaN; % Set low-count means to NaN
std_rf = std(cthkd_rf,0,3,"omitnan"); % Compute standard deviation of RHO-FFE differences
std_rf(idxn) = NaN; % Set low-count standard deviations to NaN
n_tf = sum(~isnan(cthkd_tf),3); % Count valid T2S-FFE differences
idxn = find(n_tf<3); % Find indices with fewer than 3 subjects
n_tf(idxn) = NaN; % Set low-count indices to NaN
avg_tf = mean(cthkd_tf,3,"omitnan"); % Compute mean T2S-FFE differences
avg_tf(idxn) = NaN; % Set low-count means to NaN
std_tf = std(cthkd_tf,0,3,"omitnan"); % Compute standard deviation of T2S-FFE differences
std_tf(idxn) = NaN; % Set low-count standard deviations to NaN
rflat1 = mean(avg_rf(idl_pos),"omitnan"); % Compute mean RHO-FFE lateral posterior difference
rflat2 = mean(avg_rf(idl_ctr),"omitnan"); % Compute mean RHO-FFE lateral central difference
rflat3 = mean(avg_rf(idl_ant),"omitnan"); % Compute mean RHO-FFE lateral anterior difference
tflat1 = mean(avg_tf(idl_pos),"omitnan"); % Compute mean T2S-FFE lateral posterior difference
tflat2 = mean(avg_tf(idl_ctr),"omitnan"); % Compute mean T2S-FFE lateral central difference
tflat3 = mean(avg_tf(idl_ant),"omitnan"); % Compute mean T2S-FFE lateral anterior difference
rfmed1 = mean(avg_rf(idm_pos),"omitnan"); % Compute mean RHO-FFE medial posterior difference
rfmed2 = mean(avg_rf(idm_ctr),"omitnan"); % Compute mean RHO-FFE medial central difference
rfmed3 = mean(avg_rf(idm_ant),"omitnan"); % Compute mean RHO-FFE medial anterior difference
tfmed1 = mean(avg_tf(idm_pos),"omitnan"); % Compute mean T2S-FFE medial posterior difference
tfmed2 = mean(avg_tf(idm_ctr),"omitnan"); % Compute mean T2S-FFE medial central difference
tfmed3 = mean(avg_tf(idm_ant),"omitnan"); % Compute mean T2S-FFE medial anterior difference
% 3 - 10 Colormap for Number of Subjects
cmap=[0 0.8 1; 0.2 1 0.8; 0.6 1 0.4; 1 1 0; 1 0.6 0; 1 0.4 0;1 0 0; 0.6 0 0]; % Define colormap for subject counts
f = figure; % Create new figure
t = tiledlayout(3,2); % Set up 3x2 tiled layout
title(t,'Femoral Cartilage Thickness Differences','FontSize',20,'FontWeight','bold'); % Add layout title
hf1 = nexttile(1); % Create first tile
patch(tg(quad'),zg(quad'),avg_rf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot mean RHO-FFE differences
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
colormap(hf1,parula); % Set colormap
hf1.CLim = [-1.5 1.5]; % Set color limits
colorbar; % Add colorbar
xlim(hf1,[tmin tmax]); % Set X-axis limits
ylim(hf1,[zmin zmax]); % Set Y-axis limits
hf1.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
axis off; % Hide axes
box off; % Hide box
hf2 = nexttile(2); % Create second tile
patch(tg(quad'),zg(quad'),avg_tf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot mean T2S-FFE differences
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
title('T2S - T1FFE Averages','FontSize',14,'FontWeight','bold'); % Add title
colormap(hf2,parula); % Set colormap
hf2.CLim = [-1.5 1.5]; % Set color limits
colorbar; % Add colorbar
xlim(hf2,[tmin tmax]); % Set X-axis limits
ylim(hf2,[zmin zmax]); % Set Y-axis limits
hf2.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
axis off; % Hide axes
box off; % Hide box
hf3 = nexttile(3); % Create third tile
patch(tg(quad'),zg(quad'),std_rf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot RHO-FFE standard deviations
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
title('Rho - T1FFE','Standard Deviatons','FontSize',14,'FontWeight','bold'); % Add title
hf3.CLim = [0 1.5]; % Set color limits
colormap(hf3,jet); % Set colormap
colorbar(hf3,'Ticks', 0:.25:1.5); % Add colorbar with ticks
xlim(hf3,[tmin tmax]); % Set X-axis limits
ylim(hf3,[zmin zmax]); % Set Y-axis limits
hf3.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
hf4 = nexttile(4); % Create fourth tile
patch(tg(quad'),zg(quad'),std_tf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot T2S-FFE standard deviations
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
title('T2S - T1FFE','Standard Deviatons','FontSize',14,'FontWeight','bold'); % Add title
hf4.CLim = [0 1.5]; % Set color limits
colormap(hf4,jet); % Set colormap
colorbar(hf4,'Ticks', 0:.25:1.5); % Add colorbar with ticks
xlim(hf4,[tmin tmax]); % Set X-axis limits
ylim(hf4,[zmin zmax]); % Set Y-axis limits
hf4.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
hf5 = nexttile(5); % Create fifth tile
patch(tg(quad'),zg(quad'),n_rf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot RHO-FFE subject counts
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
title('Rho - T1FFE','Subject Numbers','FontSize',14,'FontWeight','bold'); % Add title
hf5.CLim=[2 10]; % Set color limits
colormap(hf5,cmap); % Set colormap
colorbar(hf5,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0); % Add colorbar with custom labels
xlim(hf5,[tmin tmax]); % Set X-axis limits
ylim(hf5,[zmin zmax]); % Set Y-axis limits
hf5.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
hf6 = nexttile(6); % Create sixth tile
patch(tg(quad'),zg(quad'),n_tf(quad'),'FaceColor','interp','EdgeColor','interp'); % Plot T2S-FFE subject counts
hold on; % Enable hold for multiple plots
view(-90,90); % Set view to top-down
axlim = axis; % Get current axis limits
title('T2S - T1FFE','Subject Numbers','FontSize',14,'FontWeight','bold'); % Add title
hf6.CLim=[2 10]; % Set color limits
colormap(hf6,cmap); % Set colormap
colorbar(hf6,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0); % Add colorbar with custom labels
xlim(hf6,[-210 30]); % Set X-axis limits
ylim(hf6,[zmin zmax]); % Set Y-axis limits
hf6.PlotBoxAspectRatio = [0.9855 1 0.9855]; % Set aspect ratio
plot([-145 -145],[axlim(3) axlim(4)],"k",'LineWidth',1.5); % Plot vertical black line at -145
plot([axlim(1) axlim(2)],[0 0],"k",'LineWidth',1.5); % Plot horizontal black line at Z=0
plot(tg(A)-1,zg(A),"k",'LineWidth',1.5); % Plot anterior boundary
set(f,'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
exportgraphics(f,pic_nam, "Resolution", 300, 'Append', true); % Append plot to PDF
close(f); % Close figure
return % Exit the script