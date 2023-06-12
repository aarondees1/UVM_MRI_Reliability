%#######################################################################
%
%           * Tibia Cartilage THicKness Comparison Program *
%
%          M-File which reads the T1FFE and T1rho cartilage thicknesses,
%     finds a combined grid that covers both data sets and finds the
%     differences in cartilage thicknesses.
%
%          Plots, Outputs, ...
%
%     NOTES:  1.  Both grids must have integer coordinates.
%
%     17-Mar-2022 * Mack Gardner-Morse
%

%#######################################################################
%%
clear;
close all;
clc;
    grayColor = [.7 .7 .7]; 

div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
load(fullfile(rdir,'Subject Files Full.mat'));
ns=size(sd,1);
%mat_rho = fullfile(ddir_rho,'tcart08_1.mat');
%for i=1:ns


% %%
% %
% % Get Analysis Grids
% %
ffeg = load(fullfile(gdir,'FFE_tgrid08_.mat'));
rhog = load(fullfile(gdir,'RHO_tgrid08_.mat'));
t2sg = load(fullfile(gdir,'T2S_tgrid08_.mat'));
% rhog = load(fullfile('TibiaCartThk_1_11Mar2022(T1Rho)', ...
%             'tgrid08_1.mat'));
%%%
%
% Get Coordinates and Ranges for Lateral Compartment
%
xql_ffe = ffeg.xql;
yql_ffe = ffeg.yql;
%
xql_rho = rhog.xql;
yql_rho = rhog.yql;
%
xql_t2s = t2sg.xql;
yql_t2s = t2sg.yql;
%
xmin_ffe = min(xql_ffe);
xmax_ffe = max(xql_ffe);
nc_ffe = ffeg.nxl;         % Number of columns
%
xmin_rho = min(xql_rho);
xmax_rho = max(xql_rho);
nc_rho = rhog.nxl;          % Number of columns
%
xmin_t2s = min(xql_t2s);
xmax_t2s = max(xql_t2s);
nc_t2s = t2sg.nxl;
%
ymin_ffe = min(yql_ffe);
ymax_ffe = max(yql_ffe);
nr_ffe = ffeg.nyl;         % Number of rows

ymin_rho = min(yql_rho);
ymax_rho = max(yql_rho);
nr_rho = rhog.nyl;         % Number of rows

ymin_t2s = min(yql_t2s);
ymax_t2s = max(yql_t2s);
nr_t2s = t2sg.nyl;         % Number of rows
%%

%Get Combined Grid for Both Data Sets

xmin=[xmin_ffe; xmin_rho; xmin_t2s];
xmin=min(xmin);

xmax = [xmax_ffe; xmax_rho; xmax_t2s];
xmax=max(xmax);

ymin = [ymin_ffe; ymin_rho; ymin_t2s];
ymin = min(ymin);

ymax = [ymax_ffe; ymax_rho; ymax_t2s];
ymax = max(ymax);


[xgl,ygl] = meshgrid(xmin:xmax,ymin:ymax);
[nlr,nlc] = size(xgl);
xgl = xgl(:);
ygl = ygl(:);
quadl = quadconn(nlr,nlc);   % Quadrilateral connectivity for grid
%%
%
% Get Indexes into Combined Grid
%
nl = nlr*nlc;           % Number of points in combined lateral grid
idx = reshape(1:nl,nlr,nlc);
%
offstcf = round(xmin_ffe-xmin+1);      % Differences in integer grid == column index
offstrf = round(ymin_ffe-ymin+1);      % Differences in integer grid == row index
idxfl = idx(offstrf:offstrf+nr_ffe-1,offstcf:offstcf+nc_ffe-1);   % Index for T1FFE
%
offstcr = round(xmin_rho-xmin+1);      % Differences in integer grid == column index
offstrr = round(ymin_rho-ymin+1);      % Differences in integer grid == row index
idxrl  = idx(offstrr:offstrr+nr_rho-1,offstcr:offstcr+nc_rho-1);   % Index for T1rho
%
offstct = round(xmin_t2s-xmin+1);      % Differences in integer grid == column index
offstrt = round(ymin_t2s-ymin+1);      % Differences in integer grid == row index
idxtl = idx(offstrt:offstrt+nr_t2s-1,offstct:offstct+nc_t2s-1);   % Index for T2*

%%
%
% Get Coordinates and Ranges for Medial Compartment
%
xqm_ffe = ffeg.xqm;
yqm_ffe = ffeg.yqm;
%
xqm_rho = rhog.xqm;
yqm_rho = rhog.yqm;
%
xqm_t2s = t2sg.xqm;
yqm_t2s = t2sg.yqm;
%
xmin_ffe = min(xqm_ffe);
xmax_ffe = max(xqm_ffe);
nc_ffe = ffeg.nxm;         % Number of columns
xmin_rho = min(xqm_rho);
xmax_rho = max(xqm_rho);
nc_rho = rhog.nxm;         % Number of columns
xmin_t2s = min(xqm_t2s);
xmax_t2s = max(xqm_t2s);
nc_t2s = t2sg.nxm;         % Number of columns
%
ymin_ffe = min(yqm_ffe);
ymax_ffe = max(yqm_ffe);
nr_ffe = ffeg.nym;         % Number of rows
ymin_rho = min(yqm_rho);
ymax_rho = max(yqm_rho);
nr_rho = rhog.nym;         % Number of rows
ymin_t2s = min(yqm_t2s);
ymax_t2s = max(yqm_t2s);
nr_t2s = t2sg.nym;         % Number of rows

%%
%
% Get Combined Grid for Both Data Sets
%
xmin=[xmin_ffe; xmin_rho; xmin_t2s];
xmin=min(xmin);

xmax = [xmax_ffe; xmax_rho; xmax_t2s];
xmax=max(xmax);

ymin = [ymin_ffe; ymin_rho; ymin_t2s];
ymin = min(ymin);

ymax = [ymax_ffe; ymax_rho; ymax_t2s];
ymax = max(ymax);
%
[xgm,ygm] = meshgrid(xmin:xmax,ymin:ymax);
[nmr,nmc] = size(xgm);
xgm = xgm(:);
ygm = ygm(:);
quadm = quadconn(nmr,nmc);   % Quadrilateral connectivity for grid

%%
%
% Get Indexes into Combined Grid
%
nm = nmr*nmc;           % Number of points in combined medial grid
idx = reshape(1:nm,nmr,nmc);
%
offstcf = round(xmin_ffe-xmin+1);      % Differences in integer grid == column index
offstrf = round(ymin_ffe-ymin+1);      % Differences in integer grid == row index
idxfm = idx(offstrf:offstrf+nr_ffe-1,offstcf:offstcf+nc_ffe-1);   % Index for T1FFE
%
offstcr = round(xmin_rho-xmin+1);      % Differences in integer grid == column index
offstrr = round(ymin_rho-ymin+1);      % Differences in integer grid == row index
idxrm = idx(offstrr:offstrr+nr_rho-1,offstcr:offstcr+nc_rho-1);   % Index for T1rho

offstct = round(xmin_t2s-xmin+1);      % Differences in integer grid == column index
offstrt = round(ymin_t2s-ymin+1);      % Differences in integer grid == row index
idxtm = idx(offstrt:offstrt+nr_t2s-1,offstct:offstct+nc_t2s-1);   % Index for T2*

%%
%
% Read Cartilage Thicknesses
%

cthkld_rf = zeros(nlr,nlc,ns);
cthkld_tf = zeros(nlr,nlc,ns);
cthkmd_rf = zeros(nmr,nmc,ns);
cthkmd_tf = zeros(nmr,nmc,ns);

bonemaxsl=zeros(ns,3,3);
bonemaxsm=zeros(ns,3,3);
boneminsl=zeros(ns,3,3);
boneminsm=zeros(ns,3,3);

% j == 1 FFE, j == 2 RHO, j == 3 T2S
%%
for i=1:ns
    bonemaxsl(i,:,1) = max(sd(i).FFE.tibia.bone.lateral.points);
    bonemaxsm(i,:,1) = max(sd(i).FFE.tibia.bone.medial.points);
    boneminsl(i,:,1) = min(sd(i).FFE.tibia.bone.lateral.points);
    boneminsm(i,:,1) = min(sd(i).FFE.tibia.bone.medial.points);

    bonemaxsl(i,:,2) = max(sd(i).RHO.tibia.bone.lateral.points);
    bonemaxsm(i,:,2) = max(sd(i).RHO.tibia.bone.medial.points);
    boneminsl(i,:,2) = min(sd(i).RHO.tibia.bone.lateral.points);
    boneminsm(i,:,2) = min(sd(i).RHO.tibia.bone.medial.points);

    bonemaxsl(i,:,3) = max(sd(i).T2S.tibia.bone.lateral.points);
    bonemaxsm(i,:,3) = max(sd(i).T2S.tibia.bone.medial.points);
    boneminsl(i,:,3) = min(sd(i).T2S.tibia.bone.lateral.points);
    boneminsm(i,:,3) = min(sd(i).T2S.tibia.bone.medial.points);
end
bonerangel = bonemaxsl - boneminsl;
bonerangem = bonemaxsm - boneminsm;

thirdsl=bonerangel/3;
thirdsm=bonerangem/3;

reg_lines1l = boneminsl + thirdsl;
reg_lines1m = boneminsm + thirdsm;
reg_lines2l = reg_lines1l + thirdsl;
reg_lines2m = reg_lines1m + thirdsm;


reg_lat1 = squeeze(mean(reg_lines1l(:,1,:))); % Average over subjects
reg_lat1_all=mean(reg_lat1);          % Average over scans
%
reg_lat2 = squeeze(mean(reg_lines2l(:,1,:))); % Average over subjects
reg_lat2_all=mean(reg_lat2);          % Average over scans
%
reg_med1 = squeeze(mean(reg_lines1m(:,1,:))); % Average over subjects
reg_med1_all = mean(reg_med1);          % Average over scans
%
reg_med2 = squeeze(mean(reg_lines2m(:,1,:))); % Average over subjects
reg_med2_all = mean(reg_med2);          % Average over scans


%%
for i=1:ns
    fstr=sd(i).FFE.tibia.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_tcart08_thk.mat'];
    ffe = load(fullfile(tdir,fstr));

    fstr=sd(i).RHO.tibia.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_tcart08_thk.mat'];
    rho = load(fullfile(tdir,fstr));

    fstr=sd(i).T2S.tibia.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_tcart08_thk.mat'];
    t2s = load(fullfile(tdir,fstr));

    cthkfl = ffe.cthkl;
    cthkrl = rho.cthkl;
    cthktl= t2s.cthkl;
    %%
    %
    % Get Cartilage Thicknesses at Similar Coordinates
    %
    cthklf = NaN(nlr,nlc);  % NaN == missing data
    cthklr = NaN(nlr,nlc);
    cthklt = NaN(nlr,nlc);
    %
    cthklf(idxfl) = cthkfl;
    cthklr(idxrl) = cthkrl;
    cthklt(idxtl) = cthktl;


    %
    % Read Cartilage Thicknesses
    %
    cthkfm = ffe.cthkm;
    cthkrm = rho.cthkm;
    cthktm = t2s.cthkm;

    %%
    %
    % Get Cartilage Thicknesses at Similar Coordinates
    %
    cthkmf = NaN(nmr,nmc);  % NaN == missing data
    cthkmr = NaN(nmr,nmc);
    cthkmt = NaN(nmr,nmc);
    %
    cthkmf(idxfm) = cthkfm;
    cthkmr(idxrm) = cthkrm;
    cthkmt(idxtm) = cthktm;

    %
    % Calculate Differences
    %
    cthkld_rf(:,:,i) = cthklr-cthklf;
    cthkld_tf(:,:,i) = cthklt-cthklf;
    cthkmd_rf(:,:,i) = cthkmr-cthkmf;
    cthkmd_tf(:,:,i) = cthkmt-cthkmf;

    idvl_rf = ~isnan(cthkld_rf(:,:,i));  % Valid differences
    idvm_rf = ~isnan(cthkmd_rf(:,:,i));  % Valid differences
    idvl_tf = ~isnan(cthkld_tf(:,:,i));  % Valid differences
    idvm_tf = ~isnan(cthkmd_tf(:,:,i));  % Valid differences

    %%

    %%
    %
    % Plots of Cartilage Thicknesses and Thickness Differences
    %
    figure;
    tiledlayout(2,6);
    sgtitle({[fstr(1:5) ' - Tibia'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;

    %
    hf1=nexttile(1,[1 2]);
    cthk = cthklr;
    %cthk(~idvl_rf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmr;
    %cthk(~idvm_rf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
        axis equal;
    view(-90,90);
    title('T1\rho Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar;
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);    
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1(2) reg_lat1(2)],[axlim(4) 0],"k");
    plot([reg_lat2(2) reg_lat2(2)],[axlim(4) 0],"k");
    plot([reg_med1(2) reg_med1(2)],[0 axlim(3)],"k");
    plot([reg_med2(2) reg_med2(2)],[0 axlim(3)],"k");
    %
    hf2=nexttile(3,[1 2]);
    cthk = cthklt;
    %cthk(~idvl_tf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmt;
    %cthk(~idvm_tf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
        axis equal;
    view(-90,90);
    title('T2S Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar;
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);
        scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1(3) reg_lat1(3)],[axlim(4) 0],"k");
    plot([reg_lat2(3) reg_lat2(3)],[axlim(4) 0],"k");
    plot([reg_med1(3) reg_med1(3)],[0 axlim(3)],"k");
    plot([reg_med2(3) reg_med2(3)],[0 axlim(3)],"k");
    %
    hf3=nexttile(5,[1 2]);
    cthk = cthklf;
    %cthk(~idvl_rf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmf;
    %cthk(~idvm_rf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    thk_max=max(max(cthk(quadl')));

    disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
        axis equal;
    view(-90,90);
    title('T1FFE Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar;
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);
        scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1(1) reg_lat1(1)],[axlim(4) 0],"k");
    plot([reg_lat2(1) reg_lat2(1)],[axlim(4) 0],"k");
    plot([reg_med1(1) reg_med1(1)],[0 axlim(3)],"k");
    plot([reg_med2(1) reg_med2(1)],[0 axlim(3)],"k");
    %
    colormap(hf1, 'jet(12)');
    colormap(hf2, 'jet(12)');
    colormap(hf3, 'jet(12)');
    %
    hf4=nexttile(7,[1 3]);
    dthklrf = cthkld_rf(:,:,i);
    patch(xgl(quadl'),ygl(quadl'),dthklrf(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    dthkmrf = cthkmd_rf(:,:,i);
    patch(xgm(quadm'),ygm(quadm'),dthkmrf(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
        axis equal;
    view(-90,90);
    title('T1\rho - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colorbar;
    clim([-2.5 2.5]);
    xlim([-30 30]);
    ylim([-40 40]);
        scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
        axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    %
    hf5=nexttile(10,[1 3]);
    dthkltf = cthkld_tf(:,:,i);
    patch(xgl(quadl'),ygl(quadl'),dthkltf(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    dthkmtf = cthkmd_tf(:,:,i);
    patch(xgm(quadm'),ygm(quadm'),dthkmtf(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    axis equal;
    view(-90,90);
    title('T2S - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colorbar;
    clim([-2.5 2.5]);
    xlim([-30 30]);
    ylim([-40 40]);
    colormap(hf4, 'parula(10)');
    colormap(hf5, 'parula(10)');
        scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
        axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    %
    pic_nam=fullfile(rdir,[fstr(1:6) 'tcthk_cmp.ps']);
    print('-dpsc2', '-r600', '-fillpage',pic_nam);


    %%
    %
    % Statistics and Plot of Differences
    %
    dthklrf = dthklrf(idvl_rf(:));
    dthkmrf = dthkmrf(idvm_rf(:));
    dthkltf = dthkltf(idvl_tf(:));
    dthkmtf = dthkmtf(idvm_tf(:));
    dmeanlrf = mean(dthklrf);
    dmeanmrf = mean(dthkmrf);
    dmeanltf = mean(dthkltf);
    dmeanmtf = mean(dthkmtf);
    dstdlrf = std(dthklrf);
    dstdmrf = std(dthkmrf);
    dstdltf = std(dthkltf);
    dstdmtf = std(dthkmtf);

    dmaxlrf = max(dthklrf);
    dminlrf = min(dthklrf);
    dmaxmrf = max(dthkmrf);
    dminmrf = min(dthkmrf);
    dmaxltf = max(dthkltf);
    dminltf = min(dthkltf);
    dmaxmtf = max(dthkmtf);
    dminmtf = min(dthkmtf);

    coordl_rf=find(idvl_rf);
    coordl_tf=find(idvl_tf);
    coordm_rf=find(idvm_rf);
    coordm_tf=find(idvm_tf);
    %

    %
    figure;
    sgtitle({[fstr(1:5) ' - Tibia'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;
    %
    subplot(2,2,1);
    plot(dthklrf,'k.');
    hold on;
    axlim = axis;
    plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanlrf dmeanlrf],'k-','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf+3*dstdlrf dmeanlrf+3*dstdlrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf-3*dstdlrf dmeanlrf-3*dstdlrf],'r--','LineWidth',1);
    title('Rho/FFE Lateral Differences', ...
        'FontSize',16,'FontWeight','bold');
    %
    subplot(2,2,3);
    plot(dthkmrf,'k.');
    hold on;
    axlim = axis;
    plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanmrf dmeanmrf],'k-','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf+3*dstdmrf dmeanmrf+3*dstdmrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf-3*dstdmrf dmeanmrf-3*dstdmrf],'r--','LineWidth',1);
    title('Rho/FFE Medial Differences', ...
        'FontSize',16,'FontWeight','bold');

    subplot(2,2,2);
    plot(dthkltf,'k.');
    hold on;
    axlim = axis;
    plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanltf dmeanltf],'k-','LineWidth',1);
    plot(axlim(1:2),[dmeanltf+3*dstdltf dmeanltf+3*dstdltf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanltf-3*dstdltf dmeanltf-3*dstdltf],'r--','LineWidth',1);
    title('T2S/FFE Lateral Differences', ...
        'FontSize',16,'FontWeight','bold');
    %
    subplot(2,2,4);
    plot(dthkmtf,'k.');
    hold on;
    axlim = axis;
    plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanmtf dmeanmtf],'k-','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf+3*dstdmtf dmeanmtf+3*dstdmtf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf-3*dstdmtf dmeanmtf-3*dstdmtf],'r--','LineWidth',1);
    title('T2S/FFE Medial Differences', ...
        'FontSize',16,'FontWeight','bold');
    

    % Write Differences to CSV Spreadsheet
    output = fullfile(rdir,[fstr(1:6) 'Thickness_differences.xlsx']);
    col_header = {'Thickness Diff (mm)','Grid Coordinates','Mean','Std','Max','Min'};

    writematrix(dthklrf,output,'Sheet','Lateral RhoFFE','Range','A2')
    writematrix(coordl_rf,output,'Sheet','Lateral RhoFFE','Range','B2')
    writematrix(dmeanlrf,output,'Sheet','Lateral RhoFFE','Range','C2')
    writematrix(dstdlrf,output,'Sheet','Lateral RhoFFE','Range','D2')
    writematrix(dmaxlrf,output,'Sheet','Lateral RhoFFE','Range','E2')
    writematrix(dminlrf,output,'Sheet','Lateral RhoFFE','Range','F2')
    writecell(col_header,output,'Sheet','Lateral RhoFFE','Range','A1')

    writematrix(dthkltf,output,'Sheet','Lateral T2SFFE','Range','A2')
    writematrix(coordl_tf,output,'Sheet','Lateral T2SFFE','Range','B2')
    writematrix(dmeanltf,output,'Sheet','Lateral T2SFFE','Range','C2')
    writematrix(dstdltf,output,'Sheet','Lateral T2SFFE','Range','D2')
    writematrix(dmaxltf,output,'Sheet','Lateral T2SFFE','Range','E2')
    writematrix(dminltf,output,'Sheet','Lateral T2SFFE','Range','F2')
    writecell(col_header,output,'Sheet','Lateral T2SFFE','Range','A1')

    writematrix(dthkmrf,output,'Sheet','Medial RhoFFE','Range','A2')
    writematrix(coordm_rf,output,'Sheet','Medial RhoFFE','Range','B2')
    writematrix(dmeanmrf,output,'Sheet','Medial RhoFFE','Range','C2')
    writematrix(dstdmrf,output,'Sheet','Medial RhoFFE','Range','D2')
    writematrix(dmaxmrf,output,'Sheet','Medial RhoFFE','Range','E2')
    writematrix(dminmrf,output,'Sheet','Medial RhoFFE','Range','F2')
    writecell(col_header,output,'Sheet','Medial RhoFFE','Range','A1')

    writematrix(dthkmtf,output,'Sheet','Medial T2SFFE','Range','A2')
    writematrix(coordm_tf,output,'Sheet','Medial T2SFFE','Range','B2')
    writematrix(dmeanmtf,output,'Sheet','Medial T2SFFE','Range','C2')
    writematrix(dstdmtf,output,'Sheet','Medial T2SFFE','Range','D2')
    writematrix(dmaxmtf,output,'Sheet','Medial T2SFFE','Range','E2')
    writematrix(dminmtf,output,'Sheet','Medial T2SFFE','Range','F2')
    writecell(col_header,output,'Sheet','Medial T2SFFE','Range','A1')
    %
    print('-dpsc2', '-r600', '-fillpage', '-append',pic_nam);
    %
end
    avg_rflat = mean(cthkld_rf,3,"omitnan");
    n_rflat = sum(~isnan(cthkld_rf),3);
    avg_tflat = mean(cthkld_tf,3,"omitnan");
    n_tflat = sum(~isnan(cthkld_tf),3);
    avg_rfmed = mean(cthkmd_rf,3,"omitnan");
    n_rfmed = sum(~isnan(cthkmd_rf),3);
    avg_tfmed = mean(cthkmd_tf,3,"omitnan");
    n_tfmed = sum(~isnan(cthkmd_tf),3);
%%
cmap=[0 0 0.6; 0 0 1; 0 0.4 1; 0 0.8 1; 0.2 1 0.8; 0.6 1 0.4;
1 1 0; 1 0.6 0; 1 0.4 0; 1 0 0; 0.6 0 0];

newtickVals=0:10;

    figure;
    tiledlayout(2,2);
        sgtitle({'Averages - Tibia';}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;

    hf1=nexttile;
    patch(xgl(quadl'),ygl(quadl'),avg_rflat(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    patch(xgm(quadm'),ygm(quadm'),avg_rfmed(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    view(-90,90);
     axis equal;
    title('Average Rho - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colormap(hf1,parula(9));
    hf1.CLim = [-2 2];
    colorbar;
    xlim(hf1,[-30 30]);
    ylim(hf1,[-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
   
        

    hf2=nexttile;
    patch(xgl(quadl'),ygl(quadl'),avg_tflat(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    patch(xgm(quadm'),ygm(quadm'),avg_tfmed(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    view(-90,90);
     axis equal;
    title('Average T2S - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colormap(hf2,parula(9));
    hf2.CLim = [-2 2];
    colorbar;
    xlim(hf2,[-30 30]);
    ylim(hf2,[-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
 

    hf3 = nexttile;
    patch(xgl(quadl'),ygl(quadl'),n_rflat(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    patch(xgm(quadm'),ygm(quadm'),n_rfmed(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    view(-90,90);
     axis equal;
        title('Average Rho - T1FFE Subject Numbers', ...
        'FontSize',16,'FontWeight','bold');
    hf3.CLim=[0 11];
    colormap(hf3,cmap);
    colorbar(hf3,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
    xlim(hf3,[-30 30]);
    ylim(hf3,[-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);

    hf4 = nexttile;
    patch(xgl(quadl'),ygl(quadl'),n_tflat(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    patch(xgm(quadm'),ygm(quadm'),n_tfmed(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    view(-90,90);
     axis equal;
            title('Average T2S - T1FFE Subject Numbers', ...
        'FontSize',16,'FontWeight','bold');
    hf4.CLim=[0 11];
    colormap(hf4,cmap);
    colorbar(hf4,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
    xlim(hf4,[-30 30]);
    ylim(hf4,[-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
 
    print('-dpsc2', '-r600', '-fillpage', fullfile(rdir,'Average Thicknesses.ps'));

    coordl_rf = find(~isnan(avg_rflat));
    coordl_tf = find(~isnan(avg_tflat));
    coordm_rf = find(~isnan(avg_rfmed));
    coordm_tf = find(~isnan(avg_tfmed));

    output = fullfile(rdir,'Average Thickness_differences.xlsx');
    col_header = {'Thickness Diff (mm)','Grid Coordinates','Mean','Std','Max','Min'};

    writematrix(avg_rflat,output,'Sheet','Lateral RhoFFE','Range','A2')
    writematrix(coordl_rf,output,'Sheet','Lateral RhoFFE','Range','B2')
    writematrix(mean(avg_rflat),output,'Sheet','Lateral RhoFFE','Range','C2')
    writematrix(std(avg_rflat),output,'Sheet','Lateral RhoFFE','Range','D2')
    writematrix(max(avg_rflat),output,'Sheet','Lateral RhoFFE','Range','E2')
    writematrix(min(avg_rflat),output,'Sheet','Lateral RhoFFE','Range','F2')
    writecell(col_header,output,'Sheet','Lateral RhoFFE','Range','A1')

    writematrix(avg_tflat,output,'Sheet','Lateral T2SFFE','Range','A2')
    writematrix(coordl_tf,output,'Sheet','Lateral T2SFFE','Range','B2')
    writematrix(mean(avg_tflat),output,'Sheet','Lateral T2SFFE','Range','C2')
    writematrix(std(avg_tflat),output,'Sheet','Lateral T2SFFE','Range','D2')
    writematrix(max(avg_tflat),output,'Sheet','Lateral T2SFFE','Range','E2')
    writematrix(min(avg_tflat),output,'Sheet','Lateral T2SFFE','Range','F2')
    writecell(col_header,output,'Sheet','Lateral T2SFFE','Range','A1')

    writematrix(avg_rfmed,output,'Sheet','Medial RhoFFE','Range','A2')
    writematrix(coordm_rf,output,'Sheet','Medial RhoFFE','Range','B2')
    writematrix(mean(avg_rfmed),output,'Sheet','Medial RhoFFE','Range','C2')
    writematrix(std(avg_rfmed),output,'Sheet','Medial RhoFFE','Range','D2')
    writematrix(max(avg_rfmed),output,'Sheet','Medial RhoFFE','Range','E2')
    writematrix(min(avg_rfmed),output,'Sheet','Medial RhoFFE','Range','F2')
    writecell(col_header,output,'Sheet','Medial RhoFFE','Range','A1')

    writematrix(avg_tfmed,output,'Sheet','Medial T2SFFE','Range','A2')
    writematrix(coordm_tf,output,'Sheet','Medial T2SFFE','Range','B2')
    writematrix(mean(avg_tfmed),output,'Sheet','Medial T2SFFE','Range','C2')
    writematrix(std(avg_tfmed),output,'Sheet','Medial T2SFFE','Range','D2')
    writematrix(max(avg_tfmed),output,'Sheet','Medial T2SFFE','Range','E2')
    writematrix(min(avg_tfmed),output,'Sheet','Medial T2SFFE','Range','F2')
    writecell(col_header,output,'Sheet','Medial T2SFFE','Range','A1')


    return