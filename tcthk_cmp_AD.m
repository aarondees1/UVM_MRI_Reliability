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

avg_lats = zeros(3,3,ns);
avg_meds = zeros(3,3,ns);
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
    %
    % Lateral Compartment
    %
    ilat(:,3) = xgl>reg_lat2_all;                   % Anterior region
    ilat(:,2) = xgl<reg_lat2_all&xgl>reg_lat1_all;          % Central region
    ilat(:,1) = xgl<reg_lat1_all;                   % Posterior region
    %
    % Medial Compartment
    %
    imed(:,3) = xgm>reg_med2_all;                   % Anterior region
    imed(:,2) = xgm<reg_med2_all&xgm>reg_med1_all;         % Central region
    imed(:,1) = xgm<reg_med1_all;                   % Posterior region

    avg_lats(1,1,i) = mean(cthklf(ilat(:,1)),"omitnan");
    avg_lats(2,1,i) = mean(cthklr(ilat(:,1)),"omitnan");
    avg_lats(3,1,i) = mean(cthklt(ilat(:,1)),"omitnan");

    avg_lats(1,2,i) = mean(cthklf(ilat(:,2)),"omitnan");
    avg_lats(2,2,i) = mean(cthklr(ilat(:,2)),"omitnan");
    avg_lats(3,2,i) = mean(cthklt(ilat(:,2)),"omitnan");

    avg_lats(1,3,i) = mean(cthklf(ilat(:,3)),"omitnan");
    avg_lats(2,3,i) = mean(cthklr(ilat(:,3)),"omitnan");
    avg_lats(3,3,i) = mean(cthklt(ilat(:,3)),"omitnan");

    avg_meds(1,1,i) = mean(cthkmf(imed(:,1)),"omitnan");
    avg_meds(2,1,i) = mean(cthkmr(imed(:,1)),"omitnan");
    avg_meds(3,1,i) = mean(cthkmt(imed(:,1)),"omitnan");

    avg_meds(1,2,i) = mean(cthkmf(imed(:,2)),"omitnan");
    avg_meds(2,2,i) = mean(cthkmr(imed(:,2)),"omitnan");
    avg_meds(3,2,i) = mean(cthkmt(imed(:,2)),"omitnan");

    avg_meds(1,3,i) = mean(cthkmf(imed(:,3)),"omitnan");
    avg_meds(2,3,i) = mean(cthkmr(imed(:,3)),"omitnan");
    avg_meds(3,3,i) = mean(cthkmt(imed(:,3)),"omitnan");


    rflat1 = avg_lats(2,1,i)-avg_lats(1,1,i);
    rflat2 = avg_lats(2,2,i)-avg_lats(1,2,i);
    rflat3 = avg_lats(2,3,i)-avg_lats(1,3,i);
    tflat1 = avg_lats(3,1,i)-avg_lats(1,1,i);
    tflat2 = avg_lats(3,2,i)-avg_lats(1,2,i);
    tflat3 = avg_lats(3,3,i)-avg_lats(1,3,i);
    rfmed1 = avg_meds(2,1,i)-avg_meds(1,1,i);
    rfmed2 = avg_meds(2,2,i)-avg_meds(1,2,i);
    rfmed3 = avg_meds(2,3,i)-avg_meds(1,3,i);
    tfmed1 = avg_meds(3,1,i)-avg_meds(1,1,i);
    tfmed2 = avg_meds(3,2,i)-avg_meds(1,2,i);
    tfmed3 = avg_meds(3,3,i)-avg_meds(1,3,i);
    
    %%
    %
    % Plots of Cartilage Thicknesses and Thickness Differences
    %
    f = figure('WindowState','fullscreen');
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
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmr;
    %cthk(~idvm_rf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
    axis equal;
    view(-90,90);
    title('T1\rho Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf1, 'Ticks',0:6);
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    text(-15,20,num2str(avg_lats(2,1,i)));
    text(-4,20,num2str(avg_lats(2,2,i)));
    text(7,20,num2str(avg_lats(2,3,i)));
    text(-13,-20,num2str(avg_meds(2,1,i)));
    text(0,-20,num2str(avg_meds(2,2,i)));
    text(13,-20,num2str(avg_meds(2,3,i)));    


    %
   hf2=nexttile(3,[1 2]);
    cthk = cthklf;
    %cthk(~idvl_rf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmf;
    %cthk(~idvm_rf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
    axis equal;
    view(-90,90);
    title('T1FFE Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf2, 'Ticks',0:6);
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    text(-15,20,num2str(avg_lats(1,1,i)));
    text(-4,20,num2str(avg_lats(1,2,i)));
    text(7,20,num2str(avg_lats(1,3,i)));
    text(-13,-20,num2str(avg_meds(1,1,i)));
    text(0,-20,num2str(avg_meds(1,2,i)));
    text(13,-20,num2str(avg_meds(1,3,i))); 


    hf3=nexttile(5,[1 2]);
    cthk = cthklt;
    %cthk(~idvl_tf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmt;
    %cthk(~idvm_tf) = NaN;
    patch(xgm(quadm'),ygm(quadm'),cthk(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Med of ',fstr, ' = ', num2str(thk_max), newline]);
    axis equal;
    view(-90,90);
    title('T2S Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf3, 'Ticks',0:6);
    clim([0 6]);
    xlim([-30 30]);
    ylim([-40 40]);
    scatter(xgl,ygl,4,".",'CData', grayColor);
    scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    text(-15,20,num2str(avg_lats(3,1,i)));
    text(-4,20,num2str(avg_lats(3,2,i)));
    text(7,20,num2str(avg_lats(3,3,i)));
    text(-13,-20,num2str(avg_meds(3,1,i)));
    text(0,-20,num2str(avg_meds(3,2,i)));
    text(13,-20,num2str(avg_meds(3,3,i)));    
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
    colorbar(hf4,'Ticks',-2.5:.5:2.5);
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
    text(-15,20,num2str(rflat1));
    text(-4,20,num2str(rflat2));
    text(7,20,num2str(rflat3));
    text(-13,-20,num2str(rfmed1));
    text(0,-20,num2str(rfmed2));
    text(13,-20,num2str(rfmed3));
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
    colorbar(hf5,'Ticks',-2.5:.5:2.5);
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
    text(-15,20,num2str(tflat1));
    text(-4,20,num2str(tflat2));
    text(7,20,num2str(tflat3));
    text(-13,-20,num2str(tfmed1));
    text(0,-20,num2str(tfmed2));
    text(13,-20,num2str(tfmed3));
    %
    pic_nam=fullfile(rdir, 'Thickness Differences.pdf');
    exportgraphics(f, pic_nam, 'Append', true);


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
    f = figure('WindowState','fullscreen');
    tiledlayout(2,2,"Padding","loose");
    sgtitle({[fstr(1:5) ' - Tibia'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;
    %
    nexttile;
    plot(dthklrf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanlrf dmeanlrf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf+3*dstdlrf dmeanlrf+3*dstdlrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf-3*dstdlrf dmeanlrf-3*dstdlrf],'r--','LineWidth',1);
    title('Rho/FFE Lateral Differences', ...
        'FontSize',16,'FontWeight','bold');
    

    nexttile;
    plot(dthkltf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanltf dmeanltf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanltf+3*dstdltf dmeanltf+3*dstdltf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanltf-3*dstdltf dmeanltf-3*dstdltf],'r--','LineWidth',1);
    title('T2S/FFE Lateral Differences', ...
        'FontSize',16,'FontWeight','bold');

    %
    nexttile;
    plot(dthkmrf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanmrf dmeanmrf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf+3*dstdmrf dmeanmrf+3*dstdmrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf-3*dstdmrf dmeanmrf-3*dstdmrf],'r--','LineWidth',1);
    title('Rho/FFE Medial Differences', ...
        'FontSize',16,'FontWeight','bold');


    %
    nexttile;
    plot(dthkmtf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanmtf dmeanmtf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf+3*dstdmtf dmeanmtf+3*dstdmtf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf-3*dstdmtf dmeanmtf-3*dstdmtf],'r--','LineWidth',1);
    title('T2S/FFE Medial Differences', ...
        'FontSize',16,'FontWeight','bold');


    % Write Differences to CSV Spreadsheet
    output = fullfile(rdir,'Thickness Differences.xlsx');
    col_header1 = {'LATERAL RHO - FFE'};
    col_header2 = {'MEDIAL RHO - FFE'};
    col_header3 = {'LATERAL T2S - FFE'};
    col_header4 = {'MEDIAL T2S - FFE'};
    col_header5 = {'Thickness Diff (mm)','Grid Coordinates','Mean'};
    col_header6 = {'STD'};
    col_header7 = {'Max'};
    col_header8 = {'Min'};
    col_header9 = {'Posterior Average'};
    col_header10 = {'Central Average'};
    col_header11 = {'Anterior Average'};

    writecell(col_header1,output,'Sheet',fstr(1:5),'Range','A1')
    writecell(col_header2,output,'Sheet',fstr(1:5),'Range','E1')
    writecell(col_header3,output,'Sheet',fstr(1:5),'Range','I1')
    writecell(col_header4,output,'Sheet',fstr(1:5),'Range','M1')
    writecell(col_header5,output,'Sheet',fstr(1:5),'Range','A2')
    writecell(col_header5,output,'Sheet',fstr(1:5),'Range','E2')
    writecell(col_header5,output,'Sheet',fstr(1:5),'Range','I2')
    writecell(col_header5,output,'Sheet',fstr(1:5),'Range','M2')
    writecell(col_header6,output,'Sheet',fstr(1:5),'Range','C5')
    writecell(col_header6,output,'Sheet',fstr(1:5),'Range','G5')
    writecell(col_header6,output,'Sheet',fstr(1:5),'Range','K5')
    writecell(col_header6,output,'Sheet',fstr(1:5),'Range','O5')
    writecell(col_header7,output,'Sheet',fstr(1:5),'Range','C8')
    writecell(col_header7,output,'Sheet',fstr(1:5),'Range','G8')
    writecell(col_header7,output,'Sheet',fstr(1:5),'Range','K8')
    writecell(col_header7,output,'Sheet',fstr(1:5),'Range','O8')
    writecell(col_header8,output,'Sheet',fstr(1:5),'Range','C11')
    writecell(col_header8,output,'Sheet',fstr(1:5),'Range','G11')
    writecell(col_header8,output,'Sheet',fstr(1:5),'Range','K11')
    writecell(col_header8,output,'Sheet',fstr(1:5),'Range','O11')
    writecell(col_header9,output,'Sheet',fstr(1:5),'Range','C14')
    writecell(col_header9,output,'Sheet',fstr(1:5),'Range','G14')
    writecell(col_header9,output,'Sheet',fstr(1:5),'Range','K14')
    writecell(col_header9,output,'Sheet',fstr(1:5),'Range','O14')
    writecell(col_header10,output,'Sheet',fstr(1:5),'Range','C17')
    writecell(col_header10,output,'Sheet',fstr(1:5),'Range','G17')
    writecell(col_header10,output,'Sheet',fstr(1:5),'Range','K17')
    writecell(col_header10,output,'Sheet',fstr(1:5),'Range','O17')
    writecell(col_header11,output,'Sheet',fstr(1:5),'Range','C20')
    writecell(col_header11,output,'Sheet',fstr(1:5),'Range','G20')
    writecell(col_header11,output,'Sheet',fstr(1:5),'Range','K20')
    writecell(col_header11,output,'Sheet',fstr(1:5),'Range','O20')
    
    writematrix(dthklrf,output,'Sheet',fstr(1:5),'Range','A3')
    writematrix(coordl_rf,output,'Sheet',fstr(1:5),'Range','B3')
    writematrix(dmeanlrf,output,'Sheet',fstr(1:5),'Range','C3')
    writematrix(dstdlrf,output,'Sheet',fstr(1:5),'Range','C6')
    writematrix(dmaxlrf,output,'Sheet',fstr(1:5),'Range','C9')
    writematrix(dminlrf,output,'Sheet',fstr(1:5),'Range','C12')
    writematrix(rflat1,output,'Sheet',fstr(1:5),'Range','C15')
    writematrix(rflat2,output,'Sheet',fstr(1:5),'Range','C18')
    writematrix(rflat3,output,'Sheet',fstr(1:5),'Range','C21')

    writematrix(dthkltf,output,'Sheet',fstr(1:5),'Range','E3')
    writematrix(coordl_tf,output,'Sheet',fstr(1:5),'Range','F3')
    writematrix(dmeanltf,output,'Sheet',fstr(1:5),'Range','G3')
    writematrix(dstdltf,output,'Sheet',fstr(1:5),'Range','G6')
    writematrix(dmaxltf,output,'Sheet',fstr(1:5),'Range','G9')
    writematrix(dminltf,output,'Sheet',fstr(1:5),'Range','G12')
    writematrix(tflat1,output,'Sheet',fstr(1:5),'Range','G15')
    writematrix(tflat2,output,'Sheet',fstr(1:5),'Range','G18')
    writematrix(tflat3,output,'Sheet',fstr(1:5),'Range','G21')

    writematrix(dthkmrf,output,'Sheet',fstr(1:5),'Range','I3')
    writematrix(coordm_rf,output,'Sheet',fstr(1:5),'Range','J3')
    writematrix(dmeanmrf,output,'Sheet',fstr(1:5),'Range','K3')
    writematrix(dstdmrf,output,'Sheet',fstr(1:5),'Range','K6')
    writematrix(dmaxmrf,output,'Sheet',fstr(1:5),'Range','K9')
    writematrix(dminmrf,output,'Sheet',fstr(1:5),'Range','K12')
    writematrix(rfmed1,output,'Sheet',fstr(1:5),'Range','K15')
    writematrix(rfmed2,output,'Sheet',fstr(1:5),'Range','K18')
    writematrix(rfmed3,output,'Sheet',fstr(1:5),'Range','K21')

    writematrix(dthkmtf,output,'Sheet',fstr(1:5),'Range','M3')
    writematrix(coordm_tf,output,'Sheet',fstr(1:5),'Range','N3')
    writematrix(dmeanmtf,output,'Sheet',fstr(1:5),'Range','O3')
    writematrix(dstdmtf,output,'Sheet',fstr(1:5),'Range','O6')
    writematrix(dmaxmtf,output,'Sheet',fstr(1:5),'Range','O9')
    writematrix(dminmtf,output,'Sheet',fstr(1:5),'Range','O12')
    writematrix(tfmed1,output,'Sheet',fstr(1:5),'Range','O15')
    writematrix(tfmed2,output,'Sheet',fstr(1:5),'Range','O18')
    writematrix(tfmed3,output,'Sheet',fstr(1:5),'Range','O21')    

    %
    exportgraphics(f, pic_nam, 'Append', true);
    %
end

%%
n_rflat = sum(~isnan(cthkld_rf),3);
idxn = find(n_rflat<3);
n_rflat(idxn) = NaN;
avg_rflat = mean(cthkld_rf,3,"omitnan");
avg_rflat(idxn) = NaN;
std_rflat = std(cthkld_rf,0,3,"omitnan");
std_rflat(idxn) = NaN;

n_tflat = sum(~isnan(cthkld_tf),3);
idxn = find(n_tflat<3);
n_tflat(idxn) = NaN;
avg_tflat = mean(cthkld_tf,3,"omitnan");
avg_tflat(idxn) = NaN;
std_tflat = std(cthkld_tf,0,3,"omitnan");
std_tflat(idxn) = NaN;


n_rfmed = sum(~isnan(cthkmd_rf),3);
idxn = find(n_rfmed<3);
n_rfmed(idxn) = NaN;
avg_rfmed = mean(cthkmd_rf,3,"omitnan");
avg_rfmed(idxn) = NaN;
std_rfmed = std(cthkmd_rf,0,3,"omitnan");
std_rfmed(idxn) = NaN;

n_tfmed = sum(~isnan(cthkmd_tf),3);
idxn = find(n_tfmed<3);
n_tfmed(n_tfmed<3) = NaN;
avg_tfmed = mean(cthkmd_tf,3,"omitnan");
avg_tfmed(idxn) = NaN;
std_tfmed = std(cthkmd_tf,0,3,"omitnan");
std_tfmed(idxn) = NaN;


rflat1 = mean(avg_rflat(ilat(:,1)),"omitnan");
rflat2 = mean(avg_rflat(ilat(:,2)),"omitnan");
rflat3 = mean(avg_rflat(ilat(:,3)),"omitnan");
tflat1 = mean(avg_tflat(ilat(:,1)),"omitnan");
tflat2 = mean(avg_tflat(ilat(:,2)),"omitnan");
tflat3 = mean(avg_tflat(ilat(:,3)),"omitnan");
rfmed1 = mean(avg_rfmed(imed(:,1)),"omitnan");
rfmed2 = mean(avg_rfmed(imed(:,2)),"omitnan");
rfmed3 = mean(avg_rfmed(imed(:,3)),"omitnan");
tfmed1 = mean(avg_tfmed(imed(:,1)),"omitnan");
tfmed2 = mean(avg_tfmed(imed(:,2)),"omitnan");
tfmed3 = mean(avg_tfmed(imed(:,3)),"omitnan");
% 0 - 10 Colormap for Number of Subjects
%cmap=[0 0 0.6; 0 0 1; 0 0.4 1; 0 0.8 1; 0.2 1 0.8; 0.6 1 0.4;
%1 1 0; 1 0.6 0; 1 0.4 0; 1 0 0; 0.6 0 0];

% 3 - 10 Colormap for Number of Subjects
cmap=[0 0.8 1; 0.2 1 0.8; 0.6 1 0.4; 1 1 0; 1 0.6 0; 1 0.4 0;
    1 0 0; 0.6 0 0];
dim1 = [0.2 0.2 0.3 0.3];

f = figure('WindowState','fullscreen');
t = tiledlayout(3,2);
title(t,'Tibial Cartilage Thickness Differences','FontSize',20,'FontWeight','bold');
orient tall;

hf1 = nexttile;
patch(xgl(quadl'),ygl(quadl'),avg_rflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),avg_rfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('Rho - T1FFE Averages', ...
    'FontSize',14,'FontWeight','bold');
colormap(hf1,parula(10));
hf1.CLim = [-2.5 2.5];
colorbar(hf1, 'Ticks', -2.5:.5:2.5);
xlim(hf1,[-30 30]);
ylim(hf1,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
text(-15,20,num2str(rflat1));
text(-4,20,num2str(rflat2));
text(7,20,num2str(rflat3));
text(-13,-20,num2str(rfmed1));
text(0,-20,num2str(rfmed2));
text(13,-20,num2str(rfmed3));

hf2 = nexttile;
patch(xgl(quadl'),ygl(quadl'),avg_tflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),avg_tfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('T2S - T1FFE Averages', ...
    'FontSize',14,'FontWeight','bold');
colormap(hf2,parula(10));
hf2.CLim = [-2.5 2.5];
colorbar(hf2, 'Ticks', -2.5:.5:2.5);
xlim(hf2,[-30 30]);
ylim(hf2,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
text(-15,20,num2str(tflat1));
text(-4,20,num2str(tflat2));
text(7,20,num2str(tflat3));
text(-13,-20,num2str(tfmed1));
text(0,-20,num2str(tfmed2));
text(13,-20,num2str(tfmed3));

hf3 = nexttile;
patch(xgl(quadl'),ygl(quadl'),std_rflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),std_rfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('Rho - T1FFE','Standard Deviatons', ...
    'FontSize',14,'FontWeight','bold');
hf3.CLim = [0 1.5];
colormap(hf3,jet(6));
colorbar(hf3,'Ticks', 0:.25:1.5);
xlim(hf3,[-30 30]);
ylim(hf3,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");

hf4 = nexttile;
patch(xgl(quadl'),ygl(quadl'),std_tflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),std_tfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('T2S - T1FFE','Standard Deviatons', ...
    'FontSize',14,'FontWeight','bold');
hf4.CLim = [0 1.5];
colormap(hf4,jet(6));
colorbar(hf4,'Ticks', 0:.25:1.5);
xlim(hf4,[-30 30]);
ylim(hf4,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");

hf5 = nexttile;
patch(xgl(quadl'),ygl(quadl'),n_rflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),n_rfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('Rho - T1FFE','Subject Numbers', ...
    'FontSize',14,'FontWeight','bold');
%hf5.CLim=[0 11];
hf5.CLim=[2 10];
colormap(hf5,cmap);
%colorbar(hf5,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
colorbar(hf5,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
xlim(hf5,[-30 30]);
ylim(hf5,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);

hf6 = nexttile;
patch(xgl(quadl'),ygl(quadl'),n_tflat(quadl'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
patch(xgm(quadm'),ygm(quadm'),n_tfmed(quadm'),'FaceColor','interp', ...
    'EdgeColor','interp');
view(-90,90);
axis equal;
title('T2S - T1FFE','Subject Numbers', ...
    'FontSize',14,'FontWeight','bold');
%hf6.CLim=[0 11];
hf6.CLim=[2 10];
colormap(hf6,cmap);
%colorbar(hf6,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
colorbar(hf6,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
xlim(hf6,[-30 30]);
ylim(hf6,[-40 40]);
scatter(xgl,ygl,4,".",'CData', grayColor);
scatter(xgm,ygm,4,".",'CData', grayColor);

exportgraphics(f,pic_nam, 'Append', true);
%%


%Take out NaNs

coordl_rf = find(~isnan(avg_rflat));
coordl_tf = find(~isnan(avg_tflat));
coordm_rf = find(~isnan(avg_rfmed));
coordm_tf = find(~isnan(avg_tfmed));

avg_rflat = avg_rflat((~isnan(avg_rflat)));
avg_tflat = avg_tflat((~isnan(avg_tflat)));
avg_rfmed = avg_rfmed((~isnan(avg_rfmed)));
avg_tfmed = avg_tfmed((~isnan(avg_tfmed)));


output = fullfile(rdir,'Thickness Differences.xlsx');

    writecell(col_header1,output,'Sheet','Averages','Range','A1')
    writecell(col_header2,output,'Sheet','Averages','Range','E1')
    writecell(col_header3,output,'Sheet','Averages','Range','I1')
    writecell(col_header4,output,'Sheet','Averages','Range','M1')
    writecell(col_header5,output,'Sheet','Averages','Range','A2')
    writecell(col_header5,output,'Sheet','Averages','Range','E2')
    writecell(col_header5,output,'Sheet','Averages','Range','I2')
    writecell(col_header5,output,'Sheet','Averages','Range','M2')
    writecell(col_header6,output,'Sheet','Averages','Range','C5')
    writecell(col_header6,output,'Sheet','Averages','Range','G5')
    writecell(col_header6,output,'Sheet','Averages','Range','K5')
    writecell(col_header6,output,'Sheet','Averages','Range','O5')
    writecell(col_header7,output,'Sheet','Averages','Range','C8')
    writecell(col_header7,output,'Sheet','Averages','Range','G8')
    writecell(col_header7,output,'Sheet','Averages','Range','K8')
    writecell(col_header7,output,'Sheet','Averages','Range','O8')
    writecell(col_header8,output,'Sheet','Averages','Range','C11')
    writecell(col_header8,output,'Sheet','Averages','Range','G11')
    writecell(col_header8,output,'Sheet','Averages','Range','K11')
    writecell(col_header8,output,'Sheet','Averages','Range','O11')
    writecell(col_header9,output,'Sheet','Averages','Range','C14')
    writecell(col_header9,output,'Sheet','Averages','Range','G14')
    writecell(col_header9,output,'Sheet','Averages','Range','K14')
    writecell(col_header9,output,'Sheet','Averages','Range','O14')
    writecell(col_header10,output,'Sheet','Averages','Range','C17')
    writecell(col_header10,output,'Sheet','Averages','Range','G17')
    writecell(col_header10,output,'Sheet','Averages','Range','K17')
    writecell(col_header10,output,'Sheet','Averages','Range','O17')
    writecell(col_header11,output,'Sheet','Averages','Range','C20')
    writecell(col_header11,output,'Sheet','Averages','Range','G20')
    writecell(col_header11,output,'Sheet','Averages','Range','K20')
    writecell(col_header11,output,'Sheet','Averages','Range','O20')

writematrix(avg_rflat,output,'Sheet','Averages','Range','A3')
writematrix(coordl_rf,output,'Sheet','Averages','Range','B3')
writematrix(mean(avg_rflat,"omitnan"),output,'Sheet','Averages','Range','C3')
writematrix(std(avg_rflat,"omitnan"),output,'Sheet','Averages','Range','C6')
writematrix(max(avg_rflat),output,'Sheet','Averages','Range','C9')
writematrix(min(avg_rflat),output,'Sheet','Averages','Range','C12')
writematrix(rflat1,output,'Sheet','Averages','Range','C15')
writematrix(rflat2,output,'Sheet','Averages','Range','C18')
writematrix(rflat3,output,'Sheet','Averages','Range','C21')


writematrix(avg_tflat,output,'Sheet','Averages','Range','E3')
writematrix(coordl_tf,output,'Sheet','Averages','Range','F3')
writematrix(mean(avg_tflat),output,'Sheet','Averages','Range','G3')
writematrix(std(avg_tflat),output,'Sheet','Averages','Range','G6')
writematrix(max(avg_tflat),output,'Sheet','Averages','Range','G9')
writematrix(min(avg_tflat),output,'Sheet','Averages','Range','G12')
writematrix(tflat1,output,'Sheet','Averages','Range','G15')
writematrix(tflat2,output,'Sheet','Averages','Range','G18')
writematrix(tflat3,output,'Sheet','Averages','Range','G21')


writematrix(avg_rfmed,output,'Sheet','Averages','Range','I3')
writematrix(coordm_rf,output,'Sheet','Averages','Range','J3')
writematrix(mean(avg_rfmed),output,'Sheet','Averages','Range','K3')
writematrix(std(avg_rfmed),output,'Sheet','Averages','Range','K6')
writematrix(max(avg_rfmed),output,'Sheet','Averages','Range','K9')
writematrix(min(avg_rfmed),output,'Sheet','Averages','Range','K12')
writematrix(rfmed1,output,'Sheet','Averages','Range','K15')
writematrix(rfmed2,output,'Sheet','Averages','Range','K18')
writematrix(rfmed3,output,'Sheet','Averages','Range','K21')


writematrix(avg_tfmed,output,'Sheet','Averages','Range','M3')
writematrix(coordm_tf,output,'Sheet','Averages','Range','N3')
writematrix(mean(avg_tfmed),output,'Sheet','Averages','Range','O3')
writematrix(std(avg_tfmed),output,'Sheet','Averages','Range','O6')
writematrix(max(avg_tfmed),output,'Sheet','Averages','Range','O9')
writematrix(min(avg_tfmed),output,'Sheet','Averages','Range','O12')
writematrix(tfmed1,output,'Sheet','Averages','Range','O15')
writematrix(tfmed2,output,'Sheet','Averages','Range','O18')
writematrix(tfmed3,output,'Sheet','Averages','Range','O21')

return