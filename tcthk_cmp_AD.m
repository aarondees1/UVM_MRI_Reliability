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
alphabet = [char(65:90)];

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
cthkmf2 = [];
cthkmr2 = [];
cthkmt2 = [];
cthklf2 = [];
cthklr2 = [];
cthklt2 = [];
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

    cthklf2 = cat(3,cthklf2,cthklf);
    cthklr2 = cat(3,cthklr2,cthklr);
    cthklt2 = cat(3,cthklt2,cthklt);
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

    cthkmf2 = cat(3,cthkmf2,cthkmf);
    cthkmr2 = cat(3,cthkmr2,cthkmr);
    cthkmt2 = cat(3,cthkmt2,cthkmt);
    %
    % Calculate Differences
    %
    cthkld_rf(:,:,i) = cthklr-cthklf;
    cthkld_tf(:,:,i) = cthklt-cthklf;
    cthkmd_rf(:,:,i) = cthkmr-cthkmf;
    cthkmd_tf(:,:,i) = cthkmt-cthkmf;

    % idvl_rf = ~isnan(cthkld_rf(:,:,i));  % Valid differences
    % idvm_rf = ~isnan(cthkmd_rf(:,:,i));  % Valid differences
    % idvl_tf = ~isnan(cthkld_tf(:,:,i));  % Valid differences
    % idvm_tf = ~isnan(cthkmd_tf(:,:,i));  % Valid differences
    %
    % idvl = idvl_rf & idvl_tf;
    % idvm = idvm_rf & idvm_tf;
    %
    % cthklf(~idvl) = NaN;
    % cthklr(~idvl) = NaN;
    % cthklt(~idvl) = NaN;
    %
    % cthkmf(~idvm) = NaN;
    % cthkmr(~idvm) = NaN;
    % cthkmt(~idvm) = NaN;

end

%%
%
idx_cthklf = ~isnan(cthklf2);
idx_cthklr = ~isnan(cthklr2);
idx_cthklt = ~isnan(cthklt2);
idx_cthkmf = ~isnan(cthkmf2);
idx_cthkmr = ~isnan(cthkmr2);
idx_cthkmt = ~isnan(cthkmt2);

lat_cols = sum(max(idx_cthklf),3) + sum(max(idx_cthklr),3) + sum(max(idx_cthklt),3);
lat_rows = sum(max(idx_cthklf,[],2),3) + sum(max(idx_cthklr,[],2),3) + sum(max(idx_cthklt,[],2),3);
lat_cols = logical(lat_cols);
lat_rows = logical(lat_rows);

med_cols = sum(max(idx_cthkmf),3) + sum(max(idx_cthkmr),3) + sum(max(idx_cthkmt),3);
med_rows = sum(max(idx_cthkmf,[],2),3) + sum(max(idx_cthkmr,[],2),3) + sum(max(idx_cthkmt,[],2),3);
med_cols = logical(med_cols);
med_rows = logical(med_rows);

lat_col1 = find(lat_cols,1,'first');
lat_col2 = find(lat_cols,1,'last');
med_col1 = find(med_cols,1,'first');
med_col2 = find(med_cols,1,'last');

lat_row1 = find(lat_rows,1,'first');
lat_row2 = find(lat_rows,1,'last');
med_row1 = find(med_rows,1,'first');
med_row2 = find(med_rows,1,'last');

cthklf2 = cthklf2(lat_row1:lat_row2,lat_col1:lat_col2,:);
cthklr2 = cthklr2(lat_row1:lat_row2,lat_col1:lat_col2,:);
cthklt2 = cthklt2(lat_row1:lat_row2,lat_col1:lat_col2,:);
cthkmf2 = cthkmf2(med_row1:med_row2,med_col1:med_col2,:);
cthkmr2 = cthkmr2(med_row1:med_row2,med_col1:med_col2,:);
cthkmt2 = cthkmt2(med_row1:med_row2,med_col1:med_col2,:);

cthkld_rf = cthkld_rf(lat_row1:lat_row2,lat_col1:lat_col2,:);
cthkld_tf = cthkld_tf(lat_row1:lat_row2,lat_col1:lat_col2,:);
cthkmd_rf = cthkmd_rf(med_row1:med_row2,med_col1:med_col2,:);
cthkmd_tf = cthkmd_tf(med_row1:med_row2,med_col1:med_col2,:);

xgl = xgl(lat_row1:lat_row2,lat_col1:lat_col2);
ygl = ygl(lat_row1:lat_row2,lat_col1:lat_col2);
xgm = xgm(med_row1:med_row2,med_col1:med_col2);
ygm = ygm(med_row1:med_row2,med_col1:med_col2);

[nlr,nlc] = size(xgl);
xgl = xgl(:);
ygl = ygl(:);
quadl = quadconn(nlr,nlc);

[nmr,nmc] = size(xgm);
xgm = xgm(:);
ygm = ygm(:);
quadm = quadconn(nmr,nmc);

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

roi(1) = "tib_lat_pos.xlsx";
roi(2) = "tib_med_pos.xlsx";
roi(3) = "tib_lat_ctr.xlsx";
roi(4) = "tib_med_ctr.xlsx";
roi(5) = "tib_lat_ant.xlsx";
roi(6) = "tib_med_ant.xlsx";


coordl = 1:size(xgl);
coordl = coordl(:);
coordm = 1:size(xgm);
coordm = coordm(:);
% rflat1 = avg_lats(2,1,i)-avg_lats(1,1,i);
% rflat2 = avg_lats(2,2,i)-avg_lats(1,2,i);
% rflat3 = avg_lats(2,3,i)-avg_lats(1,3,i);
% tflat1 = avg_lats(3,1,i)-avg_lats(1,1,i);
% tflat2 = avg_lats(3,2,i)-avg_lats(1,2,i);
% tflat3 = avg_lats(3,3,i)-avg_lats(1,3,i);
% rfmed1 = avg_meds(2,1,i)-avg_meds(1,1,i);
% rfmed2 = avg_meds(2,2,i)-avg_meds(1,2,i);
% rfmed3 = avg_meds(2,3,i)-avg_meds(1,3,i);
% tfmed1 = avg_meds(3,1,i)-avg_meds(1,1,i);
% tfmed2 = avg_meds(3,2,i)-avg_meds(1,2,i);
% tfmed3 = avg_meds(3,3,i)-avg_meds(1,3,i);

%%
ii=0;
for i = 1:ns
   

    % Turning the thicknesses into vectors
    cthklf_v = cthklf2(:,:,i);
    cthklr_v = cthklr2(:,:,i);
    cthklt_v = cthklt2(:,:,i);
    cthkmf_v = cthkmf2(:,:,i);
    cthkmr_v = cthkmr2(:,:,i);
    cthkmt_v = cthkmt2(:,:,i);
    cthklf_v = cthklf_v(:);
    cthklr_v = cthklr_v(:);
    cthklt_v = cthklt_v(:);
    cthkmf_v = cthkmf_v(:);
    cthkmr_v = cthkmr_v(:);
    cthkmt_v = cthkmt_v(:);

    avg_lats(1,1) = mean(cthklf_v(ilat(:,1)),"omitnan");
    avg_lats(2,1) = mean(cthklr_v(ilat(:,1)),"omitnan");
    avg_lats(3,1) = mean(cthklt_v(ilat(:,1)),"omitnan");

    avg_lats(1,2) = mean(cthklf_v(ilat(:,2)),"omitnan");
    avg_lats(2,2) = mean(cthklr_v(ilat(:,2)),"omitnan");
    avg_lats(3,2) = mean(cthklt_v(ilat(:,2)),"omitnan");

    avg_lats(1,3) = mean(cthklf_v(ilat(:,3)),"omitnan");
    avg_lats(2,3) = mean(cthklr_v(ilat(:,3)),"omitnan");
    avg_lats(3,3) = mean(cthklt_v(ilat(:,3)),"omitnan");

    avg_meds(1,1) = mean(cthkmf_v(imed(:,1)),"omitnan");
    avg_meds(2,1) = mean(cthkmr_v(imed(:,1)),"omitnan");
    avg_meds(3,1) = mean(cthkmt_v(imed(:,1)),"omitnan");

    avg_meds(1,2) = mean(cthkmf_v(imed(:,2)),"omitnan");
    avg_meds(2,2) = mean(cthkmr_v(imed(:,2)),"omitnan");
    avg_meds(3,2) = mean(cthkmt_v(imed(:,2)),"omitnan");

    avg_meds(1,3) = mean(cthkmf_v(imed(:,3)),"omitnan");
    avg_meds(2,3) = mean(cthkmr_v(imed(:,3)),"omitnan");
    avg_meds(3,3) = mean(cthkmt_v(imed(:,3)),"omitnan");

    fstr=sd(i).FFE.tibia.bfnam;
    %
    % Plots of Cartilage Thicknesses and Thickness Differences
    %
    f = figure;
    t = tiledlayout(2,6);
    sgtitle({[fstr(1:5) ' - Tibia'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;
    %

    hf1=nexttile(1,[1 2]);
    cthk = cthklr_v;
    %cthk(~idvl_rf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmr_v;
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
    % scatter(xgl,ygl,4,".",'CData', grayColor);
    % scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    % text(-15,20,num2str(avg_lats(2,1)));
    % text(-4,20,num2str(avg_lats(2,2)));
    % text(7,20,num2str(avg_lats(2,3)));
    % text(-13,-20,num2str(avg_meds(2,1)));
    % text(0,-20,num2str(avg_meds(2,2)));
    % text(13,-20,num2str(avg_meds(2,3)));


    %
    hf2=nexttile(3,[1 2]);
    cthk = cthklf_v;
    %cthk(~idvl_rf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmf_v;
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
    % scatter(xgl,ygl,4,".",'CData', grayColor);
    % scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    % text(-15,20,num2str(avg_lats(1,1)));
    % text(-4,20,num2str(avg_lats(1,2)));
    % text(7,20,num2str(avg_lats(1,3)));
    % text(-13,-20,num2str(avg_meds(1,1)));
    % text(0,-20,num2str(avg_meds(1,2)));
    % text(13,-20,num2str(avg_meds(1,3)));


    hf3=nexttile(5,[1 2]);
    cthk = cthklt_v;
    %cthk(~idvl_tf) = NaN;
    patch(xgl(quadl'),ygl(quadl'),cthk(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    cthk = cthkmt_v;
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
    % scatter(xgl,ygl,4,".",'CData', grayColor);
    % scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    % text(-15,20,num2str(avg_lats(3,1)));
    % text(-4,20,num2str(avg_lats(3,2)));
    % text(7,20,num2str(avg_lats(3,3)));
    % text(-13,-20,num2str(avg_meds(3,1)));
    % text(0,-20,num2str(avg_meds(3,2)));
    % text(13,-20,num2str(avg_meds(3,3)));
    %

    % colormap(hf1, 'jet(12)');
    % colormap(hf2, 'jet(12)');
    % colormap(hf3, 'jet(12)');
    colormap(hf1, 'jet');
    colormap(hf2, 'jet');
    colormap(hf3, 'jet');
    %

    hf4=nexttile(7,[1 3]);
    dthklrf = cthkld_rf(:,:,i);
    %dthklrf(~idvl) = NaN;
    patch(xgl(quadl'),ygl(quadl'),dthklrf(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    dthkmrf = cthkmd_rf(:,:,i);
    %dthkmrf(~idvm) = NaN;
    patch(xgm(quadm'),ygm(quadm'),dthkmrf(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    axis equal;
    view(-90,90);
    title('T1\rho - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    % colorbar(hf4,'Ticks',-2.5:.5:2.5);
    colorbar;
    clim([-1.5 1.5]);
    xlim([-30 30]);
    ylim([-40 40]);
    % scatter(xgl,ygl,4,".",'CData', grayColor);
    % scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    rflat1 = mean(dthklrf(ilat(:,1)),"omitnan");
    rflat2 = mean(dthklrf(ilat(:,2)),"omitnan");
    rflat3 = mean(dthklrf(ilat(:,3)),"omitnan");
    rfmed1 = mean(dthkmrf(imed(:,1)),"omitnan");
    rfmed2 = mean(dthkmrf(imed(:,2)),"omitnan");
    rfmed3 = mean(dthkmrf(imed(:,3)),"omitnan");
    % text(-15,20,num2str(rflat1));
    % text(-4,20,num2str(rflat2));
    % text(7,20,num2str(rflat3));
    % text(-13,-20,num2str(rfmed1));
    % text(0,-20,num2str(rfmed2));
    % text(13,-20,num2str(rfmed3));


    hf5=nexttile(10,[1 3]);
    dthkltf = cthkld_tf(:,:,i);
    %dthkltf(~idvl) = NaN;
    patch(xgl(quadl'),ygl(quadl'),dthkltf(quadl'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    dthkmtf = cthkmd_tf(:,:,i);
    %dthkmtf(~idvm) = NaN;
    patch(xgm(quadm'),ygm(quadm'),dthkmtf(quadm'),'FaceColor','interp', ...
        'EdgeColor','interp');
    axis equal;
    view(-90,90);
    title('T2S - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    % colorbar(hf5,'Ticks',-2.5:.5:2.5);
    
    clim([-1.5 1.5]);
    xlim([-30 30]);
    ylim([-40 40]);
    % colormap(hf4, 'parula(10)');
    % colormap(hf5, 'parula(10)');
    colormap(hf4, 'parula');
    colormap(hf5, 'parula');
    % scatter(xgl,ygl,4,".",'CData', grayColor);
    % scatter(xgm,ygm,4,".",'CData', grayColor);
    axlim=axis;
    plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
    plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
    plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
    plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
    tflat1 = mean(dthkltf(ilat(:,1)),"omitnan");
    tflat2 = mean(dthkltf(ilat(:,2)),"omitnan");
    tflat3 = mean(dthkltf(ilat(:,3)),"omitnan");
    tfmed1 = mean(dthkmtf(imed(:,1)),"omitnan");
    tfmed2 = mean(dthkmtf(imed(:,2)),"omitnan");
    tfmed3 = mean(dthkmtf(imed(:,3)),"omitnan");
    % text(-15,20,num2str(tflat1));
    % text(-4,20,num2str(tflat2));
    % text(7,20,num2str(tflat3));
    % text(-13,-20,num2str(tfmed1));
    % text(0,-20,num2str(tfmed2));
    % text(13,-20,num2str(tfmed3));

    pic_nam=fullfile(rdir, 'Tibial Thickness Differences.pdf');

    if i==1 && exist(pic_nam,'file')==2
        delete(pic_nam);
    end
    set(f, 'units','normalized','outerposition',[0 0 1 1]);
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true);
    close(f);


    %%
    %
    % Statistics and Plot of Differences
    %
    dthklrf = dthklrf(:);
    dthkmrf = dthkmrf(:);
    dthkltf = dthkltf(:);
    dthkmtf = dthkmtf(:);

    dmeanlrf = mean(dthklrf,"omitnan");
    dmeanmrf = mean(dthkmrf,"omitnan");
    dmeanltf = mean(dthkltf,"omitnan");
    dmeanmtf = mean(dthkmtf,"omitnan");
    dstdlrf = std(dthklrf,"omitnan");
    dstdmrf = std(dthkmrf,"omitnan");
    dstdltf = std(dthkltf,"omitnan");
    dstdmtf = std(dthkmtf,"omitnan");

    dmaxlrf = max(dthklrf);
    dminlrf = min(dthklrf);
    dmaxmrf = max(dthkmrf);
    dminmrf = min(dthkmrf);
    dmaxltf = max(dthkltf);
    dminltf = min(dthkltf);
    dmaxmtf = max(dthkmtf);
    dminmtf = min(dthkmtf);

    %

    %
    f = figure;
    t = tiledlayout(2,2);
    sgtitle({[fstr(1:5) ' - Tibia'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');
    orient landscape;
    %
    nexttile;
    %plot(dthklrf,'k.');
    hold on;
    %plot(axlim(1:2),[0 0],'k-');
    plot1 = plot(dthklrf(ilat(:,1)),'r.');
    plot2 = plot(dthklrf(ilat(:,2)),'.','Color','#00c04b');
    plot3 = plot(dthklrf(ilat(:,3)),'b.');
    ylim([-2 2]);
    axlim = axis;
    plot(axlim(1:2),[dmeanlrf dmeanlrf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf+3*dstdlrf dmeanlrf+3*dstdlrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanlrf-3*dstdlrf dmeanlrf-3*dstdlrf],'r--','LineWidth',1);
    title(['Rho/FFE Lateral Differences', newline, 'Max: ', num2str(dmaxlrf), ' Min: ', num2str(dminlrf)],...
        'FontSize',16,'FontWeight','bold');
    legend([plot1, plot2, plot3], {'Posterior', 'Center', 'Anterior'});


    nexttile;
    %plot(dthkltf,'k.');
    hold on;
    %plot(axlim(1:2),[0 0],'k-');
    plot1 = plot(dthkltf(ilat(:,1)),'r.');
    plot2 = plot(dthkltf(ilat(:,2)),'.','Color','#00c04b');
    plot3 = plot(dthkltf(ilat(:,3)),'b.');
    ylim([-2 2]);
    axlim = axis;
    plot(axlim(1:2),[dmeanltf dmeanltf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanltf+3*dstdltf dmeanltf+3*dstdltf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanltf-3*dstdltf dmeanltf-3*dstdltf],'r--','LineWidth',1);
    title(['T2S/FFE Lateral Differences', newline, 'Max: ', num2str(dmaxltf), ' Min: ', num2str(dminltf)],...
        'FontSize',16,'FontWeight','bold');
    legend([plot1, plot2, plot3], {'Posterior', 'Center', 'Anterior'});


    %
    nexttile;
    %plot(dthkmrf,'k.');
    hold on;
    %plot(axlim(1:2),[0 0],'k-');
    plot1 = plot(dthkmrf(imed(:,1)),'r.');
    plot2 = plot(dthkmrf(imed(:,2)),'.','Color','#00c04b');
    plot3 = plot(dthkmrf(imed(:,3)),'b.');
    ylim([-2 2]);
    axlim = axis;
    plot(axlim(1:2),[dmeanmrf dmeanmrf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf+3*dstdmrf dmeanmrf+3*dstdmrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmrf-3*dstdmrf dmeanmrf-3*dstdmrf],'r--','LineWidth',1);
    title(['Rho/FFE Medial Differences', newline, 'Max: ', num2str(dmaxmrf), ' Min: ', num2str(dminmrf)],...
        'FontSize',16,'FontWeight','bold');
    legend([plot1, plot2, plot3], {'Posterior', 'Center', 'Anterior'});


    %
    nexttile;
    %plot(dthkmtf,'k.');
    hold on;
    %plot(axlim(1:2),[0 0],'k-');
    plot1 = plot(dthkmtf(imed(:,1)),'r.');
    plot2 = plot(dthkmtf(imed(:,2)),'.','Color','#00c04b');
    plot3 = plot(dthkmtf(imed(:,3)),'b.');
    ylim([-2 2]);
    axlim = axis;
    plot(axlim(1:2),[dmeanmtf dmeanmtf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf+3*dstdmtf dmeanmtf+3*dstdmtf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanmtf-3*dstdmtf dmeanmtf-3*dstdmtf],'r--','LineWidth',1);
    title(['T2S/FFE Medial Differences', newline, 'Max: ', num2str(dmaxmtf), ' Min: ', num2str(dminmtf)],...
        'FontSize',16,'FontWeight','bold');
    legend([plot1, plot2, plot3], {'Posterior', 'Center', 'Anterior'});

    
    jj=0;
    % Write Thicknesses to CSV Spreadsheet
    for k = 1:6
        output = fullfile(rdir,roi(k));
        if i==1 && exist(output,'file')==2
            delete(output);
        end
     
        col_header1 = {'Subject'};
        col_header2 = {'ScanType'};
        col_header3 = {'Compartment'};
        col_header4 = {'Division'};
        
    
        if rem(k,2) == 1
            cthkf = cthklf_v;
            cthkr = cthklr_v;
            cthkt = cthklt_v;
            coord = coordl;
            roi_idx = ilat;
            jj = jj+1;
            cmprt = 0; % lateral
            if rem(jj,3) == 1
                division = 0; % posterior
            elseif rem(jj,3) == 2
                division = 1; % central
            elseif rem(jj,3) == 0
                division = 2; % anterior
            end
        elseif rem(k,2) == 0
            cthkf = cthkmf_v;
            cthkr = cthkmr_v;
            cthkt = cthkmt_v;
            coord = coordm;
            roi_idx = imed;
            cmprt = 1; % medial
            if rem(jj,3) == 1
                division = 0; % posterior
            elseif rem(jj,3) == 2
                division = 1; % central
            elseif rem(jj,3) == 0
                division = 2; % anterior
            end
        end

        if i == 1
            if size(roi_idx) < 16385 %Last Excel Column is 16,384
                col_header5 = repmat('pt_', size(coord(roi_idx(:,jj)),1),1);
                col_header5 = [col_header5 int2str(coord(roi_idx(:,jj)))];
                col_header5 = string(col_header5);
                col_header5 = col_header5';
            else
                error([' *** ERROR: DATA EXCEEDS TOTAL NUMBER OF COLUMNS', ...
                    ' AVAILABLE IN EXCEL SHEET']);
            end

            writecell(col_header1,output,'Range','A1')
            writecell(col_header2,output,'Range','B1')
            writecell(col_header3,output,'Range','C1')
            writecell(col_header4,output,'Range','D1')
            writematrix(col_header5,output, 'Range','E1')
        end
        
        
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+1+ii)])
        writematrix(0,output,'Range',['B' int2str(i+1+ii)]) %FFE
        writematrix(cmprt,output,'Range',['C' int2str(i+1+ii)])
        writematrix(division,output,'Range',['D' int2str(i+1+ii)])
        writematrix(cthkf(roi_idx(:,jj))',output,'Range',['E' int2str(i+1+ii)])
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+2+ii)])
        writematrix(1,output,'Range',['B' int2str(i+2+ii)]) %RHO
        writematrix(cmprt,output,'Range',['C' int2str(i+2+ii)])
        writematrix(division,output,'Range',['D' int2str(i+2+ii)])
        writematrix(cthkr(roi_idx(:,jj))',output,'Range',['E' int2str(i+2+ii)])
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+3+ii)])
        writematrix(2,output,'Range',['B' int2str(i+3+ii)]) %T2S
        writematrix(cmprt,output,'Range',['C' int2str(i+3+ii)])
        writematrix(division,output,'Range',['D' int2str(i+3+ii)])
        writematrix(cthkt(roi_idx(:,jj))',output,'Range',['E' int2str(i+3+ii)])
        %
    end
    set(f, 'units','normalized','outerposition',[0 0 1 1]);
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true);
    close(f);
    %
    ii = ii+2;
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

% avg_rflat(~idvl) = NaN;
% avg_tflat(~idvl) = NaN;
% std_rflat(~idvl) = NaN;
% std_tflat(~idvl) = NaN;

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

% avg_rfmed(~idvm) = NaN;
% avg_tfmed(~idvm) = NaN;
% std_rfmed(~idvm) = NaN;
% std_tfmed(~idvm) = NaN;

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


f = figure;
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
colormap(hf1,parula);
hf1.CLim = [-1.5 1.5];
colorbar;
% colorbar(hf1, 'Ticks', -2.5:.5:2.5);
xlim(hf1,[-30 30]);
ylim(hf1,[-40 40]);
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
% text(-15,20,num2str(rflat1));
% text(-4,20,num2str(rflat2));
% text(7,20,num2str(rflat3));
% text(-13,-20,num2str(rfmed1));
% text(0,-20,num2str(rfmed2));
% text(13,-20,num2str(rfmed3));

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
colormap(hf2,parula);
hf2.CLim = [-1.5 1.5];
colorbar;
% colorbar(hf2, 'Ticks', -2.5:.5:2.5);
xlim(hf2,[-30 30]);
ylim(hf2,[-40 40]);
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);
axlim=axis;
plot([reg_lat1_all reg_lat1_all],[axlim(4) 0],"k");
plot([reg_lat2_all reg_lat2_all],[axlim(4) 0],"k");
plot([reg_med1_all reg_med1_all],[0 axlim(3)],"k");
plot([reg_med2_all reg_med2_all],[0 axlim(3)],"k");
% text(-15,20,num2str(tflat1));
% text(-4,20,num2str(tflat2));
% text(7,20,num2str(tflat3));
% text(-13,-20,num2str(tfmed1));
% text(0,-20,num2str(tfmed2));
% text(13,-20,num2str(tfmed3));

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
colormap(hf3,jet);
colorbar(hf3,'Ticks', 0:.25:1.5);
xlim(hf3,[-30 30]);
ylim(hf3,[-40 40]);
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);
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
colormap(hf4,jet);
colorbar(hf4,'Ticks', 0:.25:1.5);
xlim(hf4,[-30 30]);
ylim(hf4,[-40 40]);
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);
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
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);

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
% scatter(xgl,ygl,4,".",'CData', grayColor);
% scatter(xgm,ygm,4,".",'CData', grayColor);

set(f, 'units','normalized','outerposition',[0 0 1 1]);
exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true);
close(f);
%%


%Take out NaNs
% 
% 
% 
% output = fullfile(rdir,'Tibial Thicknesses.xlsx');
% 
% col_header1 = {'RHO - FFE'};
% col_header2 = {'T2S - FFE'};
% col_header4 = {'Lateral Posterior Differences (mm)','Grid Coordinates',...
%     ' Medial Posterior Differences (mm)','Grid Coordinates','Lateral Central Differences (mm)',...
%     'Grid Coordinates', 'Medial Central Differences','Grid Coordinates',...
%     'Lateral Anterior Differences (mm)','Grid Coordinates','Medial Anterior Differences (mm)',...
%     'Grid Coordinates'};
% 
% writecell(col_header1,output,'Sheet','Difference Averages','Range','A1')
% writecell(col_header2,output,'Sheet','Difference Averages','Range','N1')
% writecell(col_header4,output,'Sheet','Difference Averages','Range','A2')
% writecell(col_header4,output,'Sheet','Difference Averages','Range','N2')
% 
% writematrix(avg_rflat(ilat(:,1)),output,'Sheet','Difference Averages','Range','A3')
% writematrix(coordl(ilat(:,1)),output,'Sheet','Difference Averages','Range','B3')
% writematrix(avg_rfmed(imed(:,1)),output,'Sheet','Difference Averages','Range','C3')
% writematrix(coordm(imed(:,1)),output,'Sheet','Difference Averages','Range','D3')
% writematrix(avg_rflat(ilat(:,2)),output,'Sheet','Difference Averages','Range','E3')
% writematrix(coordl(ilat(:,2)),output,'Sheet','Difference Averages','Range','F3')
% writematrix(avg_rfmed(imed(:,2)),output,'Sheet','Difference Averages','Range','G3')
% writematrix(coordm(imed(:,2)),output,'Sheet','Difference Averages','Range','H3')
% writematrix(avg_rflat(ilat(:,3)),output,'Sheet','Difference Averages','Range','I3')
% writematrix(coordl(ilat(:,3)),output,'Sheet','Difference Averages','Range','J3')
% writematrix(avg_rfmed(imed(:,3)),output,'Sheet','Difference Averages','Range','K3')
% writematrix(coordm(imed(:,3)),output,'Sheet','Difference Averages','Range','L3')
% 
% writematrix(avg_tflat(ilat(:,1)),output,'Sheet','Difference Averages','Range','N3')
% writematrix(coordl(ilat(:,1)),output,'Sheet','Difference Averages','Range','O3')
% writematrix(avg_tfmed(imed(:,1)),output,'Sheet','Difference Averages','Range','P3')
% writematrix(coordm(imed(:,1)),output,'Sheet','Difference Averages','Range','Q3')
% writematrix(avg_tflat(ilat(:,2)),output,'Sheet','Difference Averages','Range','R3')
% writematrix(coordl(ilat(:,2)),output,'Sheet','Difference Averages','Range','S3')
% writematrix(avg_tfmed(imed(:,2)),output,'Sheet','Difference Averages','Range','T3')
% writematrix(coordm(imed(:,2)),output,'Sheet','Difference Averages','Range','U3')
% writematrix(avg_tflat(ilat(:,3)),output,'Sheet','Difference Averages','Range','V3')
% writematrix(coordl(ilat(:,3)),output,'Sheet','Difference Averages','Range','W3')
% writematrix(avg_tfmed(imed(:,3)),output,'Sheet','Difference Averages','Range','X3')
% writematrix(coordm(imed(:,3)),output,'Sheet','Difference Averages','Range','Y3')
return