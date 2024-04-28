%#######################################################################
%
%           * Femur Cartilage THicKness Comparison Program *
%
%          M-File which reads the T1FFE and T1rho cartilage thicknesses,
%     finds a combined grid that covers both data sets and finds the
%     differences in cartilage thicknesses.
%
%          Plots, Outputs, ...
%
%     NOTES:  1.  Both grids must have integer coordinates.
%
%     17-July-2023 * Aaron Dees
%

%#######################################################################
%%
clear;
close all;
clc;
grayColor = [.7 .7 .7];

div = uigetdir;
rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
load(fullfile(rdir, 'Subject Files Full Femur.mat'));
ns=size(sd,1);
%mat_rho = fullfile(ddir_rho,'tcart08_1.mat');
%for i=1:ns


% %%
% %
% % Get Analysis Grids
%
ffeg = load(fullfile(gdir,'FFE_fgrid08_.mat'));
rhog = load(fullfile(gdir,'RHO_fgrid08_.mat'));
t2sg = load(fullfile(gdir,'T2S_fgrid08_.mat'));
% rhog = load(fullfile('TibiaCartThk_1_11Mar2022(T1Rho)', ...
%             'tgrid08_1.mat'));
%%%
%
% Get Coordinates and Ranges for Lateral Compartment
%
tq_ffe = ffeg.tq;
zq_ffe = ffeg.zq;
%
tq_rho = rhog.tq;
zq_rho = rhog.zq;
%
tq_t2s = t2sg.tq;
zq_t2s = t2sg.zq;
%
tmin_ffe = min(tq_ffe);
tmax_ffe = max(tq_ffe);
nc_ffe = ffeg.nt;         % Number of columns
%
tmin_rho = min(tq_rho);
tmax_rho = max(tq_rho);
nc_rho = rhog.nt;          % Number of columns
%
tmin_t2s = min(tq_t2s);
tmax_t2s = max(tq_t2s);
nc_t2s = t2sg.nt;
%
zmin_ffe = min(zq_ffe);
zmax_ffe = max(zq_ffe);
nr_ffe = ffeg.nz;         % Number of rows

zmin_rho = min(zq_rho);
zmax_rho = max(zq_rho);
nr_rho = rhog.nz;         % Number of rows

zmin_t2s = min(zq_t2s);
zmax_t2s = max(zq_t2s);
nr_t2s = t2sg.nz;         % Number of rows
%%

%Get Combined Grid for Both Data Sets

tmin=[tmin_ffe; tmin_rho; tmin_t2s];
tmin=min(tmin);

tmax = [tmax_ffe; tmax_rho; tmax_t2s];
tmax=max(tmax);

zmin = [zmin_ffe; zmin_rho; zmin_t2s];
zmin = min(zmin);

zmax = [zmax_ffe; zmax_rho; zmax_t2s];
zmax = max(zmax);


[tg,zg] = meshgrid(tmin:2:tmax,zmin:zmax);
[nr,nc] = size(tg);

%%
%
% Get Indexes into Combined Grid
%
n = nr*nc;           % Number of points in combined lateral grid
idx = reshape(1:n,nr,nc);
%
offstcf = (tmin_ffe-tmin)/2+1;      % Differences in integer grid == column index
offstrf = (zmin_ffe-zmin+1);      % Differences in integer grid == row index
idxf = idx(offstrf:offstrf+nr_ffe-1,offstcf:offstcf+nc_ffe-1);   % Index for T1FFE
%
offstcr = (tmin_rho-tmin)/2+1;      % Differences in integer grid == column index
offstrr = (zmin_rho-zmin+1);      % Differences in integer grid == row index
idxr = idx(offstrr:offstrr+nr_rho-1,offstcr:offstcr+nc_rho-1);   % Index for T1rho
%
offstct = (tmin_t2s-tmin)/2+1;      % Differences in integer grid == column index
offstrt = (zmin_t2s-zmin+1);      % Differences in integer grid == row index
idxt = idx(offstrt:offstrt+nr_t2s-1,offstct:offstct+nc_t2s-1);   % Index for T2*

%%

%%
%
% Read Cartilage Thicknesses
%

cthkd_rf = zeros(nr,nc,ns);
cthkd_tf = zeros(nr,nc,ns);
%%
cthkf = [];
cthkr = [];
cthkt = [];
%%
for i=1:ns
    fstr=sd(i).FFE.femur.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat'];
    ffe = load(fullfile(tdir,fstr));

    fstr=sd(i).RHO.femur.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat'];
    rho = load(fullfile(tdir,fstr));

    fstr=sd(i).T2S.femur.bfnam;
    fstr=[fstr(1:6) fstr(15:17) '_fcart08_thk.mat'];
    t2s = load(fullfile(tdir,fstr));

    cthkf1 = ffe.cthk;
    cthkr1 = rho.cthk;
    cthkt1 = t2s.cthk;
    %%
    %
    % Get Cartilage Thicknesses at Similar Coordinates
    %
    cthkf2 = NaN(nr,nc);  % NaN == missing data
    cthkr2 = NaN(nr,nc);
    cthkt2 = NaN(nr,nc);
    %
    cthkf2(idxf) = cthkf1;
    cthkr2(idxr) = cthkr1;
    cthkt2(idxt) = cthkt1;

    cthkf = cat(3,cthkf,cthkf2);
    cthkr = cat(3,cthkr,cthkr2);
    cthkt = cat(3,cthkt,cthkt2);
    %
    % Calculate Differences
    %
    cthkd_rf(:,:,i) = cthkr(:,:,i) - cthkf(:,:,i);
    cthkd_tf(:,:,i) = cthkt(:,:,i) - cthkf(:,:,i);


    % idv_rf = ~isnan(cthkd_rf(:,:,i));  % Valid differences
    % idv_tf = ~isnan(cthkd_tf(:,:,i));  % Valid differences
    % idv = idv_rf & idv_tf; % Total Valid differences

    % Only include thicknesses where valid differences exist
    % cthkf2(~idv) = NaN;
    % cthkr2(~idv) = NaN;
    % cthkt2(~idv) = NaN;
end
idx_cthkf = ~isnan(cthkf);
idx_cthkr = ~isnan(cthkr);
idx_cthkt = ~isnan(cthkt);

cols = sum(max(idx_cthkf),3) + sum(max(idx_cthkr),3) + sum(max(idx_cthkt),3);
rows = sum(max(idx_cthkf,[],2),3) + sum(max(idx_cthkr,[],2),3) + sum(max(idx_cthkt,[],2),3);
cols = logical(cols);
rows = logical(rows);

col1 = find(cols,1,'first');
col2 = find(cols,1,'last');
row1 = find(rows,1,'first');
row2 = find(rows,1,'last');

cthkf = cthkf(row1:row2,col1:col2,:);
cthkr = cthkr(row1:row2,col1:col2,:);
cthkt = cthkt(row1:row2,col1:col2,:);

cthkd_rf = cthkd_rf(row1:row2,col1:col2,:);
cthkd_tf = cthkd_tf(row1:row2,col1:col2,:);

tg = tg(row1:row2,col1:col2);
zg = zg(row1:row2,col1:col2);

[nr,nc] = size(tg);
tg = tg(:);
zg = zg(:);
quad = quadconn(nr,nc);   % Quadrilateral connectivity for grid

zsc = t2sg.zsc;
t0 = -145; % Angle (theta) cutoff (-150, -145, or -140)
[idm_ant,idl_ant,idm_ctr,idl_ctr,idm_pos,idl_pos] = freg_axpf2_AD(tg,zg,zsc,t0,sd,div);
A = reshape(1:size(tg),nr,nc);
A(~(idm_ant|idl_ant)) = NaN;
A = min(A,[],2);

roi(1) = "fem_lat_pos.xlsx";
roi(2) = "fem_med_pos.xlsx";
roi(3) = "fem_lat_ctr.xlsx";
roi(4) = "fem_med_ctr.xlsx";
roi(5) = "fem_lat_ant.xlsx";
roi(6) = "fem_med_ant.xlsx";

roi_idx(:,1) = idl_pos;
roi_idx(:,2) = idm_pos;
roi_idx(:,3) = idl_ctr;
roi_idx(:,4) = idm_ctr;
roi_idx(:,5) = idl_ant;
roi_idx(:,6) = idm_ant;

coord = 1:size(tg);
coord = coord(:);

ii=0;
for i=1:ns
    %%
    % Turning the thicknesses into vectors
    cthkf_v = cthkf(:,:,i);
    cthkr_v = cthkr(:,:,i);
    cthkt_v = cthkt(:,:,i);
    cthkf_v = cthkf_v(:);
    cthkr_v = cthkr_v(:);
    cthkt_v = cthkt_v(:);

    % dthkrf = cthkd_rf(:,:,i);
    % dthktf = dthktf(:);

    % Get averages for specific ROIs
    avg_lats(1,1) = mean(cthkf_v(idl_pos),"omitnan");
    avg_lats(2,1) = mean(cthkr_v(idl_pos),"omitnan");
    avg_lats(3,1) = mean(cthkt_v(idl_pos),"omitnan");

    avg_lats(1,2) = mean(cthkf_v(idl_ctr),"omitnan");
    avg_lats(2,2) = mean(cthkr_v(idl_ctr),"omitnan");
    avg_lats(3,2) = mean(cthkt_v(idl_ctr),"omitnan");

    avg_lats(1,3) = mean(cthkf_v(idl_ant),"omitnan");
    avg_lats(2,3) = mean(cthkr_v(idl_ant),"omitnan");
    avg_lats(3,3) = mean(cthkt_v(idl_ant),"omitnan");

    avg_meds(1,1) = mean(cthkf_v(idm_pos),"omitnan");
    avg_meds(2,1) = mean(cthkr_v(idm_pos),"omitnan");
    avg_meds(3,1) = mean(cthkt_v(idm_pos),"omitnan");

    avg_meds(1,2) = mean(cthkf_v(idm_ctr),"omitnan");
    avg_meds(2,2) = mean(cthkr_v(idm_ctr),"omitnan");
    avg_meds(3,2) = mean(cthkt_v(idm_ctr),"omitnan");

    avg_meds(1,3) = mean(cthkf_v(idm_ant),"omitnan");
    avg_meds(2,3) = mean(cthkr_v(idm_ant),"omitnan");
    avg_meds(3,3) = mean(cthkt_v(idm_ant),"omitnan");

    fstr=sd(i).FFE.femur.bfnam;
    f = figure;
    tiledlayout(2,6);
    sgtitle({[fstr(1:5) ' - Femur'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');

    hf1=nexttile(1,[1 2]);
    cthk = cthkr_v;
    %cthk(~idvl_rf) = NaN;
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    axlim = axis;
    view(-90,90);
    title('T1\rho Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf1, 'Ticks',0:6);
    clim([0 6]);
    xlim(hf1,[tmin tmax]);
    ylim(hf1,[zmin zmax]);
    plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
    plot([-145 -145],[axlim(3) axlim(4)],"k");
    plot([axlim(1) axlim(2)],[0 0],"k");
    plot(tg(A)-1,zg(A),"k");
    text(-170,20,num2str(avg_lats(2,1)));
    text(-120,20,num2str(avg_lats(2,2)));
    text(-45,20,num2str(avg_lats(2,3)));
    text(-170,-20,num2str(avg_meds(2,1)));
    text(-120,-20,num2str(avg_meds(2,2)));
    text(-45,-20,num2str(avg_meds(2,3)));


    hf2=nexttile(3,[1 2]);
    cthk = cthkf_v;
    %cthk(~idvl_rf) = NaN;
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    axlim = axis;
    view(-90,90);
    title('T1FFE Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf2, 'Ticks',0:6);
    clim([0 6]);
    xlim(hf2,[tmin tmax]);
    ylim(hf2,[zmin zmax]);
    plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
    plot([-145 -145],[axlim(3) axlim(4)],"k");
    plot([axlim(1) axlim(2)],[0 0],"k");
    plot(tg(A)-1,zg(A),"k");
    text(-170,20,num2str(avg_lats(1,1)));
    text(-120,20,num2str(avg_lats(1,2)));
    text(-45,20,num2str(avg_lats(1,3)));
    text(-170,-20,num2str(avg_meds(1,1)));
    text(-120,-20,num2str(avg_meds(1,2)));
    text(-45,-20,num2str(avg_meds(1,3)));

    hf3=nexttile(5,[1 2]);
    cthk = cthkt_v;
    %cthk(~idvl_tf) = NaN;
    patch(tg(quad'),zg(quad'),cthk(quad'),'FaceColor','interp', ...
        'EdgeColor','interp');
    % thk_max=max(max(cthk(quadl')));
    % disp(['Max Lat of ',fstr, ' = ', num2str(thk_max), newline]);
    hold on;
    axlim = axis;
    view(-90,90);
    title('T2S Cartilage Thicknesses','FontSize',16,'FontWeight','bold');
    colorbar(hf3, 'Ticks',0:6);
    clim([0 6]);
    xlim(hf3,[tmin tmax]);
    ylim(hf3,[zmin zmax]);
    plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
    plot([-145 -145],[axlim(3) axlim(4)],"k");
    plot([axlim(1) axlim(2)],[0 0],"k");
    plot(tg(A)-1,zg(A),"k");
    text(-170,20,num2str(avg_lats(3,1)));
    text(-120,20,num2str(avg_lats(3,2)));
    text(-45,20,num2str(avg_lats(3,3)));
    text(-170,-20,num2str(avg_meds(3,1)));
    text(-120,-20,num2str(avg_meds(3,2)));
    text(-45,-20,num2str(avg_meds(3,3)));

    colormap(hf1, 'jet(12)');
    colormap(hf2, 'jet(12)');
    colormap(hf3, 'jet(12)');

    hf4=nexttile(8,[1 2]);
    dthkrf = cthkd_rf(:,:,i);
    % dthkrf(~idv) = NaN;
    patch(tg(quad'),zg(quad'),dthkrf(quad'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    axlim = axis;
    view(-90,90);
    title('T1\rho - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colorbar(hf4,'Ticks',-2.5:.5:2.5);
    clim([-2.5 2.5]);
    xlim(hf4,[tmin tmax]);
    ylim(hf4,[zmin zmax]);
    plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
    plot([-145 -145],[axlim(3) axlim(4)],"k");
    plot([axlim(1) axlim(2)],[0 0],"k");
    plot(tg(A)-1,zg(A),"k");
    rflat1 = mean(dthkrf(idl_pos),"omitnan");
    rflat2 = mean(dthkrf(idl_ctr),"omitnan");
    rflat3 = mean(dthkrf(idl_ant),"omitnan");
    rfmed1 = mean(dthkrf(idm_pos),"omitnan");
    rfmed2 = mean(dthkrf(idm_ctr),"omitnan");
    rfmed3 = mean(dthkrf(idm_ant),"omitnan");
    text(-170,20,num2str(rflat1));
    text(-120,20,num2str(rflat2));
    text(-45,20,num2str(rflat3));
    text(-170,-20,num2str(rfmed1));
    text(-120,-20,num2str(rfmed2));
    text(-45,-20,num2str(rfmed3));

    hf5=nexttile(10,[1 2]);
    dthktf = cthkd_tf(:,:,i);
    % dthktf(~idv) = NaN;
    patch(tg(quad'),zg(quad'),dthktf(quad'),'FaceColor','interp', ...
        'EdgeColor','interp');
    hold on;
    axlim = axis;
    view(-90,90);
    title('T2S - T1FFE Cartilage Thickness Differences', ...
        'FontSize',16,'FontWeight','bold');
    colorbar(hf5,'Ticks',-2.5:.5:2.5);
    clim([-2.5 2.5]);
    xlim(hf5,[tmin tmax]);
    ylim(hf5,[zmin zmax]);
    colormap(hf4, 'parula(10)');
    colormap(hf5, 'parula(10)');
    plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
    plot([-145 -145],[axlim(3) axlim(4)],"k");
    plot([axlim(1) axlim(2)],[0 0],"k");
    plot(tg(A)-1,zg(A),"k");
    tflat1 = mean(dthktf(idl_pos),"omitnan");
    tflat2 = mean(dthktf(idl_ctr),"omitnan");
    tflat3 = mean(dthktf(idl_ant),"omitnan");
    tfmed1 = mean(dthktf(idm_pos),"omitnan");
    tfmed2 = mean(dthktf(idm_ctr),"omitnan");
    tfmed3 = mean(dthktf(idm_ant),"omitnan");
    text(-170,20,num2str(tflat1));
    text(-120,20,num2str(tflat2));
    text(-45,20,num2str(tflat3));
    text(-170,-20,num2str(tfmed1));
    text(-120,-20,num2str(tfmed2));
    text(-45,-20,num2str(tfmed3));

    pic_nam=fullfile(rdir, 'Femoral Thickness Differences.pdf');

    if i==1 && exist(pic_nam,'file')==2
        delete(pic_nam);
    end
    set(f,'units','normalized','outerposition',[0 0 1 1]);
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true);
    close(f);

    %%
    %
    % Statistics and Plot of Differences
    %
    dthkrf = dthkrf(:);
    dthktf = dthktf(:);
    dmeanrf = mean(dthkrf,"omitnan");
    dmeantf = mean(dthktf,"omitnan");
    dstdrf = std(dthkrf,"omitnan");
    dstdtf = std(dthktf,"omitnan");

    dmaxrf = max(dthkrf);
    dminrf = min(dthkrf);
    dmaxtf = max(dthktf);
    dmintf = min(dthktf);

    %
    f = figure;
    tiledlayout(2,1);
    sgtitle({[fstr(1:5) ' - Femur'];}, ...
        'FontSize',16,'FontWeight','bold','Interpreter','none');

    %
    nexttile;
    plot(dthkrf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeanrf dmeanrf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeanrf+3*dstdrf dmeanrf+3*dstdrf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeanrf-3*dstdrf dmeanrf-3*dstdrf],'r--','LineWidth',1);
    title('Rho/FFE Differences', ...
        'FontSize',16,'FontWeight','bold');


    nexttile;
    plot(dthktf,'k.');
    hold on;
    axlim = axis;
    %plot(axlim(1:2),[0 0],'k-');
    plot(axlim(1:2),[dmeantf dmeantf],'b--','LineWidth',1);
    plot(axlim(1:2),[dmeantf+3*dstdtf dmeantf+3*dstdtf],'r--','LineWidth',1);
    plot(axlim(1:2),[dmeantf-3*dstdtf dmeantf-3*dstdtf],'r--','LineWidth',1);
    title('T2S/FFE Differences', ...
        'FontSize',16,'FontWeight','bold');

    jj = 0;
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
            jj = jj+1;
            cmprt = 0;
            if rem(jj,3) == 1
                division = 0;
            elseif rem(jj,3) == 2
                division = 1;
            elseif rem(jj,3) == 0
                division = 2;
            end
        elseif rem(k,2) == 0
            cmprt = 1;
            if rem(jj,3) == 1
                division = 0;
            elseif rem(jj,3) == 2
                division = 1;
            elseif rem(jj,3) == 0
                division = 2;
            end
        end

        if i == 1
            if sum(roi_idx(:,k)) < 16385 %Last Excel Column is 16,384
                col_header5 = repmat('pt_', size(coord(roi_idx(:,k)),1),1);
                col_header5 = [col_header5 int2str(coord(roi_idx(:,k)))];
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
        writematrix(0,output,'Range',['B' int2str(i+1+ii)])
        writematrix(cmprt,output,'Range',['C' int2str(i+1+ii)])
        writematrix(division,output,'Range',['D' int2str(i+1+ii)])
        writematrix(cthkf_v(roi_idx(:,k))',output,'Range',['E' int2str(i+1+ii)])
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+2+ii)])
        writematrix(1,output,'Range',['B' int2str(i+2+ii)])
        writematrix(cmprt,output,'Range',['C' int2str(i+2+ii)])
        writematrix(division,output,'Range',['D' int2str(i+2+ii)])
        writematrix(cthkr_v(roi_idx(:,k))',output,'Range',['E' int2str(i+2+ii)])
        writematrix(str2double(fstr(1:3)) ,output,'Range',['A' int2str(i+3+ii)])
        writematrix(2,output,'Range',['B' int2str(i+3+ii)])
        writematrix(cmprt,output,'Range',['C' int2str(i+3+ii)])
        writematrix(division,output,'Range',['D' int2str(i+3+ii)])
        writematrix(cthkt_v(roi_idx(:,k))',output,'Range',['E' int2str(i+3+ii)]) 
        %
    end
    set(f,'units','normalized','outerposition',[0 0 1 1]);
    exportgraphics(f, pic_nam, "Resolution", 300, 'Append', true);
    close(f);
    %
    ii = ii+2;
end

%%
n_rf = sum(~isnan(cthkd_rf),3);
idxn = find(n_rf<3);
n_rf(idxn) = NaN;
avg_rf = mean(cthkd_rf,3,"omitnan");
avg_rf(idxn) = NaN;
std_rf = std(cthkd_rf,0,3,"omitnan");
std_rf(idxn) = NaN;

n_tf = sum(~isnan(cthkd_tf),3);
idxn = find(n_tf<3);
n_tf(idxn) = NaN;
avg_tf = mean(cthkd_tf,3,"omitnan");
avg_tf(idxn) = NaN;
std_tf = std(cthkd_tf,0,3,"omitnan");
std_tf(idxn) = NaN;

rflat1 = mean(avg_rf(idl_pos),"omitnan");
rflat2 = mean(avg_rf(idl_ctr),"omitnan");
rflat3 = mean(avg_rf(idl_ant),"omitnan");
tflat1 = mean(avg_tf(idl_pos),"omitnan");
tflat2 = mean(avg_tf(idl_ctr),"omitnan");
tflat3 = mean(avg_tf(idl_ant),"omitnan");
rfmed1 = mean(avg_rf(idm_pos),"omitnan");
rfmed2 = mean(avg_rf(idm_ctr),"omitnan");
rfmed3 = mean(avg_rf(idm_ant),"omitnan");
tfmed1 = mean(avg_tf(idm_pos),"omitnan");
tfmed2 = mean(avg_tf(idm_ctr),"omitnan");
tfmed3 = mean(avg_tf(idm_ant),"omitnan");

% avg_rf(~idv) = NaN;
% avg_tf(~idv) = NaN;
% std_rf(~idv) = NaN;
% std_tf(~idv) = NaN;

% 0 - 10 Colormap for Number of Subjects
%cmap=[0 0 0.6; 0 0 1; 0 0.4 1; 0 0.8 1; 0.2 1 0.8; 0.6 1 0.4;
%1 1 0; 1 0.6 0; 1 0.4 0; 1 0 0; 0.6 0 0];

% 3 - 10 Colormap for Number of Subjects
cmap=[0 0.8 1; 0.2 1 0.8; 0.6 1 0.4; 1 1 0; 1 0.6 0; 1 0.4 0;
    1 0 0; 0.6 0 0];


f = figure;
t = tiledlayout(3,2);
title(t,'Femoral Cartilage Thickness Differences','FontSize',20,'FontWeight','bold');


hf1 = nexttile(1);
patch(tg(quad'),zg(quad'),avg_rf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('Rho - T1FFE Averages', ...
    'FontSize',14,'FontWeight','bold');
colormap(hf1,parula(10));
hf1.CLim = [-2.5 2.5];
colorbar(hf1, 'Ticks', -2.5:.5:2.5);
xlim(hf1,[tmin tmax]);
ylim(hf1,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf1.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");
text(-170,20,num2str(rflat1));
text(-120,20,num2str(rflat2));
text(-45,20,num2str(rflat3));
text(-170,-20,num2str(rfmed1));
text(-120,-20,num2str(rfmed2));
text(-45,-20,num2str(rfmed3));

hf2 = nexttile(2);
patch(tg(quad'),zg(quad'),avg_tf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('T2S - T1FFE Averages', ...
    'FontSize',14,'FontWeight','bold');
colormap(hf2,parula(10));
hf2.CLim = [-2.5 2.5];
colorbar(hf2, 'Ticks', -2.5:.5:2.5);
xlim(hf2,[tmin tmax]);
ylim(hf2,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf2.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");
text(-170,20,num2str(tflat1));
text(-120,20,num2str(tflat2));
text(-45,20,num2str(tflat3));
text(-170,-20,num2str(tfmed1));
text(-120,-20,num2str(tfmed2));
text(-45,-20,num2str(tfmed3));

hf3 = nexttile(3);
patch(tg(quad'),zg(quad'),std_rf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('Rho - T1FFE','Standard Deviatons', ...
    'FontSize',14,'FontWeight','bold');
hf3.CLim = [0 1.5];
colormap(hf3,jet(6));
colorbar(hf3,'Ticks', 0:.25:1.5);
xlim(hf3,[tmin tmax]);
ylim(hf3,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf3.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");

hf4 = nexttile(4);
patch(tg(quad'),zg(quad'),std_tf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('T2S - T1FFE','Standard Deviatons', ...
    'FontSize',14,'FontWeight','bold');
hf4.CLim = [0 1.5];
colormap(hf4,jet(6));
colorbar(hf4,'Ticks', 0:.25:1.5);
xlim(hf4,[tmin tmax]);
ylim(hf4,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf4.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");

hf5 = nexttile(5);
patch(tg(quad'),zg(quad'),n_rf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('Rho - T1FFE','Subject Numbers', ...
    'FontSize',14,'FontWeight','bold');
%hf5.CLim=[0 11];
hf5.CLim=[2 10];
colormap(hf5,cmap);
%colorbar(hf5,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
colorbar(hf5,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
xlim(hf5,[tmin tmax]);
ylim(hf5,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf5.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");

hf6 = nexttile(6);
patch(tg(quad'),zg(quad'),n_tf(quad'),'FaceColor','interp', ...
    'EdgeColor','interp');
hold on;
view(-90,90);
axlim = axis;
title('T2S - T1FFE','Subject Numbers', ...
    'FontSize',14,'FontWeight','bold');
%hf6.CLim=[0 11];
hf6.CLim=[2 10];
colormap(hf6,cmap);
%colorbar(hf6,'Ticks', 0.5:10.5, 'TickLabels', ["0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
colorbar(hf6,'Ticks', 2.5:9.5, 'TickLabels', ["3" "4" "5" "6" "7" "8" "9" "10"], 'TickLength', 0);
xlim(hf6,[-210 30]);
ylim(hf6,[zmin zmax]);
plot(tg,zg,".",'MarkerSize',1,'Color', grayColor);
hf6.PlotBoxAspectRatio = [0.9855 1 0.9855];
plot([-145 -145],[axlim(3) axlim(4)],"k");
plot([axlim(1) axlim(2)],[0 0],"k");
plot(tg(A)-1,zg(A),"k");

set(f,'units','normalized','outerposition',[0 0 1 1]);
exportgraphics(f,pic_nam, "Resolution", 300, 'Append', true);
close(f);
%%
%
% col_header1 = {'RHO - FFE'};
% col_header2 = {'T2S - FFE'};
% col_header4 = {'Lateral Posterior Differences (mm)','Grid Coordinates',...
%     ' Medial Posterior Differences (mm)','Grid Coordinates','Lateral Central Differences (mm)',...
%     'Grid Coordinates', 'Medial Central Differences','Grid Coordinates',...
%     'Lateral Anterior Differences (mm)','Grid Coordinates','Medial Anterior Differences (mm)',...
%     'Grid Coordinates'};
% writecell(col_header1,output,'Sheet','Difference Averages','Range','A1')
% writecell(col_header2,output,'Sheet','Difference Averages','Range','N1')
% writecell(col_header4,output,'Sheet','Difference Averages','Range','A2')
% writecell(col_header4,output,'Sheet','Difference Averages','Range','N2')
%
% writematrix(avg_rf(idl_pos),output,'Sheet','Difference Averages','Range','A3')
% writematrix(coord(idl_pos),output,'Sheet','Difference Averages','Range','B3')
% writematrix(avg_rf(idm_pos),output,'Sheet','Difference Averages','Range','C3')
% writematrix(coord(idm_pos),output,'Sheet','Difference Averages','Range','D3')
% writematrix(avg_rf(idl_ctr),output,'Sheet','Difference Averages','Range','E3')
% writematrix(coord(idl_ctr),output,'Sheet','Difference Averages','Range','F3')
% writematrix(avg_rf(idm_ctr),output,'Sheet','Difference Averages','Range','G3')
% writematrix(coord(idm_ctr),output,'Sheet','Difference Averages','Range','H3')
% writematrix(avg_rf(idl_ant),output,'Sheet','Difference Averages','Range','I3')
% writematrix(coord(idl_ant),output,'Sheet','Difference Averages','Range','J3')
% writematrix(avg_rf(idm_ant),output,'Sheet','Difference Averages','Range','K3')
% writematrix(coord(idm_ant),output,'Sheet','Difference Averages','Range','L3')
%
% writematrix(avg_tf(idl_pos),output,'Sheet','Difference Averages','Range','N3')
% writematrix(coord(idl_pos),output,'Sheet','Difference Averages','Range','O3')
% writematrix(avg_tf(idm_pos),output,'Sheet','Difference Averages','Range','P3')
% writematrix(coord(idm_pos),output,'Sheet','Difference Averages','Range','Q3')
% writematrix(avg_tf(idl_ctr),output,'Sheet','Difference Averages','Range','R3')
% writematrix(coord(idl_ctr),output,'Sheet','Difference Averages','Range','S3')
% writematrix(avg_tf(idm_ctr),output,'Sheet','Difference Averages','Range','T3')
% writematrix(coord(idm_ctr),output,'Sheet','Difference Averages','Range','U3')
% writematrix(avg_tf(idl_ant),output,'Sheet','Difference Averages','Range','V3')
% writematrix(coord(idl_ant),output,'Sheet','Difference Averages','Range','W3')
% writematrix(avg_tf(idm_ant),output,'Sheet','Difference Averages','Range','X3')
% writematrix(coord(idm_ant),output,'Sheet','Difference Averages','Range','Y3')


return