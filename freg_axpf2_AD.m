function [idm_ant,idl_ant,idm_ctr,idl_ctr,idm_pos,idl_pos] = freg_axpf2_AD(tq,zq,zsc,tmin,sd,div)
%#######################################################################
%
%         * Femur REGions using AXial Plane Female 2 Program *
%
%          M-File which reads the female femur "standard" grid and
%     divides the femur into six regions based on an axial plane view of
%     the trochlea.  The regions are plotted and the indices to each
%     region are saved in the standard grid MAT File:  fb_rngf.mat.
%
%
%     NOTES:  1.  Femur MAT files *sc.mat and tibia MAT files *CS.mat
%             must be in the directory "BM_Female_Femur_Mat_Files".  The
%             MAT file fb_rngf.mat must be in the current directory.
%
%             2.  The M-files fem_plan.m, plnorm.m, plsect.m, plsect2.m,
%             and plt_datsl.m must be in the current path or directory.
%
%     28-Jul-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Plots?
%
iplt = false;
% iplt = false;
% iprt = true;
iprt = false;
%
deg = char(176);        % Degree symbol character for X axis labels
%
% Convert Radians to Degrees and Coordinate Tolerance
%
tol = 1e-10;            % Tolerance on grid coordinates
% tmin = -145;            % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0;                 % Y cutoff (-1, 0, or 1)
%
% Plot Colors and Symbols
%
clr = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 0.5 0; ...
    0.5 0 0; 0 0 0.5; 0.5 0 0.5; 0.5 0.5 0.5; 0 0.5 0.5];
clr = [clr; flipud(clr)];
sym = ['o'; 'x'; '+'; 's'; '^'; 'v'; 'p'; '*'; '<'; '>'; 'd'; 'h'; ...
    's'];
sym = [sym; sym];


rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');

ns=size(sd,1);
%
%
% Data Directory and File Name
%

%%
% Get Standard Grid
%
nq = size(tq,1);
%
% Posterior and Medial/Lateral Cutoffs
%

%
% Get Logical Indices to the Posterior and Medial/Lateral Regions
%
ipos = tq<=tmin;        % Posterior
idl = zq>=y0;           % Lateral
idm = zq<y0;            % Medial
%
% Get Unique Z Values and Index into Theta Values (in Columns)
%
zu = unique(zq);
nu = size(zu,1);        % Number of unique Z values in grid
%
izl = zu>=y0;           % Logical index to lateral Z values
nl = sum(izl);          % Number of unique lateral Z values
%
izm = ~izl;             % Logical index to medial Z values
nm = sum(izm);          % Number of unique medial Z values
%
nid = reshape(1:nq,nu,nq/nu)';         % Number of grid points (unique Z values in rows/unique Theta values in columns)
%
% Initialize Figure with First Anterior Points (Boundary) for All
% Subjects
%
if iplt
    hf = figure;
    orient landscape;
    plot(tq,zq,'k.','MarkerSize',6,'LineWidth',0.5);
    hold on;
    %
    hb = gobjects(ns,1);    % Handles to anterior plots
    %
end
%
% Loop through MAT Files
%
idal = zeros(nl,ns);    % Index to first lateral anterior theta values
idam = zeros(nm,ns);    % Index to first medial anterior theta values
%

for i = 1:ns           % Loop through MAT files

    fstr=sd(i).FFE.femur.bfnam;
    fstr=[fstr(1:6) fstr(15:17)];


    fnams(i,:)=fstr;
    load(fullfile(bdir,[fstr '_femurCS.mat']),'datfb');
    %
    % Get Trochlea Planes
    %
    [~,~,~,~,ppl,pnl,ppm,pnm] = fem_plan_AD(datfb,tmin,iplt,fstr);
    %
    if iplt
        fh = get(0,'Child');              % Get figure handles
        [~,idmx] = max(cell2mat(get(fh,'Number'))); % Find maximum figure number
        delete(fh(idmx));  % Delete posterior-central plane plot
    end
    %
    % Get Trochlea Points for Each Slice
    %
    nsl = size(datfb,1);
    ip = NaN(nsl,3);
    % idx = [];
    for ks = 1:nsl
        ipchk = plsect2(ppl,pnl,datfb{ks});
        if isempty(ipchk)||ipchk(:,2)<y0
            ipchk = plsect2(ppm,pnm,datfb{ks});
            if isempty(ipchk)||ipchk(:,2)>y0
                continue;
            else
                ip(ks,:) = ipchk;
            end
        else
            ip(ks,:) = ipchk;
        end
    end
    %
    if iplt
        plot3(ip(:,1),ip(:,2),ip(:,3),'ro');
    end
    %
    % Convert to Polar Coordinates
    %
    [theta,r,z] = cart2pol(ip(:,1),ip(:,3),ip(:,2));
    %
    % Trap for Angles Greater Than 90 degrees (pi/2 radians and Z-Axis)
    %
    id = find(theta>pi/2);
    if ~isempty(id)
        theta(id) = theta(id)-2*pi;
    end
    %
    theta = rad2deg(theta);              % Convert to degrees
    ipp = [theta r z];   % Combine coordinates in an array
    %
    % Scale "Z" Coordinates
    %
    zqs = zsc(i,1)*zq;
    %
    % Get Unique Scaled Z Values in the Lateral and Medial Compartments
    %
    zus = zqs(nid(1,:));
    zul = zus(izl);      % Lateral compartment
    zum = zus(izm);      % Medial compartment
    %
    % Use Linear Interpolation to Get Theta Cutoffs
    %
    idxl = ipp(:,3)>=y0;
    idxm = ipp(:,3)<y0;
    tcutl = interp1(ipp(idxl,3),ipp(idxl,1),zul,'linear','extrap');
    tcutm = interp1(ipp(idxm,3),ipp(idxm,1),zum,'linear','extrap');
    %
    % Find Anterior/Trochlea Grid Points on Scaled Grid
    %
    if iplt
        figure;
        orient landscape;
        plot(tq,zqs,'k.','MarkerSize',6,'LineWidth',0.5);
        hold on;
        plot(ipp(:,1),ipp(:,3),'rs');
        plot(tcutl,zul,'bo');
        plot(tcutm,zum,'mo');
    end
    %
    for kz = 1:nu
        idx = nid(:,kz);
        zz = zqs(idx(1));
        if zz<0
            im = zz==zum;
            idc = tq(idx)>tcutm(im);
            idc = idx(idc);
            idam(im,i) = idc(1);
        else
            il = zz==zul;
            idc = tq(idx)>tcutl(il);
            idc = idx(idc);
            idal(il,i) = idc(1);
        end
        if iplt
            plot(tq(idc),zqs(idc),'go');
            plot(tq(idc(1)),zqs(idc(1)),'bs','MarkerSize',7);
        end
    end
    %
    if iplt
        xlabel(['\Theta (' deg ')'],'FontSize',12,'FontWeight','bold');
        ylabel('Z (mm)','FontSize',12,'FontWeight','bold');
        title({'Anterior Grid Points'; fstr},'FontSize',16, ...
            'FontWeight','bold','interpreter','none');
    end
    %
    % Plot Boundary on Figure for All Subjects
    %
    if iplt
        figure(hf);
        id = idal(:,i);     % Lateral boundary
        hb(i) = plot(tq(id),zq(id),['k' sym(i)],'Color',clr(i,:), ...
            'LineWidth',0.5);
        id = idam(:,i);     % Medial boundary
        plot(tq(id),zq(id),['k' sym(i)],'Color',clr(i,:), ...
            'LineWidth',0.5);
    end
    %
end
%
% Get Average First Anterior Points of All Subjects
%
idala = round(mean(idal,2));
idama = round(mean(idam,2));
%
% Get Anterior/Trochlea Grid Points of All Subjects
%
zul = zu(izl);          % Lateral Z values
zum = zu(izm);          % Medial Z values
%
ial = cell(nl,1);       % Index to anterior lateral points
iam = cell(nm,1);       % Index to anterior medial points
%
for kz = 1:nu           % Loop through Z values
    idx = nid(:,kz);
    zz = zq(idx(1));
    if zz<0
        im = zz==zum;
        idz = idx>=idama(im);
        iam{im} = nid(idz,kz);
    else
        il = zz==zul;
        idz = idx>=idala(il);
        ial{il} = nid(idz,kz);
    end
end
%
idl_ant = false(nq,1);  % Logical index to no grid points
idm_ant = idl_ant;
idl_ant(cell2mat(ial)) = true;         % Lateral anterior grid points
idm_ant(cell2mat(iam)) = true;         % Medial anterior grid points
iant = idl_ant|idm_ant; % Anterior/trochlea grid points
%
% Plot Average Anterior Grid Points
%
if iplt
    %
    figure(hf);
    %
    % Plot Anterior/Trochlea Points
    %
    plot(tq(iant),zq(iant),'r.','MarkerSize',6,'LineWidth',0.5);
    %
    % Plot Legend, Axis Labels, and Title
    %
    hl = legend(hb,fnams(:,1:5));
    set(hl,'Interpreter','none','Location','East');
    xlabel(['\Theta (' deg ')'],'FontSize',12,'FontWeight','bold');
    ylabel('Z (mm)','FontSize',12,'FontWeight','bold');
    title('Trochlea Grid Points Shown in Red','FontSize',16, ...
        'FontWeight','bold');
    if iprt
        print('-dpsc2','-r300','freg_axf2.ps');
    end
end
%
% Get Remaining Indices to the Six (6) Regions
%
ictr = true(nq,1);
ictr(iant) = false;
ictr(ipos) = false;
%
idl_pos = idl&ipos;     % Lateral posterior region
idl_ctr = idl&ictr;     % Lateral center region
%
idm_pos = idm&ipos;     % Medial posterior region
idm_ctr = idm&ictr;     % Medial center region
%
% Save Indices to Standard Grid MAT File
%
% save fb_rngf.mat -append idl_ant idl_ctr idl_pos idm_ant idm_ctr ...
%     idm_pos;
%
% Plot Regions
%
if iplt
    figure;
    orient landscape;
    plot(tq(idl_ant),zq(idl_ant),'ro','MarkerSize',4,'LineWidth',0.5);
    hold on;
    plot(tq(idm_ant),zq(idm_ant),'mo','MarkerSize',4,'LineWidth',0.5);
    plot(tq(idl_ctr),zq(idl_ctr),'gs','MarkerSize',4,'LineWidth',0.5);
    plot(tq(idm_ctr),zq(idm_ctr),'bs','MarkerSize',4,'LineWidth',0.5);
    plot(tq(idl_pos),zq(idl_pos),'y^','MarkerSize',4,'LineWidth',0.5);
    plot(tq(idm_pos),zq(idm_pos),'c^','MarkerSize',4,'LineWidth',0.5);
    axis tight;
    xlabel({['\Theta (' deg ')'];['\leftarrow posterior / anterior', ...
        ' \rightarrow']},'FontSize',12,'FontWeight','bold');
    ylabel({'Z (mm)';'\leftarrow medial / lateral \rightarrow'}, ...
        'FontSize',12,'FontWeight','bold');
    title('Average of All Femurs','FontSize',16,'FontWeight','bold');
    if iprt
        print('-dpsc2','-r300','-append','freg_axf2.ps');
    end
end
% 
% If need regions for all scan types
%
% if j == 1
% ffe_idm_ant = idm_ant;
% ffe_idl_ant = idl_ant;
% ffe_idm_ctr = idm_ctr;
% ffe_idl_ctr = idl_ctr;
% ffe_idm_pos = idm_pos;
% ffe_idl_pos = idl_pos;
% elseif j == 2
%     rho_idm_ant = idm_ant;
%     rho_idl_ant = idl_ant;
%     rho_idm_ctr = idm_ctr;
%     rho_idl_ctr = idl_ctr;
%     rho_idm_pos = idm_pos;
%     rho_idl_pos = idl_pos;
% elseif j == 3
%     t2s_idm_ant = idm_ant;
%     t2s_idl_ant = idl_ant;
%     t2s_idm_ctr = idm_ctr;
%     t2s_idl_ctr = idl_ctr;
%     t2s_idm_pos = idm_pos;
%     t2s_idl_pos = idl_pos;
% end

%
% Difference between ROI labeled points based on different scan types
%

% FR_med_ant = find(abs(ffe_idm_ant-rho_idm_ant)>0);
% FR_lat_ant = find(abs(ffe_idl_ant-rho_idl_ant)>0);
% FR_med_ctr = find(abs(ffe_idm_ctr-rho_idm_ctr)>0);
% FR_lat_ctr = find(abs(ffe_idl_ctr-rho_idl_ctr)>0);
% FR_med_pos = find(abs(ffe_idm_pos-rho_idm_pos)>0);
% FR_lat_pos = find(abs(ffe_idl_pos-rho_idl_pos)>0);
%
% FT_med_ant = find(abs(ffe_idm_ant-t2s_idm_ant)>0);
% FT_lat_ant = find(abs(ffe_idl_ant-t2s_idl_ant)>0);
% FT_med_ctr = find(abs(ffe_idm_ctr-t2s_idm_ctr)>0);
% FT_lat_ctr = find(abs(ffe_idl_ctr-t2s_idl_ctr)>0);
% FT_med_pos = find(abs(ffe_idm_pos-t2s_idm_pos)>0);
% FT_lat_pos = find(abs(ffe_idl_pos-t2s_idl_pos)>0);
%
% RT_med_ant = find(abs(rho_idm_ant-t2s_idm_ant)>0);
% RT_lat_ant = find(abs(rho_idl_ant-t2s_idl_ant)>0);
% RT_med_ctr = find(abs(rho_idm_ctr-t2s_idm_ctr)>0);
% RT_lat_ctr = find(abs(rho_idl_ctr-t2s_idl_ctr)>0);
% RT_med_pos = find(abs(rho_idm_pos-t2s_idm_pos)>0);
% RT_lat_pos = find(abs(rho_idl_pos-t2s_idl_pos)>0);
return