function [idm_ant,idl_ant,idm_ctr,idl_ctr,idm_pos,idl_pos] = freg_axpf2_AD(tq,zq,zsc,tmin,sd,div) % Define function to divide femur into six regions
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
iplt = false; % Set plotting flag to false
% iplt = false;
% iprt = true;
iprt = false; % Set printing flag to false
%
deg = char(176); % Define degree symbol character
%
% Convert Radians to Degrees and Coordinate Tolerance
%
tol = 1e-10; % Set coordinate tolerance
% tmin = -145;            % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0; % Set Y cutoff
%
% Plot Colors and Symbols
%
clr = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 0.5 0; ...
    0.5 0 0; 0 0 0.5; 0.5 0 0.5; 0.5 0.5 0.5; 0 0.5 0.5]; % Define color array
clr = [clr; flipud(clr)]; % Append flipped colors
sym = ['o'; 'x'; '+'; 's'; '^'; 'v'; 'p'; '*'; '<'; '>'; 'd'; 'h'; ...
    's']; % Define symbol array
sym = [sym; sym]; % Append symbols
%
rdir = fullfile(div,'Results'); % Define results directory path
bdir = fullfile(rdir,'Bone'); % Define bone directory path
%
ns = size(sd,1); % Get number of subjects
%
% Data Directory and File Name
%
%%
% Get Standard Grid
%
nq = size(tq,1); % Get number of grid points
%
% Posterior and Medial/Lateral Cutoffs
%
%
% Get Logical Indices to the Posterior and Medial/Lateral Regions
%
ipos = tq<=tmin; % Identify posterior points
idl = zq>=y0; % Identify lateral points
idm = zq<y0; % Identify medial points
%
% Get Unique Z Values and Index into Theta Values (in Columns)
%
zu = unique(zq); % Get unique Z values
nu = size(zu,1); % Get number of unique Z values
%
izl = zu>=y0; % Identify lateral Z values
nl = sum(izl); % Count lateral Z values
%
izm = ~izl; % Identify medial Z values
nm = sum(izm); % Count medial Z values
%
nid = reshape(1:nq,nu,nq/nu)'; % Create grid index matrix
%
% Initialize Figure with First Anterior Points (Boundary) for All
% Subjects
%
if iplt % Check if plotting is enabled
    hf = figure; % Create new figure
    orient landscape; % Set figure orientation
    plot(tq,zq,'k.','MarkerSize',6,'LineWidth',0.5); % Plot grid points
    hold on; % Enable hold for multiple plots
    %
    hb = gobjects(ns,1); % Initialize plot handles array
    %
end
%
% Loop through MAT Files
%
idal = zeros(nl,ns); % Initialize lateral anterior indices
idam = zeros(nm,ns); % Initialize medial anterior indices
%
for i = 1:ns % Loop through subjects
%
    fstr = sd(i).FFE.femur.bfnam; % Get femur file name
    fstr = [fstr(1:6) fstr(15:17)]; % Format file name
%
    fnams(i,:) = fstr; % Store file name
    load(fullfile(bdir,[fstr '_femurCS.mat']),'datfb'); % Load femur data
    %
    % Get Trochlea Planes
    %
    [~,~,~,~,ppl,pnl,ppm,pnm] = fem_plan_AD(datfb,tmin,iplt,fstr); % Get trochlea planes
    %
    if iplt % Check if plotting is enabled
        fh = get(0,'Child'); % Get figure handles
        [~,idmx] = max(cell2mat(get(fh,'Number'))); % Find maximum figure number
        delete(fh(idmx)); % Delete posterior-central plane plot
    end
    %
    % Get Trochlea Points for Each Slice
    %
    nsl = size(datfb,1); % Get number of slices
    ip = NaN(nsl,3); % Initialize intersection points
    for ks = 1:nsl % Loop through slices
        ipchk = plsect2(ppl,pnl,datfb{ks}); % Check lateral plane intersection
        if isempty(ipchk)||ipchk(:,2)<y0 % Check if intersection is valid
            ipchk = plsect2(ppm,pnm,datfb{ks}); % Check medial plane intersection
            if isempty(ipchk)||ipchk(:,2)>y0 % Check if intersection is valid
                continue; % Skip invalid intersection
            else
                ip(ks,:) = ipchk; % Store medial intersection
            end
        else
            ip(ks,:) = ipchk; % Store lateral intersection
        end
    end
    %
    if iplt % Check if plotting is enabled
        plot3(ip(:,1),ip(:,2),ip(:,3),'ro'); % Plot intersection points
    end
    %
    % Convert to Polar Coordinates
    %
    [theta,r,z] = cart2pol(ip(:,1),ip(:,3),ip(:,2)); % Convert to polar coordinates
    %
    % Trap for Angles Greater Than 90 degrees (pi/2 radians and Z-Axis)
    %
    id = find(theta>pi/2); % Find angles exceeding 90 degrees
    if ~isempty(id) % Check if angles need adjustment
        theta(id) = theta(id)-2*pi; % Adjust angles
    end
    %
    theta = rad2deg(theta); % Convert angles to degrees
    ipp = [theta r z]; % Combine polar coordinates
    %
    % Scale "Z" Coordinates
    %
    zqs = zsc(i,1)*zq; % Scale Z coordinates
    %
    % Get Unique Scaled Z Values in the Lateral and Medial Compartments
    %
    zus = zqs(nid(1,:)); % Get scaled Z values
    zul = zus(izl); % Get lateral Z values
    zum = zus(izm); % Get medial Z values
    %
    % Use Linear Interpolation to Get Theta Cutoffs
    %
    idxl = ipp(:,3)>=y0; % Identify lateral points
    idxm = ipp(:,3)<y0; % Identify medial points
    tcutl = interp1(ipp(idxl,3),ipp(idxl,1),zul,'linear','extrap'); % Interpolate lateral theta cutoffs
    tcutm = interp1(ipp(idxm,3),ipp(idxm,1),zum,'linear','extrap'); % Interpolate medial theta cutoffs
    %
    % Find Anterior/Trochlea Grid Points on Scaled Grid
    %
    if iplt % Check if plotting is enabled
        figure; % Create new figure
        orient landscape; % Set figure orientation
        plot(tq,zqs,'k.','MarkerSize',6,'LineWidth',0.5); % Plot scaled grid
        hold on; % Enable hold
        plot(ipp(:,1),ipp(:,3),'rs'); % Plot trochlea points
        plot(tcutl,zul,'bo'); % Plot lateral cutoffs
        plot(tcutm,zum,'mo'); % Plot medial cutoffs
    end
    %
    for kz = 1:nu % Loop through Z values
        idx = nid(:,kz); % Get grid indices
        zz = zqs(idx(1)); % Get scaled Z value
        if zz<0 % Check for medial compartment
            im = zz==zum; % Identify medial Z
            idc = tq(idx)>tcutm(im); % Find anterior points
            idc = idx(idc); % Get indices
            idam(im,i) = idc(1); % Store first medial anterior index
        else % Handle lateral compartment
            il = zz==zul; % Identify lateral Z
            idc = tq(idx)>tcutl(il); % Find anterior points
            idc = idx(idc); % Get indices
            idal(il,i) = idc(1); % Store first lateral anterior index
        end
        if iplt % Check if plotting is enabled
            plot(tq(idc),zqs(idc),'go'); % Plot anterior points
            plot(tq(idc(1)),zqs(idc(1)),'bs','MarkerSize',7); % Plot first anterior point
        end
    end
    %
    if iplt % Check if plotting is enabled
        xlabel(['\Theta (' deg ')'],'FontSize',12,'FontWeight','bold'); % Set X-axis label
        ylabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Set Y-axis label
        title({'Anterior Grid Points'; fstr},'FontSize',16, ...
            'FontWeight','bold','interpreter','none'); % Set title
    end
    %
    % Plot Boundary on Figure for All Subjects
    %
    if iplt % Check if plotting is enabled
        figure(hf); % Switch to main figure
        id = idal(:,i); % Get lateral boundary indices
        hb(i) = plot(tq(id),zq(id),['k' sym(i)],'Color',clr(i,:), ...
            'LineWidth',0.5); % Plot lateral boundary
        id = idam(:,i); % Get medial boundary indices
        plot(tq(id),zq(id),['k' sym(i)],'Color',clr(i,:), ...
            'LineWidth',0.5); % Plot medial boundary
    end
    %
end
%
% Get Average First Anterior Points of All Subjects
%
idala = round(mean(idal,2)); % Compute average lateral anterior indices
idama = round(mean(idam,2)); % Compute average medial anterior indices
%
% Get Anterior/Trochlea Grid Points of All Subjects
%
zul = zu(izl); % Get lateral Z values
zum = zu(izm); % Get medial Z values
%
ial = cell(nl,1); % Initialize lateral anterior indices cell
iam = cell(nm,1); % Initialize medial anterior indices cell
%
for kz = 1:nu % Loop through Z values
    idx = nid(:,kz); % Get grid indices
    zz = zq(idx(1)); % Get Z value
    if zz<0 % Check for medial compartment
        im = zz==zum; % Identify medial Z
        idz = idx>=idama(im); % Find anterior points
        iam{im} = nid(idz,kz); % Store medial anterior indices
    else % Handle lateral compartment
        il = zz==zul; % Identify lateral Z
        idz = idx>=idala(il); % Find anterior points
        ial{il} = nid(idz,kz); % Store lateral anterior indices
    end
end
%
idl_ant = false(nq,1); % Initialize lateral anterior logical index
idm_ant = idl_ant; % Initialize medial anterior logical index
idl_ant(cell2mat(ial)) = true; % Set lateral anterior points
idm_ant(cell2mat(iam)) = true; % Set medial anterior points
iant = idl_ant|idm_ant; % Combine anterior points
%
% Plot Average Anterior Grid Points
%
if iplt % Check if plotting is enabled
    %
    figure(hf); % Switch to main figure
    %
    % Plot Anterior/Trochlea Points
    %
    plot(tq(iant),zq(iant),'r.','MarkerSize',6,'LineWidth',0.5); % Plot anterior points
    %
    % Plot Legend, Axis Labels, and Title
    %
    hl = legend(hb,fnams(:,1:5)); % Create legend
    set(hl,'Interpreter','none','Location','East'); % Set legend properties
    xlabel(['\Theta (' deg ')'],'FontSize',12,'FontWeight','bold'); % Set X-axis label
    ylabel('Z (mm)','FontSize',12,'FontWeight','bold'); % Set Y-axis label
    title('Trochlea Grid Points Shown in Red','FontSize',16, ...
        'FontWeight','bold'); % Set title
    if iprt % Check if printing is enabled
        print('-dpsc2','-r300','freg_axf2.ps'); % Print figure
    end
end
%
% Get Remaining Indices to the Six (6) Regions
%
ictr = true(nq,1); % Initialize center region index
ictr(iant) = false; % Exclude anterior points
ictr(ipos) = false; % Exclude posterior points
%
idl_pos = idl&ipos; % Identify lateral posterior region
idl_ctr = idl&ictr; % Identify lateral center region
%
idm_pos = idm&ipos; % Identify medial posterior region
idm_ctr = idm&ictr; % Identify medial center region
%
% Save Indices to Standard Grid MAT File
%
% save fb_rngf.mat -append idl_ant idl_ctr idl_pos idm_ant idm_ctr ...
%     idm_pos;
%
% Plot Regions
%
if iplt % Check if plotting is enabled
    figure; % Create new figure
    orient landscape; % Set figure orientation
    plot(tq(idl_ant),zq(idl_ant),'ro','MarkerSize',4,'LineWidth',0.5); % Plot lateral anterior
    hold on; % Enable hold
    plot(tq(idm_ant),zq(idm_ant),'mo','MarkerSize',4,'LineWidth',0.5); % Plot medial anterior
    plot(tq(idl_ctr),zq(idl_ctr),'gs','MarkerSize',4,'LineWidth',0.5); % Plot lateral center
    plot(tq(idm_ctr),zq(idm_ctr),'bs','MarkerSize',4,'LineWidth',0.5); % Plot medial center
    plot(tq(idl_pos),zq(idl_pos),'y^','MarkerSize',4,'LineWidth',0.5); % Plot lateral posterior
    plot(tq(idm_pos),zq(idm_pos),'c^','MarkerSize',4,'LineWidth',0.5); % Plot medial posterior
    axis tight; % Adjust axis limits
    xlabel({['\Theta (' deg ')'];['\leftarrow posterior / anterior', ...
        ' \rightarrow']},'FontSize',12,'FontWeight','bold'); % Set X-axis label
    ylabel({'Z (mm)';'\leftarrow medial / lateral \rightarrow'}, ...
        'FontSize',12,'FontWeight','bold'); % Set Y-axis label
    title('Average of All Femurs','FontSize',16,'FontWeight','bold'); % Set title
    if iprt % Check if printing is enabled
        print('-dpsc2','-r300','-append','freg_axf2.ps'); % Print figure
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
return % Exit the function