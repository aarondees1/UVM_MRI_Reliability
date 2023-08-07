%#######################################################################
%
%                * Tibia CARilage Thickness 8 Program *
%
%          M-File which runs through the left or right tibia
%     within selected subject directories and defines a 1 mm by 1 mm
%     grid, scales and projects the grid onto the bone surface mesh,
%     calculates the intersection of the bone normals with the
%     cartilage surface mesh and calculates cartilage thicknesses.  The
%     same grid is used for both the sagittal and coronal cartilage
%     digitizations, so similar points can be averaged across the
%     digitizations.
%
%     NOTES:  1.  The M-files car_thk8.m, gridproj.m, line_fit.m,
%             nod2tri.m, nod_norm.m, pts2lin.m, tri_fix2.m,
%             tri_norm.m, tsect4.m and xprod.m must be in the current
%             path or directory.
%
%             2.  Bone and cartilage MAT files must already exist in
%             the selected subject directories.
%
%             3.  This program produces a MAT file with the scaling for
%             the tibias within the selected directories: tgrid08_*.mat,
%             where * is the number of selected subject directories.
%             Final thicknesses are saved in the MAT file:
%             tcart08_*.mat.  The MAT files are saved in the directory:
%             TibiaCartThk_*_ddmmmyyyy, where * is the number of
%             selected subject directories and ddmmmyyyy is the date in
%             day, month (first three letters) and year format.
%
%             4.  This M-file outputs PostScript files,
%             tcart08_coronal.ps, tcart08_sagittal.ps,
%             tcart08_coronal_thk.ps and tcart08_sagittal_thk.ps with
%             plots of the tibial bony and cartilage surfaces and tibial
%             cartilage thicknesses in the results directory.
%
%             5. Select "Visit 2" for MRIR data directory, then "Rho" or
%             "FFE" for the subdirectory. Make sure to change Flagged areas
%             in the code for a RHO or and FFE scan.
%
%     20-Nov-2019 * Mack Gardner-Morse
%

%#######################################################################
clear;
close all;
clc;
%
% Plot Parameter
%
% iplt = false;           % No plots
iplt = true;            % Plots
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots

csnam = '_tibiaCS.mat';
cnam = '_tibiaCart.mat';
%
% Get All Subject Subdirectories
%
%
%
div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
if ~isfolder(tdir)
    mkdir(tdir);
end
if ~isfolder(gdir)
    mkdir(gdir);
end
load(fullfile(rdir,'Subject Bone and Cartilage Files.mat'));
ns=size(sd,1);

%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1);    % 0 (false) for left/1 (true) for right
kids = repmat(blanks(ns)',1,5);        % Knee IDs
tribl = cell(ns,1);     % Sagittal bone lateral triangular mesh connectivity
tribm = cell(ns,1);     % Sagittal bone medial triangular mesh connectivity
xyzbl = cell(ns,1);     % Sagittal lateral bone point coordinates
xyzbm = cell(ns,1);     % Sagittal medial bone point coordinates
xyz_min_l = zeros(ns,3);   % Minimum individual lateral tibia coordinates
xyz_max_l = zeros(ns,3);   % Maximum individual lateral tibia coordinates
xyz_min_m = zeros(ns,3);   % Minimum individual medial tibia coordinates
xyz_max_m = zeros(ns,3);   % Maximum individual medial tibia coordinates
xyzs = zeros(ns,3);     % Range of proximal tibia outline
%

scan(1,:)='FFE';
scan(2,:)='RHO';
scan(3,:)='T2S';
%% Load Bone Data and Generate min/max and range of Grid
for i = 1:ns
    i_str=int2str(i);

    for j=1:3
        if j==1
            fstr=sd(i).FFE.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end

        bone = load(fullfile(rdir,'Bone',[fstr csnam]));

        if isfield(bone,'xyzpto')

            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.trimesh = bone.tril;']);
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.medial.trimesh = bone.trim;']);
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points = bone.xyzl;']);
            eval(['sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points = bone.xyzm;']);
            eval(['sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline = bone.xyzpto;']);

        else
            error([' *** ERROR in tcart08:  Not able to find proximal ', ...
                'tibial outline.  Please run tibias08b again.']);
        end

        cart = load(fullfile(rdir,'Cartilage',[fstr cnam]));

        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.trimesh = cart.tril;']);
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.medial.trimesh = cart.trim;']);
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.points = cart.xyzl;']);
        eval(['sd(' i_str ').' scan(j,:) '.tibia.cart.medial.points = cart.xyzm;']);

        % if ~exist(cart,'file')
        %     error([' *** ERROR in tcart08:  Not able to find tibia ', ...
        %         'cartilage file:  ',cart,'\n']);
        % end
    end
end
for j=1:3
    for i=1:ns
        i_str=int2str(i);
        % Get Minimum and Maximum Coordinates
        eval(['xyz_min_l(' i_str ',:) = min(sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points);']);     % Lateral Minimums
        eval(['xyz_max_l(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points);']);     % Lateral Maximums

        eval(['xyz_min_m(' i_str ',:) = min(sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points);']);     % Medial Minimums
        eval(['xyz_max_m(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points);']);     % Medial Maximums

        % Get Range of Proximal Tibia Outline
        %
        eval(['xyzs(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline)-min(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline);']);
        %
    end
    % Plot and Calculate Grid
    xyz_min_l=min(xyz_min_l);
    xyz_max_l=max(xyz_max_l);
    xyz_min_m=min(xyz_min_m);
    xyz_max_m=max(xyz_max_m);
    %
    % Generate Scaling (Based on Proximal Tibia Outline) and Uniform Grid
    %

    xyzs_avg = mean(xyzs);
    sc = xyzs./repmat(xyzs_avg,ns,1);
    scx = sc(:,1);          % X scale
    scy = sc(:,2);          % Y scale
    % sc=0;
    % scx=0;
    % scy=0;
    %
    %
    % Get Range and Make a Uniform Rectangular Grid in X and Y for the
    % Lateral Compartment
    %

    x = (floor(xyz_min_l(1)):ceil(xyz_max_l(1)))';    % X range in 1 mm increments
    nxl = size(x,1);
    y = (floor(xyz_min_l(2)):ceil(xyz_max_l(2)))';    % Y range in 1 mm increments
    nyl = size(y,1);
    [X,Y] = meshgrid(x,y);  % Rectangle grid of square quadrilaterals
    xql = X(:);             % Grid points as a vector
    yql = Y(:);             % Grid points as a vector
    nl = size(xql,1);       % Number of grid points in lateral compartment
    %
    quadl = (1:nyl-1)';
    quadl = [quadl quadl+nyl quadl+nyl+1 quadl+1];
    quadl = repmat(quadl,nxl-1,1);
    addcol = repmat(nyl*(0:nxl-2),(nyl-1)*4,1);
    addcol = reshape(addcol,4,(nxl-1)*(nyl-1))';
    quadl = quadl+addcol;
    %
    nel = size(quadl,1);    % Number of lateral quadrilateral elements

    %
    % Get Range and Make a Uniform Rectangular Grid in X and Y for the
    % Medial Compartment
    %

    x = (floor(xyz_min_m(1)):ceil(xyz_max_m(1)))';    % X range in 1 mm increments
    nxm = size(x,1);
    y = (floor(xyz_min_m(2)):ceil(xyz_max_m(2)))';    % Y range in 1 mm increments
    nym = size(y,1);
    [X,Y] = meshgrid(x,y);  % Rectangle grid of square quadrilaterals
    xqm = X(:);             % Grid points as a vector
    yqm = Y(:);             % Grid points as a vector
    nm = size(xqm,1);       % Number of grid points in medial compartment
    %
    quadm = (1:nym-1)';
    quadm = [quadm quadm+nym quadm+nym+1 quadm+1];
    quadm = repmat(quadm,nxm-1,1);
    addcol = repmat(nym*(0:nxm-2),(nym-1)*4,1);
    addcol = reshape(addcol,4,(nxm-1)*(nym-1))';
    quadm = quadm+addcol;
    %
    nem = size(quadm,1);    % Number of medial quadrilateral elements
    %
    % Triangular Grids
    %
    trigl = [quadl(:,1:3) quadl(:,1) quadl(:,3:4)]';
    trigl = reshape(trigl,3,2*nel)';
    %
    trigm = [quadm(:,1:3) quadm(:,1) quadm(:,3:4)]';
    trigm = reshape(trigm,3,2*nem)';
    %
    % Save Analysis Grids
    %
    for i=1:ns
        i_str=int2str(i);
        gnam = [scan(j,:) '_tgrid08_.mat'];
        eval(['sd(' i_str ').' scan(j,:) '.tibia.grid = gnam;']);
        gnam = fullfile(gdir,gnam);   % Grid MAT file
    end
    %
    % if exist(rdir,'dir')
    %   ierr = menu('Overwrite Previous Results','Yes','No')-1;
    %   if ierr
    %     return;
    %   end
    % else
    %   mkdir(rdir);
    % end
    %
    kid=fstr(1:5);
    save(gnam,'nel','nem','nl','nm','nxl','nxm','nyl','nym','kid', ...
        'quadl','quadm','sc','scx','scy','trigl','trigm','xql','xqm', ...
        'yql','yqm','xyz_min_l','xyz_max_l','xyz_min_m','xyz_max_m');


    %% Plot And Calculate Cartilage Thicknesses

    %
    % Initialize Cartilage Arrays
    %

    % Initialize Cartilage Thickness Arrays

    cthkl = zeros(nl,1);       % Lateral sagittal cartilage thicknesses
    xyzil = zeros(nl,3,1);     % Lateral sagittal cartilage intersection
    %
    cthkm = zeros(nm,1);       % Medial sagittal cartilage thicknesses
    xyzim = zeros(nm,3,1);     % Medial sagittal cartilage intersection
    %
    zgl = zeros(nl,1);     % Most superior Z coordinates on bony grid surface
    zgm = zeros(nm,1);     % Most superior Z coordinates on bony grid surface

    for i=1:ns
        i_str=int2str(i);
        eval(['trils = sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.trimesh;']);    % Lateral sagittal triangular mesh
        eval(['trims = sd(' i_str ').' scan(j,:) '.tibia.cart.medial.trimesh;']);    % Medial sagittal triangular mesh
        eval(['xyzls = sd(' i_str ').' scan(j,:) '.tibia.cart.lateral.points;']);    % Lateral sagittal coordinates
        eval(['xyzms = sd(' i_str ').' scan(j,:) '.tibia.cart.medial.points;']);  % Medial sagittal coordinates
        %
        if j==1
            fstr=sd(i).FFE.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.tibia.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end

        % Loop through Subjects (Tibias)
        %
        % Scale Grid
        %
        xgl = xql*scx(i);
        ygl = yql*scy(i);
        xgm = xqm*scx(i);
        ygm = yqm*scy(i);


        %
        % Get Bony Tibia Data
        %
        % Lateral coordinates
        eval(['xyzl = sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.points;']);
        % Medial coordinates
        eval(['xyzm = sd(' i_str ').' scan(j,:) '.tibia.bone.medial.points;']);
        % Lateral mesh
        eval(['tril = sd(' i_str ').' scan(j,:) '.tibia.bone.lateral.trimesh;']);
        % Medial mesh
        eval(['trim = sd(' i_str ').' scan(j,:) '.tibia.bone.medial.trimesh;']);
        %
        % Project Grid onto Bony Surface Mesh
        %
        if j==1
            dist = 2.4;             % Twice slice thickness in FFE images
        elseif j==2
            dist = 6;           % Twice slice thickness in T1rho images
        elseif j==3
            dist = 4;            % Twice slice thickness in T2S images
        end


        zgl(:,1) = gridproj(tril,xyzl,xgl,ygl,1,dist);
        zgm(:,1) = gridproj(trim,xyzm,xgm,ygm,1,dist);
        %
        xyzgl = [xgl ygl zgl(:,1)];
        xyzgm = [xgm ygm zgm(:,1)];
        %
        % Improve the Mesh
        %
        trigli = tri_fix2(trigl,xyzgl);
        trigmi = tri_fix2(trigm,xyzgm);
        %%

        %
        % Plot Cartilage and Bone Surfaces
        %
        if iplt


            % Plot Sagittal Cartilage and Gridded Bone Surfaces
            %
   
            hf2 = figure;
      
            hb2l = trimesh(trigli,xgl,ygl,zgl(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','k');
            hold on;
            hb2m = trimesh(trigmi,xgm,ygm,zgm(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','k');
            %
            hc2l = trimesh(trils,xyzls(:,1),xyzls(:,2), ...
                xyzls(:,3),'LineWidth',0.25, ...
                'FaceColor','none','EdgeColor','b');
            hc2m = trimesh(trims,xyzms(:,1),xyzms(:,2), ...
                xyzms(:,3),'LineWidth',0.25, ...
                'FaceColor','none','EdgeColor','b');
            xlabel('X','FontSize',12,'FontWeight','bold');
            ylabel('Y','FontSize',12,'FontWeight','bold');
            zlabel('Z','FontSize',12,'FontWeight','bold');
            title({fstr; 'Sagittal Cartilage and Bone Meshes'}, ...
                'Interpreter','none','FontSize',16,'FontWeight','bold');
            view(-60,12);
            axis equal;


            psnam = fullfile(tdir,[fstr '_tcart08_sagittal_thk.pdf']);
            orient landscape;
            set(hf2, 'units','normalized','outerposition',[0 0 1 1]);
            exportgraphics(hf2, psnam, "Resolution", 300);
            close(hf2);

            %      pause
        end
        % %

        %
        % Calculate and Save Sagittal Cartilage Thicknesses
        %
        [ct,bi] = car_thk8(trils,xyzls,trigli,xyzgl);   % Lateral
        cthkl(:,1) = ct;    % Save thicknesses
        xyzil(:,:,1) = bi;  % Save cartilage intersection points
        %
        [ct,bi] = car_thk8(trims,xyzms,trigmi,xyzgm);   % Medial
        cthkm(:,1) = ct;    % Save thicknesses
        xyzim(:,:,1) = bi;  % Save cartilage intersection points
        %
        % Plot the Cartilage Thicknesses
        %
        if iplt
            % %

            % Plot Sagittal Cartilage Thicknesses
            %

            hf4 = figure;
        
            hb4l = trimesh(trigli,xgl,ygl,zgl(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','b');
            hold on;
            hb4m = trimesh(trigmi,xgm,ygm,zgm(:,1), ...
                'LineWidth',0.25,'FaceColor','none', ...
                'EdgeColor','b');
            %
            ctl = squeeze(cthkl(:,1));
            idx = find(~isnan(ctl));
            it = nod2tri(idx,trigli,2);
            hs4l = trisurf(trigli(it,:),xgl,ygl,zgl(:,1),ctl, ...
                'FaceColor','interp', ...
                'EdgeColor','b','LineWidth',0.25);
            ctm = squeeze(cthkm(:,1));
            idx = find(~isnan(ctm));
            it = nod2tri(idx,trigmi,2);
            hs4m = trisurf(trigmi(it,:),xgm,ygm,zgm(:,1),ctm, ...
                'FaceColor','interp', ...
                'EdgeColor','b','LineWidth',0.25);
            colormap jet;
            view(-60,12);
            axis equal;
            %
            ct = [ctl; ctm];
            if min(ct)<0
                caxis([min(ct) max(ct)]);
            else
                caxis([0 max(ct)]);
            end
            hcb4 = colorbar;
            xlabel('X','FontSize',12,'FontWeight','bold');
            ylabel('Y','FontSize',12,'FontWeight','bold');
            zlabel('Z','FontSize',12,'FontWeight','bold');
            title({fstr; 'Sagittal Cartilage Thicknesses'}, ...
                'Interpreter','none','FontSize',16,'FontWeight','bold');
            orient landscape;
            set(hf4, 'units','normalized','outerposition',[0 0 1 1]);
            exportgraphics(hf4, psnam, "Resolution", 300, 'Append',true);
            close(hf4);
        end                  % End if iplt
        %

        % %% Save all Figures to Combined PDF
        % h = findobj('type','figure');
        % n = length(h);
        % for i=1:n
        %     exportgraphics(figure(i), 'output.pdf', 'Append', true);
        % end

        %% Save Data and Close Plots

        %
        % Save Tibia Data
        %
        kid=fstr(1:5);
        ileg=false;
        rnam = fullfile(tdir,[fstr '_tcart08_thk.mat']);    % Results MAT file
        save(rnam,'cthkl','cthkm','ileg','kid', ...
            'tribl','tribm','trils','trims', ...
            'xyzbl','xyzbm','xyzil','xyzim', ...
            'xyzls','xyzms','zgl','zgm');
        %
        %close all;              % Close cartilage thickness plots
        %
        % if i ~= ns
        %     clearvars -except k ns bpnams bfnams cpnams cfnams scans conds ilegs kids...
        %         iplt iprt
        % end
    end
end
main_file=fullfile(rdir, 'Subject Files Full.mat');
save (main_file, 'sd');
return