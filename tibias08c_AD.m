%#######################################################################
%
%                    * TIBIAS 08 Cartilage Program *
%
%          Program to read, transform to the bony tibia coordinate
%     system and plot digitized cartilage points from a tibia.  Plots
%     medial and lateral compartment cartilage data for user selected
%     knee.  Plots both coronal and sagittal data.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi4.m, rotxyz.m, tri_area.m, tri_fix2.m
%             and tri_norm.m must be in the current directory or path.
%
%             3.  The tibia coordinate system with the coordinate
%             transformation matrix and origin are read from the
%             directory of the knee data from the Matlab MAT file,
%             kid_tibiaCS.mat, where kid is the knee identifier
%             (xxx_R/L, xxx is the knee number and either R or L
%             for the right or left knee).
%
%             4.  This M-file outputs a PostScript file,
%             kid_tibias08c.ps, with plots of the tibial cartilage into
%             the directory of the knee data.
%
%             5.  The transformed data and mesh are saved in a Matlab
%             MAT file kid_tibiaCart.mat.
%
%     27-Sep-2019 * Mack Gardner-Morse
%

%#######################################################################
clear;
close all;
clc;
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2;              % Minimum distance between distinct points
iflag = true;           % Print message if duplicates found
%
% Tibial Coordinate and Output PS and MAT File Names
%
csnam = '_tibiaCS.mat';                % Tibia coordinate system file
psnam_suffix = '_tibias08c.ps';               % Output PS file
mnam_suffix = '_tibiaCart.mat';               % Output MAT file
%

%
% Compartment Colors and Labels
%
tclrs = [0 0 0.7; 0 0.5 0];            % Deep blue and dark green
%
tcmpt = ['lateral'
    'medial '];
%
% Get Sagittal Cartilage CSV File Name
%
div = uigetdir;
rdir = fullfile(div,'Results');

load(fullfile(rdir,'Subject Bone Files ONLY.mat'));
ns=size(sd,1);
ffe_paths=strings(ns,1);
rho_paths=strings(ns,1);
t2s_paths=strings(ns,1);
for i=1:ns
    ffe_paths(i)=fullfile(sd(i).folder,'Visit2','FFE');
    rho_paths(i)=fullfile(sd(i).folder,'Visit2','RHO');
    t2s_paths(i)=fullfile(sd(i).folder,'Visit2','T2S');

    %Get FFE Tibia Bone Names
    cnams=dir(fullfile(ffe_paths(i),'*_L_SAGAR_TIB*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).FFE.tibia.cfnam=char(cnams);

    %Get RHO Tibia Bone Names
    cnams=dir(fullfile(rho_paths(i),'Tibia','*_L_SAGAR_TIB*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).RHO.tibia.cfnam=char(cnams);

    %Get T2S Tibia Bone Names
    cnams=dir(fullfile(t2s_paths(i),'Tibia','*_L_SAGAR_TIB*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).T2S.tibia.cfnam=char(cnams);

end
ffe_paths=ffe_paths';
rho_paths=rho_paths';
t2s_paths=t2s_paths';
clear idc;

%%
cdir = fullfile(rdir,'Cartilage');
if ~isfolder(cdir)
    mkdir(cdir);
end

for i=1:ns
    for j=1:3

        % Get Tibial Coordinate System
        if j==1
            fstr=sd(i).FFE.tibia.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(rdir,'Bone',[fstr csnam]);
            cart_path=fullfile(ffe_paths(i),sd(i).FFE.tibia.cfnam);

        elseif j==2
            fstr=sd(i).RHO.tibia.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(rdir,'Bone',[fstr csnam]);
            cart_path=fullfile(rho_paths(i),'Tibia',sd(i).RHO.tibia.cfnam);
        
        elseif j==3
            fstr=sd(i).T2S.tibia.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(rdir,'Bone',[fstr csnam]);
            cart_path=fullfile(t2s_paths(i),'Tibia',sd(i).T2S.tibia.cfnam);

        end
        tcs = load(fullfile(rdir,'Bone',[fstr csnam]),'xyzc','xyzr');
        xyzc = tcs.xyzc;        % Origin
        xyzr = tcs.xyzr;        % Rotation matrix
        %
        % Full PS File Name
        %
        psnam = [fstr psnam_suffix];
        psnam = fullfile(cdir,psnam);
        %
        % Loop through Digitizations (Coronal [k==1] and Sagittal [k==2])
        %

        % Get Raw Tibia Cartilage Data

        roi = rd_roi6(cart_path);
        roinams = upper(char(roi.name));
        idc(2) = find(strcmp('MACS',cellstr(roinams)));  % Medial compartment
        idc(1) = find(strcmp('LACS',cellstr(roinams)));  % Lateral compartment
        dtxt = 'Sagittal Digitization';
        %
        % Raw Data Figure
        %
        hf1 = figure;
        orient landscape;
        view(3);
        hold on;
        xlabel('X (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
        title({[fstr ' - MRI CS FIXEDPTS']; dtxt; ['Blue - Lateral, Green - ', ...
            'Medial']},'FontSize',16,'FontWeight','bold', ...
            'Interpreter','none');
        %
        % Transformed Data Figure
        %
        hf2 = figure;
        orient landscape;
        view(3);
        hold on;
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
        title({[fstr ' - Tibial CS FIXEDPTS']; dtxt}, ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        % Loop through Compartments (Lateral [l==1] and Medial [l==2])
        %
        for l = 1:2          % Lateral = 1 and medial = 2
            %
            % Get Data
            %
            dat1 = roi(idc(l)).data';
            nsl1 = size(dat1,1);             % Number of slices
            %
            % Check for Duplicates, Direction of Digitization and Transform Data to
            % Tibia Coordinate System
            %
            datt = cell(nsl1,1);             % Transformed data
            %
            for n = 1:nsl1
                xyz = dat1{n};
                if isempty(xyz)
                    cmprt = deblank(tcmpt(l,:));
                    error([' *** ERROR in tibias08c:  No coordinates for ', ...
                        dtxt ' slice ' int2str(n) ' in ' cmprt, ...
                        'compartment!']);
                end
                xyz = fix_pts_AD(xyz,tol,iflag,fstr); % Check for duplicates
                npts = size(xyz,1);
                [~,imx] = max(xyz(:,2));
                [~,imn] = min(xyz(:,2));
                if imn<imx
                    xyz = xyz(1:npts,:);        % Anterior to posterior
                else
                    xyz = xyz(npts:-1:1,:);    % Anterior to posterior
                end
                xyzt = xyz-repmat(xyzc,npts,1);    % Center data
                xyzt = xyzt*xyzr;                  % Transform data to tibia CS
                datt{n} = xyzt;
                %
                % Plot Raw Data
                %
                figure(hf1);
                %
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','Color',tclrs(l,:), ...
                    'MarkerSize',8,'LineWidth',1);
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
                figure(hf2);
                %
                plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color', ...
                    tclrs(l,:),'MarkerSize',8,'LineWidth',1);
                if xyzt(1,1)<0
                    text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                else
                    text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                end
            end
            %
            % Get Surface Triangulation
            %
            trit = mk_tri6(datt);
            xyzt = cell2mat(datt);
            trit = tri_fix2(trit,xyzt);
            %
            % Save Data into Compartment Specific Variables
            %
            eval(['dat' tcmpt(l,1) ' = datt;']);
            eval(['tri' tcmpt(l,1) ' = trit;']);
            eval(['xyz' tcmpt(l,1) ' = xyzt;']);

            clear datt trit xyzt;
        end

        %
        % Save Plots
        %
        print(hf1,'-dpsc2','-r300','-bestfit','-append',psnam);
        print(hf2,'-dpsc2','-r300','-bestfit','-append',psnam);
        %
        % Plot Transformed Cartilage Surface Data
        %
        hf3 = figure;
        orient landscape;
        %

        trisurf(tril,xyzl(:,1),xyzl(:,2),xyzl(:,3),'FaceColor', ...
            'interp','EdgeColor',tclrs(1,:),'LineWidth',1);
        hold on;
        trisurf(trim,xyzm(:,1),xyzm(:,2),xyzm(:,3),'FaceColor', ...
            'interp','EdgeColor',tclrs(2,:),'LineWidth',1);

        axis equal;
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
        title({[fstr ' - Tibial CS FIXEDPTS']; dtxt}, ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        print('-dpsc2','-r300','-image','-bestfit','-append',psnam);
        %
       % close([hf1 hf2]);
        %
        % Save Data into a Matlab MAT File for Further Processing
        %
        mnam = [fstr mnam_suffix];
        mnam = fullfile(cdir,mnam);
        kid=fstr(1:5);
        ileg=false;
        %
        save(mnam,'datl',"mnam",'datm','kid','ileg', ...
            'tril','trim','xyzl', 'xyzm');
        %
    end
end
main_file=fullfile(rdir, 'Subject Bone And Cartilage Files.mat');
save (main_file, 'sd');
return 