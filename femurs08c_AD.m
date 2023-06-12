%#######################################################################
%
%                    * FEMURS 08 Cartilage Program *
%
%          Program to read, transform to the bony femur coordinate
%     system and plot digitized cartilage points from a femur.  Plots
%     medial and lateral compartment cartilage data for user selected
%     knee.  Plots both coronal and sagittal data.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi6.m, rotxyz.m, tri_area.m, tri_fix2.m
%             and tri_norm.m must be in the current directory or path.
%
%             3.  The femur coordinate system with the coordinate
%             transformation matrix and origin are read from the
%             directory of the knee data from the Matlab MAT file,
%             kid_femurCS.mat, where kid is the knee identifier
%             (xxx_R/L, xxx is the knee number and either R or L
%             for the right or left knee).
%
%             4.  The angles between slice ends in the triangulation is
%             controlled by the angle input (in degrees) into
%             mk_tri4f.m.  Typical angles are 5-15 degrees.
%
%             5.  This M-file outputs a PostScript file,
%             kid_femurs08c.ps, with plots of the femoral cartilage
%             into the directory of the knee data.
%
%             6.  The transformed data and mesh are saved in a Matlab
%             MAT file kid_femurCart.mat.
%
%     17-Aug-2021 * Mack Gardner-Morse
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
% Femur Coordinate and Output PS and MAT File Names
%
csnam = '_femurCS.mat';                % Femur coordinate system file
psnam_suffix = '_femurs08c.ps';               % Output PS file
mnam_suffix = '_femurCart.mat';               % Output MAT file
%

%
% Regional Colors and Labels
%
fclrs = [0 0 0.7; 0 0.5 0; 0.7 0 0];   % Deep blue, dark green and red
%
fregs = ['lateral '
    'medial  '
    'trochlea'];
%
% Get Sagittal Cartilage CSV File Name
%
div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');
cdir = fullfile(rdir, 'Cartilage');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');

load(fullfile(rdir,'Full Tibia Femur Bones ONLY.mat'));
ns=size(sd,1);
ffe_paths=strings(ns,1);
rho_paths=strings(ns,1);
t2s_paths=strings(ns,1);
for i=1:ns
    ffe_paths(i)=fullfile(sd(i).folder,'Visit2','FFE');
    rho_paths(i)=fullfile(sd(i).folder,'Visit2','RHO');
    t2s_paths(i)=fullfile(sd(i).folder,'Visit2','T2S');

    %Get FFE Femur Bone Names
    cnams=dir(fullfile(ffe_paths(i),'*_L_SAGAR_FEM*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).FFE.femur.cfnam=char(cnams);

    %Get RHO Femur Bone Names
    cnams=dir(fullfile(rho_paths(i),'Femur','*_L_SAGAR_FEM*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).RHO.femur.cfnam=char(cnams);

    %Get T2S Femur Bone Names
    cnams=dir(fullfile(t2s_paths(i),'Femur','*_L_SAGAR_FEM*.csv'));
    cnams={cnams.name};
    cnams=cnams(~contains(cnams,'dup','IgnoreCase',true)& ...
        ~contains(cnams,'MGG')&~contains(cnams,'LD'));
    idx = contains(cnams,'_RO');        % Check for _RO files
    if any(idx)
        idc = contains(cnams,'SAGAR');
        idx = idx|~idc;
        cnams = cnams(idx);
    end
    sd(i).T2S.femur.cfnam=char(cnams);

end
ffe_paths=ffe_paths';
rho_paths=rho_paths';
t2s_paths=t2s_paths';
clear idc;

if ~isfolder(cdir)
    mkdir(cdir);
end

for i=1:ns
    for j=1:3
        if j==1
            fstr=sd(i).FFE.femur.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(bdir,[fstr csnam]);
            cart_path=fullfile(ffe_paths(i),sd(i).FFE.femur.cfnam);

        elseif j==2
            fstr=sd(i).RHO.femur.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(bdir,[fstr csnam]);
            cart_path=fullfile(rho_paths(i),'Femur',sd(i).RHO.femur.cfnam);

        elseif j==3
            fstr=sd(i).T2S.femur.cfnam;
            fstr=[fstr(1:6) fstr(17:19)];
            coord_path=fullfile(bdir,[fstr csnam]);
            cart_path=fullfile(t2s_paths(i),'Femur',sd(i).T2S.femur.cfnam);

        end
        % Get Femoral Coordinate System
        fcs = load(fullfile(bdir, [fstr csnam]),'xyzc','xyzr');
        xyzc = fcs.xyzc;        % Origin
        xyzr = fcs.xyzr;        % Rotation matrix
        %
        % Full PostScript File Name
        %
        psnam = [fstr psnam_suffix];
        psnam = fullfile(cdir,psnam);
        %
        % Get Raw Femur Sagittal Cartilage Data
        %
        roi = rd_roi6(cart_path);
        roinams = {roi.name}';
        ids(3) = find(startsWith(roinams,'TRO'));   % Trochlea
        ids(2) = find(startsWith(roinams,'MED'));   % Medial condyle
        ids(1) = find(startsWith(roinams,'LAT'));   % Lateral condyle
        %
        datlc = roi(ids(1)).data';
        datmc = roi(ids(2)).data';
        dattc = roi(ids(3)).data';
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
        title({[fstr ' - MRI CS']; ['Blue - Lateral, Green - Medial']}, ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
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
        title([fstr ' - Femoral CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        % Loop through Regions
        % (Lateral [l==1], Medial [l==2] and Trochlea [l==3])
        %
        for l = 1:3             % Lateral = 1, medial = 2 and trochlea = 3
            %
            % Get Data
            %
            eval(['dat1 = dat' fregs(l,1) 'c;']);
            nsl1 = size(dat1,1);                % Number of slices
            %
            % Check for Duplicates, Direction of Digitization and Transform Data to
            % Femur Coordinate System
            %
            datt = cell(nsl1,1);                % Transformed data
            %
            for n = 1:nsl1
                xyz = dat1{n};
                if isempty(xyz)
                    reg = deblank(fregs(l,:));
                    error([' *** ERROR in femurs08c:  No coordinates for ', ...
                        'slice ' int2str(n) ' in ' reg ' region!']);
                end
                xyz = fix_pts_AD(xyz,tol,iflag,fstr);    % Check for duplicates
                npts = size(xyz,1);
                [~,imx] = max(xyz(:,2));
                [~,imn] = min(xyz(:,2));
                if imn<imx
                    xyz = xyz(1:npts,:);           % Anterior to posterior
                else
                    xyz = xyz(npts:-1:1,:);       % Anterior to posterior
                end
                xyzt = xyz-repmat(xyzc,npts,1);  % Center data
                xyzt = xyzt*xyzr;                % Transform data to femur CS
                datt{n} = xyzt;
                %
                % Plot Raw Data
                %
                figure(hf1);
                %
                plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','Color',fclrs(l,:), ...
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
                    fclrs(l,:),'MarkerSize',8,'LineWidth',1);
                if xyzt(1,1)<0
                    text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                else
                    text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                end
                %
            end
            %
            % Save Data into Compartment Specific Variables
            %
            eval(['dat' fregs(l,1) 'ct = datt;']);
            %
            clear datt;
            %
        end
        %
        % Save Plots
        %
        print(hf1,'-dpsc2','-r300','-bestfit','-append',psnam);
        print(hf2,'-dpsc2','-r300','-bestfit','-append',psnam);
        %
        % Combine Regions
        %
        datfc = comb_dat(datlct,datmct,dattct);
        %
        % Get Surface Triangulation
        %
        trifc = mk_tri4f(datfc,15);  % Adjust angle (in degrees) to smooth mesh
        xyzfc = cell2mat(datfc);
        trifc = tri_fix2(trifc,xyzfc);
        %
        % Plot Transformed Cartilage Surface Data
        %
        hf3 = figure;
        orient landscape;
        %
        trisurf(trifc,xyzfc(:,1),xyzfc(:,2),xyzfc(:,3),'FaceColor', ...
            'interp','EdgeColor',fclrs(1,:),'LineWidth',1);
        axis equal;
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
        title([fstr ' - Femoral CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        print('-dpsc2','-r300','-image','-bestfit','-append',psnam);
        %
        % Save Data into a Matlab MAT File for Further Processing
        %
        mnam = [fstr mnam_suffix];
        mnam = fullfile(cdir,mnam);
        kid=fstr(1:5);
        ileg=false;
        %
        save(mnam,'datlc','datlct','datmc','datmct','dattc','dattct','kid', ...
            'ileg','trifc','xyzfc');
        %
    end
end
main_file=fullfile(rdir, 'Full Tibia Full Femur Subject Files.mat');
save (main_file, 'sd');
return