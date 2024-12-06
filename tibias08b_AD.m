%#######################################################################
%
%                      * TIBIAS 08 Bone Program *
%
%          Program to read and plot digitized points from a tibia in the
%     tibia coordinate system.  Plots medial and lateral compartment
%     subchondral bone data for user selected knee.  The tibial
%     coordinate system for the knee is saved to a Matlab MAT file.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  See tibia_cs8.m for more information on the tibia
%             coordinate system.
%
%             3. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi6.m, rotxyz.m, tibia_cs8.m,
%             tri_area.m, tri_fix2.m and tri_norm.m must be in the
%             current directory or path.
%
%             4.  This file outputs a PostScript file, tibias08b.ps,
%             with plots of the tibia coordinate system and the
%             subchondral bone into the directory of the knee data.
%
%             5.  The subchondral bone data in the tibia coordinate
%             system and the coordinate transformation matrix and
%             origin are saved in the directory of the knee data in a
%             Matlab MAT file, kid_tibiaCS.mat, where kid is the knee
%             identifier (xxx_R/L, xxx is the knee number and either R
%             or L for the right or left knee).
%
%             6.  The aspect ratio of the tibia width to the tibia
%             height is out to a CSV file, kid_AspectRatio.csv.
%
%             7.  This is an updated version of tibias07b.m.
%
%     07-May-2019 * Mack Gardner-Morse
%

%#######################################################################
%%
clear;
close all;
clc;
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2;              % Minimum distance between distinct points
iflag = true;           % Print message if duplicates found

% Output PS and MAT File Names
%
psnam_suffix = '_tibias08b.pdf';
mnam_suffix = '_tibiaCS.mat';
%
% Compartment Colors and Labels
%
tclrs = [0 0 0.7; 0 0.5 0];            % Deep blue and dark green
%
tcmpt = ['lateral'
    'medial '];

div = uigetdir;
ddir = fullfile(div);
sd = dir(ddir);         % All subdirectories
snams = {sd.name}';
id = startsWith(snams,{'.'; '..'});    % Current and parent directories
snams = snams(~id);     % Just subject directories
sd=sd(~id);
sd=rmfield(sd,{'bytes','date','datenum','isdir'});
idx = listdlg('ListString',snams,'ListSize',[200 600],'Name', ...
    'Tibias','PromptString', ...
    {'    Select all subdirectories for'; ...
    'calculating cartilage thicknesses.'});
%
if isempty(idx)
    return;
end
sd=sd(idx);
snams = string(snams(idx));    % Make subject names into character array

%
ns = size(snams,1);     % Number of subjects
nss = int2str(ns);
snam_paths=strings(ns,1);
ffe_paths=strings(ns,1);
rho_paths=strings(ns,1);
t2s_paths=strings(ns,1);


for i=1:ns
    snam_paths(i)=fullfile(ddir,snams(i));
    ffe_paths(i)=fullfile(snam_paths(i),'Visit 2','FFE');
    rho_paths(i)=fullfile(snam_paths(i),'Visit 2','RHO');
    t2s_paths(i)=fullfile(snam_paths(i),'Visit 2','T2S');
    sd(i).folder=char(snam_paths(i));

    %Get FFE Tibia Bone Names
    bnams=dir(fullfile(ffe_paths(i),'*_L_SAG_TIB*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).FFE.tibia.bfnam=char(bnams);

    %Get RHO Tibia Bone Names
    bnams=dir(fullfile(rho_paths(i),'Tibia','*_L_SAG_TIB*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).RHO.tibia.bfnam=char(bnams);

    %Get T2S Tibia Bone Names
    bnams=dir(fullfile(t2s_paths(i),'Tibia','*_L_SAG_TIB*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).T2S.tibia.bfnam=char(bnams);

    %Get AX Tibia Files
    ax_list=dir(fullfile(ffe_paths(i),'*_L_AX_TIB*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).FFE.tibia.axfnam=char(ax_list);

    ax_list=dir(fullfile(rho_paths(i),'Tibia','*_L_AX_TIB*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).RHO.tibia.axfnam=char(ax_list);


    ax_list=dir(fullfile(t2s_paths(i),'Tibia','*_L_AX_TIB*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).T2S.tibia.axfnam=char(ax_list);


end
ffe_paths=ffe_paths';
rho_paths=rho_paths';
t2s_paths=t2s_paths';

if isequal(snams,0)
    return;
end
%%


% Setup Tibial Coordinate System Figure
%
% Make Results Folder
%
%dstr = datetime('today');            % Get date as a string
rdir = fullfile(ddir,'Results','Bone');
if ~isfolder(rdir)
    mkdir(rdir);
end
%
for i=1:ns
    for j=1:3
        if j==1
            ax_path=fullfile(ffe_paths(i),sd(i).FFE.tibia.axfnam);
            fstr=sd(i).FFE.tibia.bfnam;
            bone_path=fullfile(ffe_paths(i),sd(i).FFE.tibia.bfnam);
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            ax_path=fullfile(rho_paths(i),'Tibia',sd(i).RHO.tibia.axfnam);
            fstr=sd(i).RHO.tibia.bfnam;
            bone_path=fullfile(rho_paths(i),'Tibia',fstr);
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            ax_path=fullfile(t2s_paths(i),'Tibia',sd(i).T2S.tibia.axfnam);
            fstr=sd(i).T2S.tibia.bfnam;
            bone_path=fullfile(t2s_paths(i),'Tibia',fstr);
            fstr=[fstr(1:6) fstr(15:17)];

        end

        hf1 = figure;
        orient tall;
        view(3);
        hold on;
        %
        % Get Tibial Coordinate System, Proximal Tibia Outline (PTO) and Aspect
        % Ratio from the PTO
        %
        % NOTE:  PTO is in the tibia coordinate system.
        %

        [xyzc,xyzr,aspect,widt,height,xyzpto,prox_tib_ml,prox_tib_ap] = tibia_cs8(ax_path, ...
            false,true);
        %
        % Finish Plot
        %
        title({fstr; 'Tibia Coordinate System'},'FontSize',16, ...
            'FontWeight','bold', 'Interpreter','none');
        axis equal;
        %
        psnam = [fstr psnam_suffix];
        psnam = fullfile(rdir,psnam);
        %
        set(hf1, 'units','normalized','outerposition',[0 0 1 1]);
        exportgraphics(hf1, psnam, "Resolution", 300);
        close(hf1);
        %
        % Output Aspect Ratio to Screen and a CSV File
        %
        % diff_ht = widt-height;  % Width and height should be equal
        % % nchange_sl = round(diff/sl_thick);     % Number of slices to move
        % %
        % % fprintf(1,'\n\n Aspect Ratio = %.3f\n',aspect)
        % % if nchange_sl>0
        % %   fprintf(1,'   Go Distally %.0f Slices\n\n',abs(nchange_sl));
        % % elseif nchange_sl<0
        % %   fprintf(1,'   Go Proximally %.0f Slices\n\n',abs(nchange_sl));
        % % else
        % %   fprintf(1,'   Aspect Ratio is Close!\n\n');
        % % end
        % %
        % fid = fopen(fullfile(pnams(i),[kids(i) '_AspectRatio.csv']),'w');
        % %
        %
        % fprintf(fid,'"kid","AR","Width","Height"\n');
        % fprintf(fid,',,"(mm)","(mm)"\n');
        % fprintf(fid,'%s,%g,%g,%g\n',fnama,aspect,widt,height);
        % %
        % fclose(fid);
        %
        % Get Raw Tibia Bone Data
        %
        roi = rd_roi6(bone_path);
        %
        roinams = upper(char(roi.name));
        idc(2) = find(strcmp('MTCS',cellstr(roinams)));  % Medial compartment
        idc(1) = find(strcmp('LTCS',cellstr(roinams)));  % Lateral compartment
        %
        % Raw Data Figure
        %
        hf2 = figure;
        orient landscape;
        view(3);
        hold on;
        xlabel('X (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
        title({[fstr ' - MRI CS']; 'Blue - Lateral, Green - Medial'}, ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        % Transformed Data Figure
        %
        hf3 = figure;
        orient landscape;
        view(3);
        hold on;
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
        title([fstr ' - Tibial CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %%
        %
        % Loop through Compartments
        %

        for l = 1:2

            %
            % Get Data
            %
            dat1 = roi(idc(l)).data';
            nsl1 = size(dat1,1);              % Number of slices
            %
            % Check for Duplicates, Direction of Digitization and Transform Data to
            % Tibia Coordinate System
            %
            datt = cell(nsl1,1); % Transformed data
            %
            for n = 1:nsl1
                xyz = dat1{n};
                if isempty(xyz)
                    cmprt = deblank(tcmpt(l,:));
                    error([' *** ERROR in tibias08b:  No coordinates for ', ...
                        'slice ' int2str(n) ' in ' cmprt ' compartment!']);
                end
                xyz = fix_pts_AD(xyz,tol,iflag);    % Check for duplicates
                npts = size(xyz,1);
                [~,imx] = max(xyz(:,2));
                [~,imn] = min(xyz(:,2));
                if imn<imx
                    xyz = xyz(1:npts,:);           % Anterior to posterior
                else
                    xyz = xyz(npts:-1:1,:);       % Anterior to posterior
                end
                xyzt = xyz-repmat(xyzc,npts,1); % Center data
                xyzt = xyzt*xyzr;               % Transform data to tibia CS
                datt{n} = xyzt;
                %
                %
                % Plot Raw Data
                %
                figure(hf2);
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
                figure(hf3);
                %
                plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color',tclrs(l,:), ...
                    'MarkerSize',8,'LineWidth',1);
                if xyzt(1,1)<0
                    text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                else
                    text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
                        'Color','k','FontSize',10);
                end
            end

            %
            % Save Data into Compartment Specific Variables
            %

            %
            % Get Surface Triangulation
            %
            trit = mk_tri6(datt);
            xyzt = cell2mat(datt);
            trit = tri_fix2(trit,xyzt); 


            eval(['dat' tcmpt(l,1) ' = datt;']);
            eval(['tri' tcmpt(l,1) ' = trit;']);
            eval(['xyz' tcmpt(l,1) ' = xyzt;']);

            %
            %   clear datt trit xyzt;
            %
        end

        %%
        %
        % Save Plots
        %
        set(hf2, 'units','normalized','outerposition',[0 0 1 1]);
        set(hf3, 'units','normalized','outerposition',[0 0 1 1]);
        exportgraphics(hf2, psnam, "Resolution", 300, 'Append', true);
        exportgraphics(hf3, psnam, "Resolution", 300, 'Append', true);
        close(hf2);
        close(hf3);
        %
        % Plot Transformed Bone Surface Data
        %
        

        hf4 = figure;
        orient landscape;
        %
        trisurf(tril,xyzl(:,1),xyzl(:,2),xyzl(:,3),'FaceColor','interp', ...
            'EdgeColor',tclrs(1,:),'LineWidth',1);
        hold on;
        trisurf(trim,xyzm(:,1),xyzm(:,2),xyzm(:,3),'FaceColor','interp', ...
            'EdgeColor',tclrs(2,:),'LineWidth',1);
        axis equal;
        xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
        ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
        zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
        title([fstr ' - Tibial CS'], ...
            'FontSize',16,'FontWeight','bold','Interpreter','none');
        %
        set(hf4, 'units','normalized','outerposition',[0 0 1 1]);
        exportgraphics(hf4, psnam, "Resolution", 300, 'Append', true);
        close(hf4);
        %
        % Save Data into a Matlab MAT File for Further Processing
        %
        mnam = [fstr mnam_suffix];
        mnam = fullfile(rdir,mnam);
        kid=fstr(1:5);
        ileg=false;
        %
        save(mnam, 'datl','datm','tril','trim','xyzc','xyzr', ...
            'xyzl','xyzm','xyzpto','prox_tib_ml','prox_tib_ap');

        %
    end
end
main_file=fullfile(ddir,'Results', 'Subject Bone Files ONLY.mat');
save (main_file, 'sd');
return 