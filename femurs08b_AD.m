%#######################################################################
%
%                      * FEMURS 08 Bone Program *
%
%          Program to read and plot digitized points from a femur in the
%     femoral coordinate system.  Plots lateral and medial condyles and
%     trochlear subchondral bone data for user selected knee.  The
%     femoral coordinate system for the knee is saved to a Matlab MAT
%     file.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  See f_cs_14.m for more information on the femoral
%             coordinate system.
%
%             3.  The angles between slice ends in the triangulation is
%             controlled by the angle input (in degrees) into
%             mk_tri4f.m.  Typical angles are 5-15 degrees.
%
%             4. The Matlab M-files comb_dat.m, cyl_fit.m, cyl_plt.m,
%             dist2cyl.m, f_cs_14.m, fix_pts.m, mk_tri4f.m, plane_fit.m,
%             plt_datsl.m, pt2line.m, pts2lin.m, rd_roi6.m, sl_dir.m,
%             tri_area.m, tri_fix2.m, tri_norm.m and xzpl2pol.m must be
%             in the current directory or path.
%
%             5.  This file outputs a PostScript file, femur08b.ps,
%             with plots of the femoral coordinate system and the
%             subchondral bone into the directory of the knee data.
%
%             6.  The subchondral bone data in the femoral coordinate
%             system and the coordinate transformation matrix and
%             origin are saved in the directory of the knee data in a
%             Matlab MAT file, kid_femurCS.mat, where kid is the knee
%             identifier (xxx_R/L, xxx is the knee number and either R
%             or L for the right or left knee).
%
%     12-Aug-2021 * Mack Gardner-Morse
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
% Output PS and MAT File Names
%
psnam_suffix = '_femur08b.ps';
mnam_suffix = '_femurCS.mat';
%
% Regional Colors and Labels
%
fclrs = [0 0 0.7; 0 0.5 0; 0.7 0 0];   % Deep blue, dark green and red
%
fregs = ['lateral '
         'medial  '
         'trochlea'];
%
% Get Sagittal Bone CSV File Name
%


div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
load(fullfile(rdir,'Subject Files Full.mat'));
ns=size(sd,1);
nss = int2str(ns);
ffe_paths=strings(ns,1);
rho_paths=strings(ns,1);
t2s_paths=strings(ns,1);
%
%
% Check for Axial Plane File
%
for i=1:ns
    ffe_paths(i)=fullfile(sd(i).folder,'Visit2','FFE');
    rho_paths(i)=fullfile(sd(i).folder,'Visit2','RHO');
    t2s_paths(i)=fullfile(sd(i).folder,'Visit2','T2S');
   
    %Get FFE Femur Bone Names
    bnams=dir(fullfile(ffe_paths(i),'*_L_SAG_FEM*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).FFE.femur.bfnam=char(bnams);

    %Get RHO Femur Bone Names
    bnams=dir(fullfile(rho_paths(i),'Femur','*_L_SAG_FEM*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).RHO.femur.bfnam=char(bnams);

    %Get T2S Femur Bone Names
    bnams=dir(fullfile(t2s_paths(i),'Femur','*_L_SAG_FEM*.csv'));
    bnams={bnams.name};
    bnams=bnams(~contains(bnams,'dup','IgnoreCase',true)& ...
        ~contains(bnams,'MGG')&~contains(bnams,'LD'));
    sd(i).T2S.femur.bfnam=char(bnams);

    %Get AX Femur Files
    ax_list=dir(fullfile(ffe_paths(i),'*_L_AX_FEM*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).FFE.femur.axfnam=char(ax_list);

    ax_list=dir(fullfile(rho_paths(i),'Femur','*_L_AX_FEM*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).RHO.femur.axfnam=char(ax_list);


    ax_list=dir(fullfile(t2s_paths(i),'Femur','*_L_AX_FEM*.csv'));
    ax_list={ax_list.name};
    ax_list=ax_list(~contains(ax_list,'LD'));
    sd(i).T2S.femur.axfnam=char(ax_list);


end
ffe_paths=ffe_paths';
rho_paths=rho_paths';
t2s_paths=t2s_paths';


%
%
% if isempty(fnamas)
%   error(' *** ERROR in femurs08b:  Unable to find axial file!');
% end
% %
% fnamas = {fnamas.name}';     % Get axial file names in a cell array
% %
% % Find Femur Axial File
% %
% ifem = contains(upper(fnamas),'FEM');  % Femur axial file
% fnama = char(fnamas(ifem,:));
% %
% if size(fnama,1)~=1
%   error([' *** ERROR in femurs08b:  Unable to find an unique ', ....
%         'axial file!']);
% end
% 
% %
% % Get Leg
% %
% ileg = strcmpi(kid(5),'R');
%


% Setup Femoral Coordinate System Figure
%

if ~isfolder(bdir)
    mkdir(rdir);
end
for i=1:ns
    for j=1:3
        if j==1
            ax_path=fullfile(ffe_paths(i),sd(i).FFE.femur.axfnam);
            fstr=sd(i).FFE.femur.bfnam;
            bone_path=fullfile(ffe_paths(i),sd(i).FFE.femur.bfnam);
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            ax_path=fullfile(rho_paths(i),'Femur',sd(i).RHO.femur.axfnam);
            fstr=sd(i).RHO.femur.bfnam;
            bone_path=fullfile(rho_paths(i),'Femur',fstr);
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            ax_path=fullfile(t2s_paths(i),'Femur',sd(i).T2S.femur.axfnam);
            fstr=sd(i).T2S.femur.bfnam;
            bone_path=fullfile(t2s_paths(i),'Femur',fstr);
            fstr=[fstr(1:6) fstr(15:17)];

        end
hf1 = figure;
orient landscape;
view(3);
hold on;
%
% Get MRI Sagittal Femur Bone Data
%
roi = rd_roi6(bone_path);
%
roinams = {roi.name}';
%roinams = upper(char(roi.name));
ids(3) = find(startsWith(roinams,'TRO'));   % Trochlea
ids(2) = find(startsWith(roinams,'MED'));   % Medial condyle
ids(1) = find(startsWith(roinams,'LAT'));   % Lateral condyle
%
datlb = roi(ids(1)).data';
datmb = roi(ids(2)).data';
dattb = roi(ids(3)).data';
%
% Get MRI Axial Femoral Shaft Outline (FS) and Center Point (CP) Data
%
roia = rd_roi6(ax_path);
%
roinama = {roia.name}';
ida(2) = find(startsWith(roinama,'PROX'));  % Proximal femoral shaft 
ida(1) = find(startsWith(roinama,'CENT'));  % Center point
%
datcp = roia(ida(1)).data;
datfs = roia(ida(2)).data';
%
% Get Femoral Coordinate System
%
[xyzc,xyzr,xyzl,r] = f_cs_14(datlb,datmb,datfs,datcp,true,hf1);
%
% Finish Plot
%
title({fstr; 'Femoral Coordinate System'},'FontSize',16, ...
      'FontWeight','bold', 'Interpreter','none');
axis equal;
%
psnam = [fstr psnam_suffix];
psnam = fullfile(bdir,psnam);
%
print('-dpsc2','-r300','-bestfit',psnam);
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
title({[fstr ' - MRI CS']; ['Blue - Lateral, Green - Medial, ', ...
      'Red - Trochlea']},'FontSize',16,'FontWeight','bold', ...
      'Interpreter','none');
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
title({[fstr ' - Femoral CS']; ['Blue - Lateral, Green - Medial, ', ...
      'Red - Trochlea']},'FontSize',16,'FontWeight','bold', ...
      'Interpreter','none');
%
% Loop through Regions
%
for l = 1:3
%
% Get Data
%
   eval(['dat1 = dat' fregs(l,1) 'b;']);
   nsl1 = size(dat1,1);              % Number of slices
%
% Check for Duplicates, Direction of Digitization and Transform Data to
% Femoral Coordinate System
%
   datt = cell(nsl1,1); % Transformed data
%
   for n = 1:nsl1
      xyz = dat1{n};
      if isempty(xyz)
        reg = deblank(fregs(l,:));
        error([' *** ERROR in femurs08b:  No coordinates for ', ...
               'slice ' int2str(n) ' in ' reg ' region!']);
      end
      xyz = fix_pts(xyz,tol,iflag);    % Check for duplicates
      npts = size(xyz,1);
      [~,imx] = max(xyz(:,2));
      [~,imn] = min(xyz(:,2));
      if imn<imx
        xyz = xyz(1:npts,:);           % Anterior to posterior
      else
         xyz = xyz(npts:-1:1,:);       % Anterior to posterior
      end
      xyzt = xyz-repmat(xyzc,npts,1); % Center data
      xyzt = xyzt*xyzr;               % Transform data to femur CS
      datt{n} = xyzt;
% 
% Plot Raw Data
%
      figure(hf2);
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
      figure(hf3);
%
      plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color',fclrs(l,:), ...
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
% Save Data into Region Specific Variables
%
   eval(['dat' fregs(l,1) 'bt = datt;']);
%
   clear datt;
%
end
%
% Save Plots
%
print(hf2,'-dpsc2','-r300','-bestfit','-append',psnam);
print(hf3,'-dpsc2','-r300','-bestfit','-append',psnam);
%
% Combine Regions
%
datfb = comb_dat(datlbt,datmbt,dattbt);
%
% Get Surface Triangulation
%
trifb = mk_tri4f(datfb,15);  % Adjust angle (in degrees) to smooth mesh
xyzfb = cell2mat(datfb);
trifb = tri_fix2(trifb,xyzfb);
%
% Plot Transformed Bone Surface Data
% 
hf4 = figure;
orient landscape;
%
trisurf(trifb,xyzfb(:,1),xyzfb(:,2),xyzfb(:,3),'FaceColor','interp', ...
        'EdgeColor',fclrs(1,:),'LineWidth',1);
axis equal;
xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
title([fstr ' - Femoral CS'], ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
print('-dpsc2','-r300','-image','-bestfit','-append',psnam);
%
%Save Data into a Matlab MAT File for Further Processing
%
        mnam = [fstr mnam_suffix];
        mnam = fullfile(bdir,mnam);
        kid=fstr(1:5);
        ileg=false;
%
save(mnam,'datlb','datlbt','datmb','datmbt','dattb','dattbt','datfb','kid', ...
          'ileg','r','trifb','xyzc','xyzfb','xyzl','xyzr');

    end
end
main_file=fullfile(ddir,'Results', 'Full Tibia Femur Bones ONLY.mat');
save (main_file, 'sd');
return