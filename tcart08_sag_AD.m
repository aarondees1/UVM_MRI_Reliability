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
%
% Plot Parameter
%
% iplt = false;           % No plots
iplt = true;            % Plots
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots
%%
mnam = '_tibiaCS.mat';
%
% Get All Subject Subdirectories
%
% 
% 
div = uigetdir;
ddir = fullfile(div); 
rdir = fullfile(div,'Results');
load(fullfile(rdir,'Subject Bone and Cartillage Files.mat'));
ns=size(sd,1);

%
% Initialize Coordinate and Bone Arrays
%
%ileg = false(ns,1);    % 0 (false) for left/1 (true) for right
kid = repmat(blanks(ns)',1,5);        % Knee IDs
tribl = cell(ns,1);     % Sagittal bone lateral triangular mesh connectivity
tribm = cell(ns,1);     % Sagittal bone medial triangular mesh connectivity
xyzbl = cell(ns,1);     % Sagittal lateral bone point coordinates
xyzbm = cell(ns,1);     % Sagittal medial bone point coordinates
xyzmnl = zeros(ns,3);   % Minimum individual lateral tibia coordinates
xyzmxl = zeros(ns,3);   % Maximum individual lateral tibia coordinates
xyzmnm = zeros(ns,3);   % Minimum individual medial tibia coordinates
xyzmxm = zeros(ns,3);   % Maximum individual medial tibia coordinates
xyzs = zeros(ns,3);     % Range of proximal tibia outline
%

for k = 1:ns
%% Load Bone Data and Generate min/max and range of Grid

%
% Loop through Subjects (Tibias)
%

switch scans{k}
    case 'RHO'
        dist = 3;             % Twice slice thickness in T1rho images
    case 'FFE'
        dist = 2.4;           % Twice slice thickness in FFE images
    case 'T2S'
        dist = 4;
end



% Get Subject Number
%
   %snam = snams(k,:);   % Subject directory

%
% Get Bone Data
%
   bone = load(bfnams{k});
   if isfield(bone,'xyzpto')
     ilegs(k) = bone.ileg;
     kids{k} = bone.kid;
     tribl{k} = bone.tril;
     tribm{k} = bone.trim;
     xyzbl{k} = bone.xyzl;
     xyzbm{k} = bone.xyzm;
     xyzpt = bone.xyzpto;
   else
     error([' *** ERROR in tcart08:  Not able to find proximal ', ...
            'tibial outline.  Please run tibias08b again.']);
   end


% Initialize Overall Minimums and Maximums
%
xyzmnl = zeros(1,3);     % Minimums
xyzmxl = zeros(1,3);     % Maximums

xyzmnm = zeros(1,3);     % Minimums
xyzmxm = zeros(1,3);     % Maximums

for i=1:ns

    fstr=sd(i).FFE_bone;
    fstr=[fstr(1:6) fstr(15:17)];
    load(fullfile(rdir,'Bone',[fstr mnam]));
    %
    % Get Minimum and Maximum Coordinates
%
   xyzmnl(k,:) = min(xyzbl{k});
   xyzmxl(k,:) = max(xyzbl{k});
%
   xyzmnm(k,:) = min(xyzbm{k});
   xyzmxm(k,:) = max(xyzbm{k});

   itstl = tzrn<tzrmn;
   tzrmn(itstl) = tzrn(itstl);
   itstl = tzrx>tzrmx;
   tzrmx(itstl) = tzrx(itstl);
%
% Get Range of Proximal Tibia Outline
%
   xyzs(k,:) = max(xyzpt)-min(xyzpt);
%
end
% Plot and Calculate Grid

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

x = (floor(xyzmnl(1)):ceil(xyzmxl(1)))';    % X range in 1 mm increments
nxl = size(x,1);
y = (floor(xyzmnl(2)):ceil(xyzmxl(2)))';    % Y range in 1 mm increments
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

x = (floor(xyzmnm(1)):ceil(xyzmxm(1)))';    % X range in 1 mm increments
nxm = size(x,1);
y = (floor(xyzmnm(2)):ceil(xyzmxm(2)))';    % Y range in 1 mm increments
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
gnam = fullfile(bpnams{k},[kids{k} '_' scans{k} '_tgrid08_.mat']);   % Grid MAT file
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

save(gnam,'nel','nem','nl','nm','nxl','nxm','nyl','nym','kid', ...
     'quadl','quadm','sc','scx','scy','trigl','trigm','xql','xqm', ...
     'yql','yqm','xyzmnl','xyzmxl','xyzmnm','xyzmxm');


i=1;
cfnams=cell(size(bfnams,1));
cpnams=cell(size(bpnams,1));
while i<=size(bfnams,1)
    cfnams{i}=[bfnams{i}(1:13),'tibiaCart.mat'];
        i=i+1;   
end

end
%% Plot And Calculate Cartilage Thicknesses

%
% Initialize Cartilage Arrays
%

trils = cell(1,1);    % Lateral sagittal triangular mesh
trims = cell(1,1);    % Medial sagittal triangular mesh
xyzls = cell(1,1);    % Lateral sagittal coordinates
xyzms = cell(1,1);    % Medial sagittal coordinates
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
%
% Loop through Subjects (Tibias)
%

%
% Get Subject Number
%
   %snam = snams(k,:);   % Subject directory
  % kid = kids(k,:);     % Knee ID
%
% Get Cartilage Data
%
%   cfnam = fullfile(ddir,snam,[kid '_tibiaCart.mat']);     % Cartilage MAT file name for FFE
   %cfnam = fullfile(char(pnams{k}),char(cfnams{k}));     % Cartilage MAT file name For Rho
   if ~exist(cfnams{k},'file')
     error([' *** ERROR in tcart08:  Not able to find tibia ', ...
            'cartilage file:  ',cfnam,'\n']);
   end
   cart = load(cfnams{k});
%
% Get and Save Cartilage Data
%

   
   trils{1} = cart.tril;
 
   trims{1} = cart.trim;

   xyzls{1} = cart.xyzl;

   xyzms{1} = cart.xyzm;
%
% Reverse X-Axis for Right Knees
%
   if strcmpi(ilegs(k),'R')
     
     xyzls{1}(:,1) = -xyzls{1}(:,1);
 
     xyzms{1}(:,1) = -xyzms{1}(:,1);
   end
%
% Scale Grid
%
%    xgl = xql*scx(k);
%    ygl = yql*scy(k);
%    xgm = xqm*scx(k);
%    ygm = yqm*scy(k);
   xgl = xql;           % No scaling
   ygl = yql;           % No scaling
   xgm = xqm;           % No scaling
   ygm = yqm;           % No scaling

%
% Get Bony Tibia Data
%
   xyzl = xyzbl{1};     % Lateral coordinates
   xyzm = xyzbm{1};     % Medial coordinates
%
   tril = tribl{1};     % Lateral triangular mesh connectivity
   trim = tribm{1};     % Medial triangular mesh connectivity
%
% Project Grid onto Bony Surface Mesh
 if strcmpi(scans{k},'RHO')
   dist = 6;             % Twice slice thickness in T1rho images
 else
   dist = 2.4;           % Twice slice thickness in FFE images
 end
%
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
     if exist('hf2','var')
       figure(hf2);
       clf;
     else
       hf2 = figure;
     end
     hb2l = trimesh(trigli,xgl,ygl,zgl(:,1), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','k');
     hold on;
     hb2m = trimesh(trigmi,xgm,ygm,zgm(:,1), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','k');
%
     hc2l = trimesh(trils{1},xyzls{1}(:,1),xyzls{1}(:,2), ...
                    xyzls{1}(:,3),'LineWidth',0.25, ...
                    'FaceColor','none','EdgeColor','b');
     hc2m = trimesh(trims{1},xyzms{1}(:,1),xyzms{1}(:,2), ...
                    xyzms{1}(:,3),'LineWidth',0.25, ...
                    'FaceColor','none','EdgeColor','b');
     xlabel('X','FontSize',12,'FontWeight','bold');
     ylabel('Y','FontSize',12,'FontWeight','bold');
     zlabel('Z','FontSize',12,'FontWeight','bold');
%      title({kid; 'Sagittal Cartilage and Bone Meshes'}, ...
%            'Interpreter','none','FontSize',16,'FontWeight','bold');
     view(-60,12);
     axis equal;
     if iprt
       psnam = fullfile([bfnams{k}(1:13) 'tcart08_sagittal.ps']);
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
%      pause
   end
% %

%
% Calculate and Save Sagittal Cartilage Thicknesses
%
   [ct,bi] = car_thk8(trils{1},xyzls{1},trigli,xyzgl);   % Lateral
   cthkl(:,1) = ct;    % Save thicknesses
   xyzil(:,:,1) = bi;  % Save cartilage intersection points
%
   [ct,bi] = car_thk8(trims{1},xyzms{1},trigmi,xyzgm);   % Medial
   cthkm(:,1) = ct;    % Save thicknesses
   xyzim(:,:,1) = bi;  % Save cartilage intersection points
%
% Plot the Cartilage Thicknesses
%
   if iplt
% %

% Plot Sagittal Cartilage Thicknesses
%
     if exist('hf4','var')
       figure(hf4);
       clf;
     else
       hf4 = figure;
     end
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
%      title({kid; 'Sagittal Cartilage Thicknesses'}, ...
%     'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       psnam = fullfile([bfnams{k}(1:13) 'tcart08_sagittal_thk.ps']);
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
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
scan=scans{k};
kid=kids{k};
cond=conds{k};
ileg=ilegs(k);
rnam = fullfile([bfnams{k}(1:13) 'tcart08_thk.mat']);    % Results MAT file
save(rnam,'cthkl','cthkm','ileg','kid', 'scan', 'cond', ...
     'tribl','tribm','trils','trims', ...
     'xyzbl','xyzbm','xyzil','xyzim', ...
     'xyzls','xyzms','zgl','zgm','xyzpt');
%
%close all;              % Close cartilage thickness plots
%
if k ~= ns
clearvars -except k ns bpnams bfnams cpnams cfnams scans conds ilegs kids...
    iplt iprt
end

return