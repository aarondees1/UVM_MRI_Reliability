clear;
close all;
clc;
iprt =true;

csnam = '_femurCS.mat';
cnam = '_femurCart.mat';
%
%
% Get All Subject Subdirectories
%
div = uigetdir;
ddir = fullfile(div);
rdir = fullfile(div,'Results');
bdir = fullfile(rdir,'Bone');
cdir = fullfile(rdir, 'Cartilage');
tdir = fullfile(rdir,'Thickness');
gdir = fullfile(rdir,'Grids');
load(fullfile(rdir, 'Full Tibia Full Femur Subject Files.mat'));
ns=size(sd,1);
%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1);    % 0 (false) for left/1 (true) for right
kids = repmat(blanks(ns)',1,5);        % Knee IDs
tribfs = cell(ns,1);    % Sagittal bone triangular mesh connectivity
xyzbfs = cell(ns,1);    % Sagittal bone point coordinates
xyzmnl = zeros(ns,3);   % Minimum individual lateral femur coordinates
xyzmxl = zeros(ns,3);   % Maximum individual lateral femur coordinates
xyzmnm = zeros(ns,3);   % Minimum individual medial femur coordinates
xyzmxm = zeros(ns,3);   % Maximum individual medial femur coordinates
xyzs = zeros(ns,3);     % Range of proximal femur outline
rf = zeros(ns,3);       % Radii of femur condyles

scan(1,:)='FFE';
scan(2,:)='RHO';
scan(3,:)='T2S';
%

% Initialize Overall Minimums and Maximums
%

%
hpi = pi/2;             % Half pi
dpi = 2*pi;             % Double pi
rad2deg = 180/pi;       % Radians to degrees
%
% Loop through MAT Files
%
%% Load Bone Data and Generate min/max and range of Grid
tpw=zeros(ns,3);
zsc=zeros(ns,3);
for j=1:3
    tzrmn = zeros(1,3);     % Minimums
    tzrmx = zeros(1,3);     % Maximums
    for i = 1:ns
        i_str=int2str(i);

        if j==1
            fstr=sd(i).FFE.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end


        bone = load(fullfile(rdir,'Bone',[fstr csnam]));
        cart = load(fullfile(rdir,'Cartilage',[fstr cnam]));

        trifb = mk_tri4f(bone.datfb);       %taken from fcomba.m for scaled data to solve 'xyzsb' and 'trisb' AD 3/3/23
        xyzfb = cell2mat(bone.datfb);
        trifb = tri_fix2(trifb,xyzfb);

        datfc = comb_dat(cart.datlct,cart.datmct,cart.dattct);
        trifc = mk_tri4f(datfc);
        xyzfc = cell2mat(datfc);
        trifc = tri_fix2(trifc,xyzfc);

        eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.trimesh = trifb;']);
        eval(['sd(' i_str ').' scan(j,:) '.femur.bone.full.points = xyzfb;']);

        eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.trimesh = trifc;']);
        eval(['sd(' i_str ').' scan(j,:) '.femur.cart.full.points = xyzfc;']);
        eval(['xyzs(' i_str ',:) = max(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline)-min(sd(' i_str ').' scan(j,:) '.tibia.prox_tib_outline);']);
        %
        % Plot Bone Data
        % %
        %    figure;
        %    view(3);
        %    hold on;
        %    plt_datsl(bone.datfb,'b.-');
        %    axis equal;
        %    xlabel('X','FontSize',12,'FontWeight','bold');
        %    ylabel('Y','FontSize',12,'FontWeight','bold');
        %    zlabel('Z','FontSize',12,'FontWeight','bold');
        %    title(fstr,'FontSize',24,'FontWeight','bold', ...
        %           'Interpreter','none');

        % Transform to Cylindrical Coordinate System

        %
        [th,r,z] = cart2pol(xyzfb(:,1),xyzfb(:,3),xyzfb(:,2));
        id = find(th>hpi);                  % Greater than pi/2
        th(id) = th(id)-dpi;                % Minus 2*pi
        tzr = [th,z,r];
        %
        % Get Minimums and Maximums
        %
        tzrn = min(tzr);
        tzrx = max(tzr);

        itst = tzrn<tzrmn;
        tzrmn(itst) = tzrn(itst);
        itst = tzrx>tzrmx;
        tzrmx(itst) = tzrx(itst);


    end

    tpw(:,j) = xyzs(:,2);
    tpw_mn = mean(tpw(:,j));
    zsc(:,j) = tpw(:,j)./tpw_mn;
    %
    % Get Range and Make a Rectangular Grid
    %
    tzrmn(1) = tzrmn(1)*180/pi;            % Convert to degrees
    tzrmx(1) = tzrmx(1)*180/pi;            % Convert to degrees
    tzmn = floor(tzrmn(1:2));              % Get degrees/millimeters closest to -infinity
    tzmx = ceil(tzrmx(1:2));               % Get degrees/millimeters closest to +infinity
    %
    tzmn(1) = floor(tzmn(1)/2)*2;
    tzmx(1) = ceil(tzmx(1)/2)*2;

    t = (tzmn(1):2:tzmx(1))';              % Theta range in 2 degrees increments
    nt = size(t,1);
    % z = (tzmn(2):2:tzmx(2))';              % Z range in 2 mm increments
    z = (tzmn(2):tzmx(2)+2)';              % Z range in 1 mm increments
    nz = size(z,1);
    [T,Z] = meshgrid(t,z);  % Rectangle grid of square quadrilaterals
    tq = T(:);              % Grid points as a vector
    zq = Z(:);              % Grid points as a vector
    nq = size(tq,1);
    %
    % FROM QUAD CONNECT
    quad = (1:nz-1)';
    quad = [quad quad+nz quad+nz+1 quad+1];
    quad = repmat(quad,nt-1,1);
    addcol = repmat(nz*(0:nt-2),(nz-1)*4,1);
    addcol = reshape(addcol,4,(nt-1)*(nz-1))';
    quad = quad+addcol;
    %
    % Triangular Grid
    %
    ne = size(quad,1);      % Number of elements
    %
    trig = [quad(:,1:3) quad(:,1) quad(:,3:4)]';
    trig = reshape(trig,3,2*ne)';


    for i=1:ns
        i_str=int2str(i);
        gnam = [scan(j,:) '_fgrid08_.mat'];
        eval(['sd(' i_str ').' scan(j,:) '.femur.grid = gnam;']);
        gnam = fullfile(gdir,gnam);   % Grid MAT file
    end
    %
    % Save Data
    %
    save(gnam, 'nq', 'nt', 'nz', 'ne', 'trig','quad', 'tq', 'zq', 'zsc', 'tzrmn', 'tzrmx');

    clear quad addcol nq tq nz T Z nt t tzrmn tzrmx tzrn tzmn tzmx tzrx xyzs
end

for j=1:3
    for i = 1:ns           % Loop through MAT files

        i_str=int2str(i);

        if j==1
            fstr=sd(i).FFE.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==2
            fstr=sd(i).RHO.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        elseif j==3
            fstr=sd(i).T2S.femur.bfnam;
            fstr=[fstr(1:6) fstr(15:17)];

        end
        gnam = [scan(j,:) '_fgrid08_.mat'];
        load(fullfile(gdir,gnam));
        cthk = zeros(nq,1);               % Cartilage thicknesses
        xyzic = zeros(nq,3,1);            % Cartilage intersection
        rq = zeros(nq,1);                 % Radial coordinates
        zqs = zeros(nq,1);                % Polar Z coordinates (scaled)

        eval(['trifb = sd(' i_str ').' scan(j,:) '.femur.bone.full.trimesh;']);
        eval(['xyzfb = sd(' i_str ').' scan(j,:) '.femur.bone.full.points;']);

        eval(['trifc = sd(' i_str ').' scan(j,:) '.femur.cart.full.trimesh;']);
        eval(['xyzfc = sd(' i_str ').' scan(j,:) '.femur.cart.full.points;']);

        %
        % Get Cylindrical Coordinates
        %
        [t,r,z] = cart2pol(xyzfb(:,1),xyzfb(:,3),xyzfb(:,2));
        id = find(t>pi/2);
        t(id) = t(id)-2*pi;
        td = t*rad2deg;      % Angles in degrees
        %
        % Check Triangle Normals and Make Orientations +R
        % Not Necessary For Thickness Calculation, Add this Step to MK_TRI4F.M?
        %
        %    [nt,nz,nr] = tri_norm(trisb,[td z r]);
        %
        %    idn = find(nr<0);
        %    if ~isempty(idn)
        %      trisb(idn,:) = trisb(idn,[1 3 2]);
        %    end
        %
        % Scale Polar "Z" Coordinates
        %

        zs = zsc(i,j)*zq;
        %
        % Interpolate to Scaled "Standard" Grid
        %
        dist=12;
        rg = gridproj(trifb,[td z r],tq,zs,1,dist);
        %
        % Back to Cartesian Coordinates
        %
        [xg,zg,yg] = pol2cart(tq/rad2deg,rg,zs);
        xyzg = [xg yg zg];
        %
        % Improve the Mesh
        %
        trig = tri_fix2(trig,xyzg);
        %
        % Calculate Cartilage Thicknesses
        %
        [ct,bi] = car_thk8(trifc,xyzfc,trig,xyzg);
        %car_thk5(trib,xyzb,tric,xyzc,iplt,fstr,cmprt,hf1);
        % Save Thicknesses and Print Plot
        %
        cthk(:,1) = ct';
        xyzic(:,:,1) = bi;
        % Plot the Cartilage Thicknesses
        %


        hf2 = figure;
        
        hb2 = trimesh(trig,xyzg(:,1),xyzg(:,2),xyzg(:,3), ...
            'LineWidth',0.25,'FaceColor','none', ...
            'EdgeColor','k');
        hold on;
        %
        hc2 = trimesh(trifc,xyzfc(:,1),xyzfc(:,2), ...
            xyzfc(:,3),'LineWidth',0.25, ...
            'FaceColor','none','EdgeColor','b');
        xlabel('X','FontSize',12,'FontWeight','bold');
        ylabel('Y','FontSize',12,'FontWeight','bold');
        zlabel('Z','FontSize',12,'FontWeight','bold');
        title({fstr; 'Sagittal Cartilage and Bone Meshes'}, ...
            'Interpreter','none','FontSize',16,'FontWeight','bold');
        view(-60,12);
        axis equal;

        psnam = fullfile(tdir,[fstr '_fcart08_sagittal_thk.pdf']);
        orient landscape;
        set(hf2, 'units','normalized','outerposition',[0 0 1 1]);
        exportgraphics(hf2, psnam, "Resolution", 300);
        close(hf2);


        % Plot Sagittal Cartilage Thicknesses
        %
        hf4 = figure;
        
        hb4 = trimesh(trig,xyzg(:,1),xyzg(:,2),xyzg(:,3), ...
            'LineWidth',0.25,'FaceColor','none', ...
            'EdgeColor','b');
        hold on;
        %
        ctf = squeeze(cthk(:,1));
        idx = find(~isnan(ctf));
        it = nod2tri(idx,trig,2);
        hs4l = trisurf(trig(it,:),xyzg(:,1),xyzg(:,2),xyzg(:,3),ctf, ...
            'FaceColor','interp', ...
            'EdgeColor','b','LineWidth',0.25);
        colormap jet;
        view(-60,12);
        axis equal;
        %

        if min(ctf)<0
            caxis([min(ctf) max(ctf)]);
        else
            caxis([0 max(ctf)]);
        end
        hcb4 = colorbar;
        xlabel('X','FontSize',12,'FontWeight','bold');
        ylabel('Y','FontSize',12,'FontWeight','bold');
        zlabel('Z','FontSize',12,'FontWeight','bold');
        title({fstr; 'Sagittal Cartilage Thicknesses'}, ...
            'Interpreter','none','FontSize',16,'FontWeight','bold');

        psnam = fullfile(tdir,[fstr '_fcart08_sagittal_thk.pdf']);
        orient landscape
        set(hf4, 'units','normalized','outerposition',[0 0 1 1]);
        exportgraphics(hf4, psnam, "Resolution", 300, 'Append',true);
        %close(hf4);
            
        
        % End if iplt
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
        rnam = fullfile(tdir,[fstr '_fcart08_thk.mat']);    % Results MAT file
        save(rnam,'cthk', 'trig', 'trifc', 'xyzg', 'xyzfc');
        %
        %close all;              % Close cartilage thickness plots
        %
        % if i ~= ns
        %     clearvars -except k ns bpnams bfnams cpnams cfnams scans conds ilegs kids...
        %         iplt iprt
        % end
    end
end
main_file=fullfile(rdir, 'Subject Files Full Femur.mat');
save (main_file, 'sd');
return