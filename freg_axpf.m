%#######################################################################
%
%          * Femur REGions using AXial Plane Female Program *
%
%          M-File which reads the female femur "standard" grid and
%     divides the femur into six regions based on an axial plane view of 
%     the trochlea.  The regions are plotted and the points in each
%     region are written to the MS-Excel spreadsheet:
%     Partitions6_f2.xlsx.
%
%
%     NOTES:  1.  Femur MAT files *sc.mat and tibia MAT files *CS.mat
%             must be in the directory "BM_Female_Femur_Mat_Files".  The
%             MAT file fb_rngf.mat must be in the current directory.
%
%             2.  The M-files in_tri2d.m, meshbnd4.m, nod2tri.m,
%             sl_info.m and xzpl2pol.m must be in the current path or
%             directory.
%
%             3.  To average the trochlea regions, run the M-file:
%             fregs_bmf.m.
%
%     26-Sep-2016 * Mack Gardner-Morse
%

%#######################################################################
%
% Plots?
%
iplt = true;
% iplt = false;
iplt2 = true;           % Debugging plots
% iplt2 = false;          % Debugging plots
% iprt = true;
iprt = false;
%
% Convert Radians to Degrees and Coordinate Tolerance
%
rad2deg = 180/pi;
tol = 1e-10;            % Tolerance on grid coordinates
%
% MS-Excel Spreadsheet Name
%
xnam = 'Partitions6_f2.xlsx';
%
% Data Directory and File Names
%
fdir = 'BM_Female_Femur_MAT_Files';
d = dir(fullfile(fdir,'*EAbc7tsc.mat'));
fnams = deblank(char(d.name));
nf = size(fnams,1);   % Number of files
%
tdir = 'BM_Female_Femur_MAT_Files';
tnams = dir(fullfile(tdir,'*CS.mat'));
tnams = char(tnams.name);
%
% Get Tibia Widths
%
nft = size(tnams,1);
tpw = zeros(nft,1);
for k = 1:nft
   fnamt = tnams(k,:);
   load(fullfile(tdir,fnamt));
   tpw(k) = max(xyzpt(:,2))-min(xyzpt(:,2));
end
tnams = tnams(:,1:5);
tpw_mn = mean(tpw);
%
% Get Standard Grid
%
load fb_rngf.mat;
%
% Posterior and Medial/Lateral Cutoffs
%
tmin = -145;            % Angle (theta) cutoff (-150, -145, or -140)
y0 = 0;                 % Y cutoff (-1, 0, or 1)
%
% Get Indices to the Posterior and Medial/Lateral Regions
%
ipos = find(tq<=tmin);  % Posterior
idl = find(zq>=y0);     % Lateral
idm = find(zq<y0);      % Medial
%
% Loop through MAT Files
%
for kf = 1:nf           % Loop through MAT files
% for kf = 9:nf           % Loop through MAT files
%
   fnam = fnams(kf,:);                 % File name
   idot = strfind(fnam,'.');
   fstr = fnam(1:idot-1);              % File name w/o extension
%
% Load Data
%
   load(fullfile(fdir,fnam));
%
% Combine Bone Data and Transform Coordinates to Femur Coordinate System
%
   datb = datsb;
   [nsl,ns,is] = sl_info(datb);
%
% Get Maximum "Peaks" of Femur in X-Direction
%
   mx = zeros(nsl,1);
   idx = zeros(nsl,1);
   my = zeros(nsl,1);
   for ks = 1:nsl
      [mx(ks),idx(ks)] = max(datb{ks}(:,1));
      my(ks) = datb{ks}(idx(ks),2);
   end
%
   [mxs ids] = sort(mx);
   idp = find(my(ids)>0);
   idp = ids(idp(end-2:end));          % Medial peak index to top three
   idn = find(my(ids)<0);
   idn = ids(idn(end-2:end));          % Lateral peak index to top three
   mxp = mean(mx(idp));                % Medial peak X
   myp = mean(my(idp));                % Medial peak Y
   mxn = mean(mx(idn));                % Lateral peak X
   myn = mean(my(idn));                % Lateral peak Y
%
% Get Trochlea Groove
%
   idxmn = [round(mean(idp)); round(mean(idn))];
   idxmn = sort(idxmn);
   idmn = find(ids>idxmn(1)&ids<idxmn(2));
   idmn = ids(idmn(1:3));
   mxm = mean(mx(idmn));
   mym = mean(my(idmn));
%
% Get Line Directions
%
   xyc = [mxm mym];
   xym = [mxp myp];
   xyl = [mxn myn];
   vm = xym-xyc;
   vl = xyl-xyc;
%
% Get Center of Trochlea
%
   xyc = [0 0];         % Origin
%
% Get Trochlea Points for Each Slice
%
   ip = NaN(nsl,2);
   idx = [];
   for ks = 1:nsl
      xys = mean(datb{ks}(:,1:2));
      mxs = [mx(ks) my(ks)];
      vx = 5*(xys-mxs);
      yrng = sort([my(ks); xys(2)+vx(2)]);
      t = 1.2*yrng(2)./vm(2);
      if t>=0
        ip(ks,:) = lsect3(xyc(1,1:2),t*vm,mxs,vx)';
        if iplt2
          h1 = plot([xyc(1); xyc(1)+t*vm(1)],[xyc(2); xyc(2)+t*vm(2)],'r-','LineW',2);
          hold on;
          h2 = plot([mxs(1); mxs(1)+vx(1)],[mxs(2); mxs(2)+vx(2)],'g-','LineW',2);
          h3 = plot(ip(ks,1),ip(ks,2),'ks','LineW',2);
          pause;
          delete([h1,h2,h3]);
        end
      end
      if isnan(ip(ks,1))
        t = 1.2*yrng(1)./vl(2);
        if t>=0
          ip(ks,:) = lsect3(xyc(1,1:2),t*vl,mxs,vx)';
          if iplt2
            h1 = plot([xyc(1); xyc(1)+t*vl(1)],[xyc(2); xyc(2)+t*vl(2)],'r-','LineW',2);
            hold on;
            h2 = plot([mxs(1); mxs(1)+4*vx(1)],[mxs(2); mxs(2)+4*vx(2)],'g-','LineW',2);
            h3 = plot(ip(ks,1),ip(ks,2),'ms','LineW',2);
            pause;
            delete([h1,h2,h3]);
          end
        end
      end
      if ~isnan(ip(ks,1))
        id = find(datb{ks}(:,1)>sign(mx(ks))*ip(ks,1));
        ids = (is(ks)+1:is(ks+1))';
        id = ids(id);
        idx = [idx; id];
      end
   end
%
% Get Femur Mesh
%
   trib = trisb;
   xyzb = xyzsb;
%
% Plot Femur in Cartesian Coordinates
%
   if iplt
     hfc = figure;
     orient landscape;
     view(3);
     hold on;
     hm = trimesh(trib,xyzb(:,1),xyzb(:,2),xyzb(:,3),'FaceColor', ...
                  'none','EdgeColor','b','LineWidth',0.5);
     ht = plot3(xyzb(idx,1),xyzb(idx,2),xyzb(idx,3),'ro', ...
                'LineWidth',0.5);
     axis equal;
     xlabel({'X (mm)';['\leftarrow posterior / anterior', ...
            ' \rightarrow']},'FontSize',12,'FontWeight','bold');
     ylabel({'Y (mm)';'\leftarrow lateral / medial \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     zlabel({'Z (mm)';'\leftarrow inferior / superior \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     title(fstr,'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       print('-dpsc2','-r300',[fstr(1:5) '.ps']);
     end
   end
%
% Convert to Polar Coordinates
%
   datbp = xzpl2pol(datb);
   trzb = cell2mat(datbp);
   trzb(:,1) = trzb(:,1)*rad2deg;
%
% Scale "Z" Coordinates
%
   tqs = tq;            % No scaling of theta
   iz = strmatch(fstr(1:5),tnams);
   zsc = tpw(iz)./tpw_mn;
   zqs = zsc*zq;
%
% Find Boundary Points on Trochlea
%
   ir = nod2tri(idx,trib,2);
   ina = in_tri2d(trib,trzb(:,[1 3]),[tqs zqs]);
   ina = find(ina);
   in = in_tri2d(trib(ir,:),trzb(:,[1 3]),[tqs zqs]);
   in = find(in);
   itta = nod2tri(ina,trig,2);
   itt = nod2tri(in,trig,2);
   bida = meshbnd4(trig(itta,:));
   bid = meshbnd4(trig(itt,:));
   idb = setdiff(bid,bida);
   idbl = intersect(idl,idb);          % Grid points in lateral compartment
   idbm = intersect(idm,idb);          % Grid points in medial compartment
%
% Get Unique Z Values in the Lateral and Medial Compartments
%
   zs = sort(zqs);
   izu = find(diff(zs)>tol);
   zu = zs([izu; izu(end)+1]);         % Unique Z values
   nu = size(zu,1);
%
   izl = find(zu>=y0);
   zul = zu(izl);
   nul = size(zul,1);
%
   izm = find(zu<y0);
   zum = zu(izm);
   num = size(zum,1);
%
% Extend Lines Based on Boundary Points
%
   vl = [zqs(idbl) ones(size(idbl))];
   lpl = vl\tqs(idbl);                 % Slope and intercept of theta as a function of Z
   imxl = find(zu>max(zqs(idbl)));
   tcol = [tqs(idbl); [zu(imxl) ones(size(imxl))]*lpl];    % Theta cutoffs
   zcol = [zqs(idbl); zu(imxl)];
%
   vm = [zqs(idbm) ones(size(idbm))];
   lpm = vm\tqs(idbm);                 % Slope and intercept of theta as a function of Z
   imxm = find(zu<min(zqs(idbm)));
   tcom = [tqs(idbm); [zu(imxm) ones(size(imxm))]*lpm];    % Theta cutoffs
   zcom = [zqs(idbm); zu(imxm)];
%
% Find Grid Points with Theta > Theta Cut Offs
%
   ztl = sortrows([zcol tcol]);
   [zlu izlu] = unique(ztl(:,1));
   ztl = ztl(izlu,:);
   idtl = [];
   for kl = 1:nul
      zuc = zul(kl);
      idlc = find(ztl(:,1)>zuc-tol&ztl(:,1)<zuc+tol);
      if isempty(idlc)
        tuc = [zuc 1]*lpl;
      else
        tuc = ztl(idlc,2);
        tuc = tuc(1);
      end
      idtl = [idtl; find(zqs>zuc-tol&zqs<zuc+tol&tqs>=tuc)];
   end
%
   ztm = sortrows([zcom tcom]);
   [zmu izmu] = unique(ztm(:,1));
   ztm = ztm(izmu,:);
   idtm = [];
   for km = 1:num
      zuc = zum(km);
      idmc = find(ztm(:,1)>zuc-tol&ztm(:,1)<zuc+tol);
      if isempty(idmc)
        tuc = [zuc 1]*lpm;
      else
        tuc = ztm(idmc,2);
        tuc = tuc(1);
      end
      idtm = [idtm; find(zqs>zuc-tol&zqs<zuc+tol&tqs>=tuc)];
   end
%
% Plot Femur in Polar Coordinates
%
   if iplt
     hfp = figure;
     orient landscape;
     view(2);
     hold on;
     hm = trimesh(trib,trzb(:,1),trzb(:,3),'Color','b', ...
                  'LineWidth',0.5);
     hg = plot(tqs,zqs,'k.','MarkerSize',6,'LineWidth',0.5);
     hql = plot(tqs(idtl),zqs(idtl),'ro','MarkerSize',4, ...
                'LineWidth',0.5);
     hqm = plot(tqs(idtm),zqs(idtm),'go','MarkerSize',4, ...
                'LineWidth',0.5);
   end
%
% Get Remaining Indices to the Six (6) Regions
%
   iant = [idtl; idtm]; % Anterior (Trochlea)
   ictr = true(nq,1);
   ictr(iant) = false;
   ictr(ipos) = false;
   ictr = find(ictr);   % Center
%
   idl_pos = intersect(idl,ipos);      % Lateral posterior region
   idl_ctr = intersect(idl,ictr);      % Lateral center region
   idl_ant = idtl;                     % Lateral anterior region
%
   idm_pos = intersect(idm,ipos);      % Medial posterior region
   idm_ctr = intersect(idm,ictr);      % Medial center region
   idm_ant = idtm;                     % Medial anterior region
%
% Plot Regions
%
   if iplt
     figure(hfp);
     plot(tqs(idl_pos),zqs(idl_pos),'y^','MarkerSize',4, ...
          'LineWidth',0.5);
     plot(tqs(idm_pos),zqs(idm_pos),'c^','MarkerSize',4, ...
          'LineWidth',0.5);
     axis tight;
     xlabel({'\Theta (degrees)';['\leftarrow posterior / anterior', ...
            ' \rightarrow']},'FontSize',12,'FontWeight','bold');
     ylabel({'Z (mm)';'\leftarrow medial / lateral \rightarrow'}, ...
            'FontSize',12,'FontWeight','bold');
     title(fstr,'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       print('-dpsc2','-r300','-append',[fstr(1:5) '.ps']);
     end
     pause;
   end
%
% Create Index for the Six Regions (1-6)
%
   id = NaN(nq,1);
   id(idl_pos) = 1;
   id(idl_ctr) = 2;
   id(idl_ant) = 3;
   id(idm_pos) = 4;
   id(idm_ctr) = 5;
   id(idm_ant) = 6;
%
% Write Index to Spreadsheet
%
%   lbl = mat2cell([repmat('Pt ',nq,1) int2str((1:nq)')],ones(nq,1));
%   xlswrite(xnam,lbl,fstr(1:5),'A1');
%   xlswrite(xnam,id,fstr(1:5),'B1');
%
end
%
return