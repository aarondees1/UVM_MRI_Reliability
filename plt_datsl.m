function p = plt_datsl(dat,cmstr,lw,colr,lbl)
% PLT_DATSL Plots ROI data one slice at a time in 3D.
%
%          P = PLT_DATSL(DAT) given a column cell array containing
%          three (3) column matrices with slice coordinate point data,
%          DAT, returns a 3D plot with a default black line for each
%          slice.
%
%          P = PLT_DATSL(DAT,CMSTR) given a cell array containing three
%          (3) column matrices with slice coordinate point data in the
%          cell array, DAT, returns a 3D plot with a marker and color
%          designation of CMSTR.
%
%          P = PLT_DATSL(DAT,CMSTR,LW) adds a designated line width to
%          the integer, LW.
%
%          P = PLT_DATSL(DAT,CMSTR,LW,COLR) adds a designated RGB color
%          to the 1x3 vector, COLR.
%
%          25 July 2012 Daniel R. Sturnick
%

% Check if cmstr is provided; if not, set default to black ('k').
if nargin<2||isempty(cmstr)
    cmstr = 'k';
end

% Check if line width (lw) is provided; if not, set default to 0.5 and no custom color.
if nargin<3||isempty(lw)
    lw = .5;
    cl = false;
end

% Check if custom color (colr) is provided; set color flag accordingly.
if nargin<4||isempty(colr)
    cl = false;
else
    cl = true;
end

% Check if label flag (lbl) is provided; if not, default to false (no labels).
if nargin<5||isempty(lbl)
    lbl = false;
end

% Verify that input data is a cell array containing slice data.
if iscell(dat)
    nslice = length(dat); % Get number of slices in the cell array.
    
    % Allocate output array for graphics handles based on label flag.
    if lbl
        p = zeros(2*nslice,1); % Double size for plot and text handles.
    else
        p = zeros(nslice,1); % Single size for plot handles only.
    end

    % Loop through each slice to plot its 3D data.
    for s = 1:nslice
        xyz = dat{s}; % Extract Nx3 matrix of (x,y,z) coordinates for current slice.
        
        % Plot slice with or without custom color based on cl flag.
        if cl
            l = plot3(xyz(:,1),xyz(:,2),xyz(:,3),cmstr,'LineWidth',lw, ...
                'MarkerSize',10,'Color',colr); % Plot with custom RGB color.
        else
            l = plot3(xyz(:,1),xyz(:,2),xyz(:,3),cmstr,'LineWidth',lw, ...
                'MarkerSize',10); % Plot with default color from cmstr.
        end

        % If labels are requested, add text at first point of slice.
        if lbl
            t = text(xyz(1,1),xyz(1,2),xyz(1,3),int2str(s),'FontSize', ...
                12,'FontWeight','bold'); % Add bold slice number label.
            id = 2*s; % Calculate index for storing handles.
            p(id-1:id) = [l; t]; % Store plot and text handles.
        else
            p(s) = l; % Store only plot handle.
        end

        % Enable hold on after first slice to overlay subsequent plots.
        if s==1
            hold on;
        end
    end
else
    % If input is not a cell array, display error and return NaN.
    fprintf(1,'\n No slices to plot\n\n');
    p = NaN;
end
return