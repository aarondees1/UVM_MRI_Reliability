function roi = rd_roi6(filenam,ipx) % Define function to read OSIRIX ROI CSV file
%RD_ROI6  Reads an OSIRIX ROI CSV file.
%
%         ROI = RD_ROI6(FILENAM) reads the OSIRIX ROI CSV file, FILENAM,
%         and returns the structure, ROI, with the names of the ROI in
%         field "name" and the X, Y and Z data points for each slice is
%         in the columns of cell arrays in field "data". The X, Y and Z
%         data is in ordered triplets.  Each triplet represent a point
%         from a slice in the ROI.  Each ROI is returned in a row in
%         the structure.  
%
%         ROI = RD_ROI6(FILENAM,IPX) if IPX is true (nonzero) the X and
%         Y pixel data points for each slice is in the columns of cell
%         arrays in field "data" in the structure ROI.  Also returns the
%         field "imageno" with the image (or slice) number.
%
%         NOTES:  1.  The data was collected from OSIRIX using the
%                 polygon tool.
%
%                 2.  If the ROI name does not start with an alphabetic
%                 letter, the letter "a" is prepended to the ROI name.
%
%         09-Aug-2021 * Mack Gardner-Morse
%
% ######################################################################
%
% Check for Inputs
%
if (nargin<1) % Check if filename input is missing
  error(' *** ERROR in RD_ROI6:  An input file name is required!'); % Throw error for missing filename
end
%
if (nargin<2) % Check if pixel flag input is missing
  ipx = false; % Set default pixel flag to false
end
%
if isempty(ipx) % Check if pixel flag is empty
  ipx = false; % Set pixel flag to false
end
%
% Open File and Read First Two Lines
%
fid = fopen(filenam,'rt'); % Open file for reading
lin = fgetl(fid); % Read first line of headers
hdrs = textscan(lin,'%s','Delimiter',','); % Parse headers
hdrs = hdrs{1}; % Extract header array
idr = find(startsWith(hdrs,'RoiName')); % Find ROI name column index
idn = find(startsWith(hdrs,'NumOfPoints')); % Find number of points column index
lin = fgetl(fid); % Read first line of data
%
% Keep Reading Data Until End of File
%
while lin~=-1 % Loop until end of file
     idx = strfind(lin,','); % Find comma positions
     idx = [idx length(lin)+1]; % Append end of line index
     if ipx % Check if pixel mode is enabled
       img_num = str2double(lin(1:idx(1)-1))+1; % Extract and adjust image number
     end
     rnam = lin(idx(idr-1)+1:idx(idr)-1); % Extract ROI name
%
% Strip Any Quotes from ROI Name
%
     idq = rnam==''''|rnam=='"'; % Identify quotes
     idq = ~idq; % Invert to keep non-quote characters
     rnam = rnam(idq); % Remove quotes from ROI name
%
% Check for Hyphens, Percents and Spaces
%
     htrap = strfind(rnam,'-'); % Find hyphens
     if ~isempty(htrap) % Check if hyphens exist
       rnam(htrap) = '_'; % Replace hyphens with underscores
     end
     htrap = strfind(rnam,'%'); % Find percent signs
     if ~isempty(htrap) % Check if percent signs exist
       rnam(htrap) = 'p'; % Replace percent signs with 'p'
     end
     htrap = strfind(rnam,' '); % Find spaces
     if ~isempty(htrap) % Check if spaces exist
       rnam(htrap) = '_'; % Replace spaces with underscores
     end
%
% Check that First Character is a Letter
%
     if ~isletter(rnam(1)) % Check if first character is not a letter
       rnam = ['a' rnam]; % Prepend 'a' to ROI name
     end
%
% Number of Data Points for this ROI
%
     npts = eval(lin(idx(idn-1)+1:idx(idn)-1)); % Extract number of points
     idp = idn+(0:5:(npts-1)*5)'; % Calculate point data indices
     if ipx % Check if pixel mode is enabled
       mat = zeros(npts,2); % Initialize matrix for pixel data
     else
       mat = zeros(npts,3); % Initialize matrix for 3D data
     end
%
% Get Point Data
%
     for k = 1:npts % Loop through points
        if ipx % Check if pixel mode is enabled
          mat(k,:) = eval(['[' lin(idx(idp(k)+3)+1:idx(idp(k)+5)-1) ']']); % Extract pixel data
        else
          mat(k,:) = eval(['[' lin(idx(idp(k))+1:idx(idp(k)+3)-1) ']']); % Extract 3D data
        end
     end
%
% Save Data for Each ROI
%
     if exist(rnam,'var') % Check if ROI variable exists
       eval(['nslice = size(' rnam ',2);']); % Get number of slices
       eval([rnam '{1,' int2str(nslice+1) '} = mat;']); % Append data to ROI
       if ipx % Check if pixel mode is enabled
         eval([rnam '{2,' int2str(nslice+1) '} = img_num;']); % Append image number
       end
     else
       eval([rnam '{1,1} = mat;']); % Initialize ROI data
       if ipx % Check if pixel mode is enabled
         eval([rnam '{2,1} = img_num;']); % Initialize image number
       end
     end
     lin = fgetl(fid); % Read next line
end
%
% Close File
%
fclose(fid); % Close file
%
% Clear Workspace
%
clear ans fid filenam hdrs htrap idn idq idr lin k idx img_num ipx ...
      nslice rnam npts idp mat; % Clear temporary variables
%
% Get ROI
%
rnam = who; % Get ROI variable names
nvar = size(rnam,1); % Get number of ROIs
ipx = eval(['size(' rnam{1} ',1)'])>1; % Check if pixel mode is active
%
% Get Data into a Cell
%
data = cell(nvar,1); % Initialize data cell array
if ipx % Check if pixel mode is enabled
  img_num = cell(nvar,1); % Initialize image number cell array
end
%
for k = 1:nvar % Loop through ROIs
  data{k} = eval([ rnam{k} '(1,:)']); % Extract data
  if ipx % Check if pixel mode is enabled
    img_num{k} = eval([ '[' rnam{k} '{2,:}]']); % Extract image numbers
  end
end
%
% Put Information into a Structure
%
roi = struct('name',rnam,'data',data); % Create ROI structure
if ipx % Check if pixel mode is enabled
  [roi.imageno] = deal(img_num{:}); % Add image numbers to structure
end
%
return % Exit the function