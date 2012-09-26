function m = loadsnotel_raw(interval, siteID, varargin)
% loadsnotel_raw.m
%
% Generates a loadsnotel-type matrix with no modification of original
% data. This can be used to compare to datasets from the same site/interval
% that have undergone some data cleaning or filtering (using
% plot_snoteltests.m for example).
%
% If varargin is used, columns can be deleted or dummy values introduced
%
% arguments:
% interval and siteID are the same as in loadsnotel.m
% varargin = cell array containing a test type ('delete' or 'dummyvalues')
% and a column number referring to the desired header/column number (1-20)

% -----------------------------------------------------------------------
% select daily files (all sensors) or hourly soil sensor data
if strcmp(interval, 'daily')
    datapath = '../rawdata/allsensors_daily/';
    %datapath = '/home/greg/data/rawdata/NRCSdata/allsensors_daily/';
    %datapath = '../rawdata/NiwotSNOTELs_90-2010/';
    
    % A complete datafile has this many sensors
    completesensorset = 20;
    % Columns to find in daily files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'WTEQ.I-1' 'PREC.I-1' 'TOBS.I-1' ...
    'TMAX.D-1' 'TMIN.D-1' 'TAVG.D-1' 'SNWD.I-1' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20' 'RDC.I-1:-2 ' ...
    'RDC.I-1:-8' 'RDC.I-1:-20' 'BATT.I-1'};

elseif strcmp(interval, 'hourly')
    datapath = '../rawdata/soilsensors_hourly/';
    % datapath = '/home/greg/data/rawdata/SNOTELdata/soilsensors_hourly/';
    
    % A complete datafile has this many sensors
    completesensorset = 9;
    % Columns to find in hourly files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20'};
else
    error('Not a valid data type (daily or hourly)')
end

testtype = 'null';
% TEST - Set testing values if varargin is present - TEST
% This will trigger deletions or a dummy value generator in file
if size(varargin,2) > 0
    testtype = varargin{1};
    testcolumn = varargin{2};
end

%--------------------------------------------------------------------------
% create a cellarray of filenames for the site
try
    filelist = textscan(ls([datapath num2str(siteID) '*']), '%s'); 
catch exception
    error('missing:SiteData', ['no data found for site ' num2str(siteID)]); 
end
files = sort(filelist{1});

%--------------------------------------------------------------------------
% Sensors, and therefore columns vary from year to year, site to site,
% therefore we must read header of each file and select desired columns.
%
% These data arrays must then be put in a standardized cell array format

% Preallocate sensor/column sorting cellarrays
t = cell(length(files), 1);
templates = {t t t t t t t};
[fileheaders collectindices desiredcols headerspresent fillcolumnorder ...
    missingdatacols testindices] = templates{:}; %TEST-add deleteindices

% First, we will loop through each file and identify what sensors are
% present (and in which columns), and which are missing(and in which col)
for i = 1:length(files)
    % read file headers to cell of strings - (1 per year of data)
    fid = fopen(files{i});
    % use fgetl(fid) to get the second line of the file
    fgetl(fid); % skip the first line
    fileheaders(i) = textscan(fgetl(fid), '%s', 'Delimiter', ',');
    fclose(fid);
    
    % TEST - mark test columns and remove header if delete test - TEST
    if strcmp(testtype, 'delete')
        testindices{i} = strncmpi(desiredheaders{testcolumn}, ...
            fileheaders{i}', length(desiredheaders{testcolumn}));
        fileheaders{i}(testindices{i}) = [];
    elseif strcmp(testtype, 'dummyvalues')
        testindices{i} = strncmpi(desiredheaders{testcolumn}, ...
            fileheaders{i}', length(desiredheaders{testcolumn}));
    end

    % Now we have the headers from file(i) - next loop through each desired
    % header, compare to file headers, and make arrays to identify the
    % presence/absence and location of each
    for j = 1:length(desiredheaders)
        explength = length(desiredheaders{j});
        % Get logical index of which column header matches desiredheader j
        desiredcolidx = strncmpi(desiredheaders{j}, fileheaders{i}', ...
            explength);
        % Collect these logical indices for each desired header{i} - each 
        % row represents one desired header with a logical index of the 
        % column location in the file
        collectindices{i} = [collectindices{i}; desiredcolidx];
        % Place column number containing each desired sensor (in order)
        % into an array
        desiredcols{i} = [desiredcols{i}, find(desiredcolidx)];
    end
    % Row-wise sum of collectindices{i} gives a logical array of the
    % presence/absence of all desired headers in the file
    headerspresent{i} = logical(sum(collectindices{i}, 2)');
    % Based on what headers are present in the input file (headerspresent),
    % fill columns of output matrix with data in this order
    fillcolumnorder{i} = [fillcolumnorder{i}; find(headerspresent{i})];
    % Columns that have missing data (no header in file) are listed in this
    % array (so they can be changed to nan later)
    missingdatacols{i} = find(~headerspresent{i});
end

clear filelist templates collectindices headerspresent i j  test; 


% -------------------------------------------------------------------------
% Using the information about each file's contents from the arrays above,
% it is now time to read in the data one file at a time. As each file is
% read, the desired columns will be added in a pre-determined order to a
% new array with a standardized format (desiredheaders)
%
% First, create the cell array to be returned (m) and an intermediate one
% for holding the data from each new file
orderedcell = cell(1, completesensorset);
m = cell(1, completesensorset);

% Read in data one file at a time
for i = 1:length(files)
    fid = fopen(files{i});
    
    % Create format string to read in original file
    form = num2str(ones(1, length(fileheaders{i})));
    forma = regexprep(form, '  1', '%f');
    
    % TEST - if delete test called add in the removed column (keeps the
    % columns in the correct order, but will be removed later)
    if strcmp(testtype, 'delete')
        format = ['%d%s%s' forma(6:end) '%f'];
    else
        format = ['%d%s%s' forma(6:end)]; % normal format string (no TEST)
    end
    
    % Load data file into raw cell array
    readfilecell = textscan(fid, format, 'Headerlines', 2, 'Delimiter', ',');
    fclose(fid);
    
    % TEST - delete or set dummy values in test column - TEST
    if strcmp(testtype, 'delete')
        readfilecell(testindices{i}) = [];
    elseif strcmp(testtype, 'dummyvalues')
        if strcmp(interval, 'hourly')
            tvec = datevec(strcat(readfilecell{2}, readfilecell{3}), ...
                'yyyy-mm-ddHH:MM');
        elseif strcmp(interval, 'daily')
            tvec = datevec(readfilecell{2}, 'yyyy-mm-dd');
        end
        dummytest = tvec(:,2)==7; % select all july values
        readfilecell{testindices{i}}(dummytest, :) = 1.00; % set to 1
    end
    
    % Read data from readfilecell to correct columns in orderedcell...
    loopindex = 1;
    for j = fillcolumnorder{i} 
        % j iterates over which NEW columns to add data to
        % desiredcols iterates over which INPUTFILE columns to use
        % remember to index desiredcols by loop iteration (with loopindex)
        orderedcell{j} = readfilecell{desiredcols{i}(loopindex)};
        
        % Then append the column to m (output array)
        m{j} = [m{j}; orderedcell{j}];
        
        % increment loopindex
        loopindex = loopindex + 1;
    end
    clear readfilecell loopindex j test form forma format;
    
    len = length(orderedcell{1});
    
    % Insert nans for missing sensor columns
    for j = missingdatacols{i}
        orderedcell{j} = nan * zeros(len, 1);
        % ...and then append it to m (output array)
        m{j} = [m{j}; orderedcell{j}];
    end
    clear orderedcell;
end
% junk = 99
