function [m,fileheaders]=test_loadsnotel(interval, siteID)
% test_loadsnotel.m
%
% Generates a loadsnotel-type matrix without any modification of original
% data. This can be used to compare to datasets from the same site.interval
% that have undergone some data cleaning or filtering (using
% plot_snoteltests.m for example)
%
% FIXME - needs methods to add testing data (dummy values in rows/files)

%fignum=0; close all;

% ------------------------------------------------------------------------
% select daily files (all sensors) or hourly soil sensor data
if strcmp(interval, 'daily')
    datapath = '../datafiles/allsensors_daily/';
    %datapath = '/home/greg/data/rawdata/SNOTELdata/allsensors_daily/';
    %datapath = '../datafiles/NiwotSNOTELs_90-2010/';
    
    % A complete datafile has this many sensors
    completesensorset = 20;
    % Columns to find in daily files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'WTEQ.I-1' 'PREC.I-1' 'TOBS.I-1' ...
    'TMAX.D-1' 'TMIN.D-1' 'TAVG.D-1' 'SNWD.I-1' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20' 'RDC.I-1:-2 ' ...
    'RDC.I-1:-8' 'RDC.I-1:-20' 'BATT.I-1'};

elseif strcmp(interval, 'hourly')
    datapath = '../datafiles/soilsensors_hourly/';
    %datapath = '/home/greg/data/rawdata/SNOTELdata/soilsensors_hourly/';
    
    % A complete datafile has this many sensors
    completesensorset = 9;
    % Columns to find in hourly files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20'};
else
    error('Not a valid data type (daily or hourly)')
end

% Set dummy flag on if siteID = 0
% This will trigger dummy generator at end of file
if siteID == 0
    siteID = 828;
    dummysite = 1;
else
    dummysite = 0;
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
templates = {t t t t t t};
[fileheaders collectindices desiredcols headerspresent fillcolumnorder ...
    missingdatacols] = templates{:};

% First, we will loop through each file and identify what sensors are
% present (and in which columns), and which are missing(and in which col)
for i = 1:length(files)
    % read file headers to cell of strings - (1 per year of data)
    fid = fopen(files{i});
    % use fgetl(fid) to get the second line of the file
    fgetl(fid); % skip the first line
    fileheaders(i) = textscan(fgetl(fid), '%s', 'Delimiter', ',');
    fclose(fid);
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
    format = ['%d%s%s' forma(6:end)];
    
    % Load data file into raw cell array
    readfilecell = textscan(fid, format, 'Headerlines', 2, 'Delimiter', ',');
    fclose(fid);
    
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



% %--------------------------------------------------------------------------
% % CREATE DUMMYSITE if called for -  for QA testing of code
% %
% if strcmp(interval, 'hourly') && dummysite == 1;
%     % find all July time periods
%     test = tvec(:,2) == 7;
%     % get the actual T data for this site (2 in)
%     tmp = m{7};
%     % modify the July periods to be exactly 4 degress
%     tmp(test,:) = 4;
%     % replace T data in m with dummy data
%     m{7} = tmp;
% end

junk = 99;