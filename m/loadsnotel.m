function m = loadsnotel(interval, siteID)
% loadsnotel.m
%
% if headers are needed use:
% function [m,fileheaders]=loadsnotel(interval, siteID)
%
% Parses headers and loads desired fields from NWCC (SNOTEL/SCAN) yearly 
% datafiles into one array, removes bad data, returns array + headers.
%
% arguments:
% interval = 'hourly' or 'daily' ... load hourly or daily data files
% siteID = 3-4 digit int ... site identifier
% ------------------------------------------------------------------------

siteID = uint16(siteID); % Convert input to int just in case

% select daily files (all sensors) or hourly soil sensor data 
if strcmpi(interval, 'daily')
    datapath = '../rawdata/allsensors_daily/';
    %datapath = '/home/greg/data/rawdata/NRCSdata/allsensors_daily/';
    %datapath = '../rawdata/other/NiwotSNOTELs_90-2010/';
    
    % A complete datafile has this many sensors
    completesensorset = 20;
    % Columns to find in daily files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'WTEQ.I-1' 'PREC.I-1' 'TOBS.I-1' ...
    'TMAX.D-1' 'TMIN.D-1' 'TAVG.D-1' 'SNWD.I-1' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20' 'RDC.I-1:-2 ' ...
    'RDC.I-1:-8' 'RDC.I-1:-20' 'BATT.I-1'};
    % Load list of bad data years for all sites
    % badDataList = sortrows(csvread([datapath 'baddata.txt'], 1, 0));
    % badYears = badDataList((badDataList(:, 1)==siteID), 2);


elseif strcmpi(interval, 'hourly')
    datapath = '../rawdata/soilsensors_hourly/';
    % datapath = '/home/greg/data/rawdata/SNOTELdata/soilsensors_hourly/';
    
    % A complete datafile has this many sensors
    completesensorset = 9;
    % Columns to find in hourly files (note spaces after -2 sensors)
    desiredheaders = {'Site Id' 'Date' 'Time' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20'};
    % badYears = [];
else
    error('Not a valid data type (daily or hourly)')
end

% Until there is better baddata control, use this file:
badDataList = sortrows(csvread(...
    '../rawdata/allsensors_daily/baddata.txt', 1, 0));
badYears = badDataList((badDataList(:, 1)==siteID), 2);

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

%--------------------------------------------------------------------------
% CLEAN OUT BAD DATA

% Start with hourly files-----
if strcmpi(interval, 'hourly')
    
    % Mark and remove sub-hourly rows
    tvec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    subhourly = tvec(:, 5) ~= 0;
    tvec(subhourly,:) = [];
    clear test;
    
    % remove rows marked in subhourly array
    for i = 1:completesensorset
        m{i}(subhourly) = [];
    end
    
    % Create a water-year vector and add to end of matrix
    wyearvec = tvec(:, 1);
    wytest = (tvec(:,2)==10 | tvec(:,2)==11 | tvec(:,2)==12);
    wyearvec(wytest) = wyearvec(wytest) + 1;
    m{10} = wyearvec;
    
    % Test to remove badyears
%     badyeartest = ismember(wyearvec, badYears);
    
    for i = 4:completesensorset
        
        % mark the datalogger error signals (-99.9, -273.2)
        test = m{i} == -99.9 | m{i} == -273.2;
        m{i}(test) = nan;
        clear test;
        
        % Remove the bad data in badyeartest
%         m{i}(badyeartest) = nan;
        
%         % change unreasonably high and low numbers to nan
%         test = m{i} < -25 | m{i} > 100;
%         m{i}(test) = nan;
%         clear test;
        
        % Remove bad soil moisture values
        if i<7
            % Not lower than 0, not greater than 45 - set by calibration 
            % equation for Hydraprobes (see Seyfried et al 2010)
            test = m{i} > 45.0 | m{i} < 0.0;
            m{i}(test) = nan;
            clear test;
        end
        
        % Remove bad soil temperature values
        if i>6           
            % Get rid of silly values
            test = m{i} > 50 | m{i} < -25;
            m{i}(test) = nan;
            clear test
        end
    end
    
%     % Ensure that each sensor set to nans in the same places
%     sensor1_nantest = (isnan(m{4}) | isnan(m{7}));
%     m{4}(sensor1_nantest) = nan;
%     m{7}(sensor1_nantest) = nan;
%     sensor2_nantest = (isnan(m{5}) | isnan(m{8}));
%     m{5}(sensor2_nantest) = nan;
%     m{8}(sensor2_nantest) = nan;
%     sensor3_nantest = (isnan(m{6}) | isnan(m{9}));
%     m{6}(sensor3_nantest) = nan;
%     m{9}(sensor3_nantest) = nan;
    
    clear i;
    

% Filter daily files-----
elseif strcmpi(interval, 'daily')
    
    % Mark and remove sub-daily rows
    tvec = datevec(m{2}, 'yyyy-mm-dd');
    subdaily = ~strcmp(m{3}, '');
    tvec(subdaily,:) = [];
    clear test;

    % Remove rows marked in subdaily array
    for i = 1:completesensorset
        m{i}(subdaily) = [];
    end
        
    % Create a water-year vector and add to end of matrix
    wyearvec = tvec(:, 1);
    wytest = (tvec(:,2)==10 | tvec(:,2)==11 | tvec(:,2)==12);
    wyearvec(wytest) = wyearvec(wytest) + 1;
    m{21} = wyearvec;
    
    % Test to remove badyears
    badyeartest = ismember(wyearvec, badYears);
    
    for i = 4:completesensorset
        % mark the datalogger error signals (-99.9, -273.2)
        test = m{i} == -99.9 | m{i} == -273.2;
        m{i}(test) = nan;
        clear test;
        
        % Remove the bad data in badyeartest
        m{i}(badyeartest) = nan;
        
        % Set negative WTEQ and PREC values to 0
        if i<6
            test = m{i} < 0.0;
            m{i}(test) = 0;
            clear test;
        end
        
        % Set negative snow depth values to 0
        if i>9 && i<11
            test = m{i} < 0.0;
            m{i}(test) = 0.0;
            clear test;
        end
        
        % Remove high soil moisture values
        if i>10 && i<14
            test = m{i} > 100;
            m{i}(test) = nan;
            clear test;
        end
        
        % Remove High temp values for deep soil sensors
        if i>13 && i<17
            test = m{i} > 35.0;
            m{i}(test) = nan;
            clear test;
        end
        
    end
    clear i;
    
    sensor1_nantest = (isnan(m{6}) | isnan(m{7}) |isnan(m{8}) |isnan(m{9}));
    m{6}(sensor1_nantest) = nan;
    m{7}(sensor1_nantest) = nan;
    m{8}(sensor1_nantest) = nan;
    m{9}(sensor1_nantest) = nan;
    sensor2_nantest = (isnan(m{11}) | isnan(m{14}) | isnan(m{17}));
    m{11}(sensor2_nantest) = nan;
    m{14}(sensor2_nantest) = nan;
    m{17}(sensor2_nantest) = nan;
    sensor3_nantest = (isnan(m{12}) | isnan(m{15}) | isnan(m{18}));
    m{12}(sensor3_nantest) = nan;
    m{15}(sensor3_nantest) = nan;
    m{18}(sensor3_nantest) = nan;
    sensor4_nantest = (isnan(m{13}) | isnan(m{16}) | isnan(m{19}));
    m{13}(sensor4_nantest) = nan;
    m{16}(sensor4_nantest) = nan;
    m{19}(sensor4_nantest) = nan;
end