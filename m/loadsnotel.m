function [m,headers]=loadsnotel(interval, site_id)
% loadsnotel.m
%
% Parses headers and loads desired fields from NWCC (SNOTEL/SCAN) yearly 
% datafiles into one array, removes bad data, returns array + headers.
%
% arg 1 = 'hourly' or 'daily' ... load hourly or daily data files
% arg 2 = 3-4 digit int ... site identifier (0 is dummy site)
%
% v2: 2/22/2011 Greg Maurer (changed to loadsnotel.m)
% added dummy site generator and plot of original/removed datapoints
% ------------------------------------------------------------------------
% select daily files (all sensors) or hourly soil sensor data
if strcmp(interval, 'daily')
    datapath = '../datafiles/allsensors_daily/';
    %datapath = '/home/greg/data/rawdata/SNOTELdata/allsensors_daily/';
    %datapath = '../datafiles/NiwotSNOTELs_90-2010/';
    
    % A complete datafile has this many sensors
    fullset = 20;
    % Columns to find in daily files
    headexpr = {'Site Id' 'Date' 'Time' 'WTEQ.I-1' 'PREC.I-1' 'TOBS.I-1' ...
    'TMAX.D-1' 'TMIN.D-1' 'TAVG.D-1' 'SNWD.I-1' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20' 'RDC.I-1:-2 ' ...
    'RDC.I-1:-8' 'RDC.I-1:-20' 'BATT.I-1'};

elseif strcmp(interval, 'hourly')
    datapath = '../datafiles/soilsensors_hourly/';
    %datapath = '/home/greg/data/rawdata/SNOTELdata/soilsensors_hourly/';
    
    % A complete datafile has this many sensors
    fullset = 9;
    % Columns to find in hourly files
    headexpr = {'Site Id' 'Date' 'Time' 'SMS.I-1:-2 ' 'SMS.I-1:-8'...
    'SMS.I-1:-20' 'STO.I-1:-2 ' 'STO.I-1:-8' 'STO.I-1:-20'};
else
    error('Not a valid data type (daily or hourly)')
end

% Set dummy flag on if site_id = 0
% This will trigger dummy generator at end of file
if site_id == 0
    site_id = 828;
    dummysite = 1;
else
    dummysite = 0;
end

%--------------------------------------------------------------------------
% create a cellarray of filenames for the site
try
    filelist = textscan(ls([datapath num2str(site_id) '*']), '%s'); 
catch exception
    error('missing:SiteData', ['no data found for site ' num2str(site_id)]); 
end
files = sort(filelist{1});

%--------------------------------------------------------------------------
% Sensors, and therefore columns vary from year to year, site to site,
% therefore we must read header of each file and select desired columns.
%
% These data arrays must then be put in a standardized cell array format

% Preallocate sensor/column sorting cellarrays
t = cell(length(files), 1);
templates = {t t t t t};
[headers selectcols missingcols sensorcols colorder] = templates{:};

for i = 1:length(files)
    % read file headers to cell of strings - (1 per year of data)
    fid = fopen(files{i});
    fgetl(fid);
    headers(i) = textscan(fgetl(fid), '%s', 'Delimiter', ',');
    fclose(fid);
    % test the header file for desired header expression and create
    % logical array
    for j = 1:length(headexpr)
        n = length(headexpr{j});
        % Build a logical array showing missing sensors (rowsum = 0)
        % and sensor column location (column index where rowsum = 1)
        test = strncmpi(headexpr{j}, headers{i}', n);
        % append the array into the i-th row in selectcols{i}
        selectcols{i} = [selectcols{i}; test];
        % Ordered column locations of desired sensors in data file
        sensorcols{i} = [sensorcols{i}, find(test)];
    end
    % Logical arrays of missing sensors (1 per file)
    missingcols{i} = logical(sum(selectcols{i}, 2)');
    % Ordered column locations of sensors for appending to data files
    colorder{i} = [colorder{i}; find(missingcols{i})];
    % Change missingcols to only contain missing sensors
    missingcols{i} = find(~missingcols{i});
end

clear filelist templates selectcols i j  test; 


% -------------------------------------------------------------------------
% Create the cell array to be returned
incells = cell(1, fullset);
m = cell(1, fullset);

% read in data
for i = 1:length(files)
    fid = fopen(files{i});
    
    % Create format string to read in original file
    form = num2str(ones(1, length(headers{i})));
    forma = regexprep(form, '  1', '%f');
    format = ['%d%s%s' forma(6:end)];
    
    % Load data file into raw cell array
    rawcells = textscan(fid, format, 'Headerlines', 2, 'Delimiter', ',');
    fclose(fid);
    
    % Read data from rawcells to correct columns in incells...
    loopindex = 1;
    for j = colorder{i}
        % remember to index sensorcols by loop iteration (with loopindex)
        incells{j} = rawcells{sensorcols{i}(loopindex)};
        
        % ...and then append it to m (output array)
        m{j} = [m{j}; incells{j}];
        
        % increment loopindex
        loopindex = loopindex + 1;
    end
    clear rawcells loopindex j test form forma format;
    
    len = length(incells{1});
    
    % Insert nans for missing sensor columns
    for j = missingcols{i}
        incells{j} = nan * zeros(len, 1);
        % ...and then append it to m (output array)
        m{j} = [m{j}; incells{j}];
    end
    clear incells;
end

%--------------------------------------------------------------------------
% CLEAN OUT BAD DATA
%
% But first, extract test column to preserve original data
testColumn = 7; % m{7} is Ts @ 5cm
unfilteredData = m{testColumn};

% Start with hourly files-----
if strcmp(interval, 'hourly')
    
    % Mark and remove sub-hourly rows
    tvec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    datenumUnfiltered = datenum(tvec); % again, for bad data removal tests
    subhourly = tvec(:, 5) ~= 0;
    tvec(subhourly,:) = [];
    clear test;
    
    % remove rows marked in subhourly array
    for i = 1:fullset
        m{i}(subhourly) = [];
    end
    
    % Create a quality control array that collects errors for all columns
    qualitycontrol = false(length(m{1}), 1);
    
    for i = 4:fullset
        
        % mark the datalogger error signals (-99.9, -273.2)
        test = m{i} == -99.9 | m{i} == -273.2;
        qualitycontrol = qualitycontrol | test;
        m{i}(test) = nan;
        clear test;
        
        % change unreasonably high and low numbers to nan
        test = m{i} < -15 | m{i} > 100;
        m{i}(test) = nan;
        clear test;
        
        % Dont let soil moisture values go down to 0.0 in winter 
        % (its usually an error)
        if i<7
            test_val = m{i} == 0.0;
            test_winter = tvec(:,2) > 10 | tvec(:, 2) < 4;
            test = test_val & test_winter;
            %qualitycontrol = qualitycontrol | test;
            m{i}(test) = nan;
            clear test_val test_winter test;
        end
        
        % Remove bad soil temperature values
        if i>6           
            % Get rid of silly values
            test = m{i} > 45 | m{i} < -15;
            m{i}(test) = nan;
            clear test
        end
    end
    clear i;
    
    % change data to nan if flagged in qc array
%     for i = 4:fullset
%         m{i}(qualitycontrol) = nan;
%     end
%     clear i;

% Filter daily files-----
elseif strcmp(interval, 'daily')
    
    % Mark and remove sub-daily rows
    tvec = datevec(m{2}, 'yyyy-mm-dd');
    datenumUnfiltered = datenum(tvec); % again, for bad data removal tests
    subdaily = ~strcmp(m{3}, '');
    tvec(subdaily,:) = [];
    clear test;
    
    % remove rows marked in subdaily array
    for i = 1:fullset
        m{i}(subdaily) = [];
    end
    
    % Create a qc array that collects errors for all columns
    qualitycontrol = false(length(m{1}), 1);
    
    for i = 4:fullset
        % mark the datalogger error signals (-99.9, -273.2)
        test = m{i} == -99.9 | m{i} == -273.2;
        m{i}(test) = nan;
        qualitycontrol = qualitycontrol | test;
        clear test;
        
        % Remove negative WTEQ and PREC values
        if i<6
            test = m{i} < 0.0;
            m{i}(test) = nan;
            clear test;
        end
        
        % Remove negative snow depth values
        if i>9 && i<11
            test = m{i} < 0.0;
            m{i}(test) = nan;
            clear test;
        end
        
        % Remove High temp values for deep soil sensors
        if i>13 && i<17
            test = m{i} > 30.0;
            m{i}(test) = nan;
            clear test;
        end
        
    end
    clear i;
    
    % change data to nan if flagged in qc array
    for i = 4:fullset
        m{i}(qualitycontrol) = nan;
    end
end

%--------------------------------------------------------------------------
% CREATE DUMMYSITE if called for -  for QA testing of code
%
if strcmp(interval, 'hourly') && dummysite == 1;
    % find all July time periods
    test = tvec(:,2) == 7;
    % get the actual T data for this site (2 in)
    tmp = m{7};
    % modify the July periods to be exactly 4 degress
    tmp(test,:) = 4;
    % replace T data in m with dummy data
    m{7} = tmp;
end

%--------------------------------------------------------------------------
% PLOT the original and cleaned data series (including dummysite)
%
% h = figure;
% set(h, 'Name', ['Site ' num2str(site_id) ' - Bad data removed by loadsnotel.m']);
% plot(datenumUnfiltered, unfilteredData, '.r');
% hold on
% plot(datenum(tvec), m{testColumn}, '.k');

junk = 99;