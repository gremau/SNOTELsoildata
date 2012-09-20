function testfiles()
% testfiles.m
%
% Tests all SNOTEL data files in a directory (hourly or daily) and 
% determines how many measurements are missing in each file. If too many 
% datapoints are missing, the file is listed in an "exlude data" file that 
% is then read by other scripts.
%
%  * Input: all SNOTEL .csv files in a directory
%  * Output: Histograms of the data, and excludefiles.txt file
%
% Greg Maurer 9/20/2012
% ------------------------------------------------------------------------
interval = input('Hourly or daily files? ', 's');

fignum = 0;
% Load a list of daily files (all sensors) or hourly soil sensor data 
if strcmpi(interval, 'daily')
    datapath = '../rawdata/allsensors_daily/';
elseif strcmpi(interval, 'hourly')
    datapath = '../rawdata/soilsensors_hourly/';
end

% Create a list of files in the datapath
filelist = sortrows(csvread([datapath 'filelist.txt']));
% Add columns that can be used to exclude sensors
checklist = [filelist, zeros(length(filelist), 2)];

% Loop through all the sites and check for missing data
for i = 1:length(checklist)
    siteID = checklist(i,1);
    year = checklist(i,2);
    m = loadsnotel_oneyear(interval, siteID, year);
    
    if strcmpi(interval, 'hourly')
        % Check that the data period is not truncated in any way
        sampletvec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        sampletnum = datenum(sampletvec);
        % How many total days have measurements?
        measDays = floor(sampletnum);
        numDays = length(unique(measDays));
        checklist(i, 3) = numDays;

        % How many hours are measured within the file?
        %measHours = unique(sampletvec(:,4));
        %measNum = length(measHours);
        %checklist(i, 5) = measNum;
        
        % Create a full hourly time vector for the water year
        tstart = datenum(year-1,10,1);
        tend = datenum(year,9,30,23,0,0);
        fulltvec = datenum(tstart:(1/24):tend);
        
        % Find missing measurement times 
        %missingRows = find(~ismember(fulltvec, sampletvec, 'rows'));
        %missingTimes = fulltvec(missingRows,:);

        % How many days have less than 6 measurements?
        possibleDays = floor(fulltvec);
        freqDays = histc(measDays, linspace(min(possibleDays), ...
            max(possibleDays), length(unique(possibleDays))));
        badDays = sum(freqDays<6);
        checklist(i,4) = badDays;

    elseif strcmpi(interval, 'daily')
        % Check that the data period is not truncated in any way
        sampletvec = datevec(m{2}, 'yyyy-mm-dd');
        sampletnum = datenum(sampletvec);
        % How many total days have measurements?
        measDays = floor(sampletnum);
        numDays = length(unique(measDays));
        checklist(i, 3) = numDays;
    end
end

% PLOT frequency histograms of the bad files
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' Histograms of missing data']);
% Plot frequency of the number of days contained in the files
xedges = linspace(0, 366, 52); % number of bins
subplot(1,2,1);
n1 = histc(checklist(:,3), xedges); % sort into bins
bar(xedges, n1, 'r'); % plot in barchart
ylabel({'number of','occurences'});
xlabel('Days/file');
% Plot frequency of low measurement days in the files
xedges = linspace(0, 366, 52); % number of bins
subplot(1,2,2);
n2 = histc(checklist(:,4), xedges);
bar(xedges, n2, 'k');
ylabel({'number of','occurences'});
xlabel('Low measurement days (<6)');

% Remove any files in which 5% of the data is missing (18 days worth).
test = checklist(:,3) < 347 | checklist(:,4) > 18;
excludefiles = checklist(test, 1:2);

% Write output to a textfile in the data directory.
csvwrite([datapath 'excludefiles.txt'], excludefiles);
junk = 99;