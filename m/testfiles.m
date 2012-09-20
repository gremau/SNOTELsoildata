function testfiles()
% testfiles.m
%
%
% arguments:
% interval = 'hourly' or 'daily' ... load hourly or daily data files
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
sitelist = sortrows(csvread([datapath 'sitelist.txt']));

% Add columns that can be used to exclude sensors
badfilelist = [sitelist, zeros(length(sitelist), 3)];
    
for i = 1:length(badfilelist)
    siteID = badfilelist(i,1);
    year = badfilelist(i,2);
    m = loadsnotel_oneyear(interval, siteID, year);
    
    if strcmpi(interval, 'hourly')

        % Check that the data period is not truncated in any way
        sampletvec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        sampletnum = datenum(sampletvec);
        % How many days are measured?
        measDays = length(unique(floor(sampletnum)));
        badfilelist(i, 3) = measDays;
%         % Check for files that are missing lots of days
%         if measDays < 350
%             disp(['Hourly data, site = ' num2str(siteID) ...
%                 ', year = ' num2str(year)]);
%             disp(['Only ' num2str(measDays) ' days measured']);
%             badfilelist(i, 3) = 1;
%         end

        % Look at how frequent measurements are
        measHours = unique(sampletvec(:,4));
        measNum = length(measHours);
        badfilelist(i, 4) = measNum;
        
        % Create a full hourly time vector for the water year
        tstart = datenum(year-1,10,1);
        tend = datenum(year,9,30,23,0,0);
        fulltvec = datevec(tstart:(1/measNum):tend);
        % Find missing measurement times 
        missingRows = find(~ismember(fulltvec, sampletvec, 'rows'));
        missingTimes = fulltvec(missingRows,:);
        % Find percentage of missing measurement times
        measPcnt = size(missingTimes,1)/length(fulltvec);
        badfilelist(i, 5) = measPcnt;


    elseif strcmpi(interval, 'daily')
        disp(['Daily data, site = ' num2str(siteID) ...
            ', year = ' num2str(year)]);
        % Check that the data period is not truncated 
        tvec = datevec(m{2}, 'yyyy-mm-dd');
        tstart = datenum(year,1,1);
        tend = datenum(year,12,31);
        fulltvec = datevec(tstart:tend)
        % Create plots to check data
        plotdaily();
    end

end