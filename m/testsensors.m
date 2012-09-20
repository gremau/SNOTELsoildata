function testsensors()
% testsensors.m
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
sitelist = sortrows(csvread([datapath 'filelist.txt']));

% Add columns that can be used to exclude sensors
sensorchecklist = [sitelist, zeros(length(sitelist), 6)];
    
for i = 1:length(sensorchecklist)
    siteID = sensorchecklist(i,1);
    year = sensorchecklist(i,2);
    m = loadsnotel_oneyear(interval, siteID, year);
    
    if strcmpi(interval, 'hourly')
        disp(['Hourly data, site = ' num2str(siteID) ...
            ', year = ' num2str(year)]);
        % Check that the data period is not truncated in any way
        sampletvec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        sampletnum = datenum(sampletvec);
        measDays = length(unique(floor(sampletnum)));
        if measDays < 350
            disp(['Only ' num2str(measDays) ' days measured']);
        end
        % Create a full hourly time vector for the water year
        tstart = datenum(year-1,10,1);
        tend = datenum(year,9,30,23,0,0);
        fulltvec = datevec(tstart:(1/24):tend);
        % Look at how frequent measurements are
        measHours = unique(sampletvec(:,4));
        if length(measHours) < 24
            disp(['Only '  num2str(length(measHours)) ...
                ' hours are measured!']);
        end
        % Find missing measurement times 
        missingRows = find(~ismember(fulltvec, sampletvec, 'rows'));
        missingTimes = fulltvec(missingRows,:);
        if size(missingTimes,1) > (0.1*length(fulltvec))
            disp([num2str(100*(size(missingTimes, 1)/length(fulltvec))) ...
                '% of measurement times are missing!']);
        end
        % Create plots to check data
        plothourly();
        % Add bad sensors to exclude matrix
        badsensors = [badtemp badvwc];
        sensorchecklist(i,3:end) = badsensors;

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
    disp('OK!');
end

    function plothourly()
        decday_h = datenum(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        Ts5 = m{7}; % column 7 is at -2 in (5cm depth)
        Ts20 = m{8}; % column 8 is at -8 in (20cm depth)
        Ts50 = m{9}; % column 9 is at -20 in (50cm depth)
        WC5 = m{4}; % column 4 is at -2 in (5cm depth)
        WC20 = m{5}; % column 5 is at -8 in (20cm depth)
        WC50 = m{6}; % column 6 is at -20 in (50cm depth)
        Ts5filt = filterseries(m{7}, 'sigma', 2.5);
        Ts20filt = filterseries(m{8}, 'sigma', 2.5);
        Ts50filt = filterseries(m{9}, 'sigma', 2.5);
        WC5filt = filterseries(m{4}, 'sigma', 2.5);
        WC20filt = filterseries(m{5}, 'sigma', 2.5);
        WC50filt = filterseries(m{6}, 'sigma', 2.5);
        
        fignum = fignum+1;
        h = figure(fignum);
        set(h, 'Name', ['Ts at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(3, 1, 1)
        plot(decday_h, Ts5, '.r', decday_h, Ts5filt, '.k');
        title('Ts -5cm'); datetick();
        subplot(3, 1, 2)
        plot(decday_h, Ts20, '.r', decday_h, Ts20filt, '.k');
        title('Ts -20cm'); datetick();
        subplot(3, 1, 3)
        plot(decday_h, Ts50, '.r', decday_h, Ts50filt, '.k');
        title('Ts - 50cm'); datetick();
        %legend('Raw data', 'Filtered using shift', 'Location', 'NorthWest');
        tempinput = input('Bad temp sensors?  [1, 2, or 3]: ');
        badtemp = zeros(1, 3);
        badtemp(tempinput) = 1;
        
        %pause();
        
        fignum = fignum+1;
        h = figure(fignum);
        set(h, 'Name', ['VWC at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(3, 1, 1)
        plot(decday_h, WC5, '.r', decday_h, WC5filt, '.k');
        title('VWC -5cm'); datetick();
        subplot(3, 1, 2)
        plot(decday_h, WC20, '.r', decday_h, WC20filt, '.k');
        title('VWC -20cm'); datetick();
        subplot(3, 1, 3)
        plot(decday_h, WC50, '.r', decday_h, WC50filt, '.k');
        title('VWC - 50cm'); datetick();
        %legend('Raw data', 'Filtered using shift', 'Location', 'NorthWest');
        vwcinput = input('Bad vwc sensors?  [1, 2, or 3]: ');
        badvwc = zeros(1, 3);
        badvwc(vwcinput) = 1;
        close all;
        
    end

end
%csvwrite([processeddatapath 'wyear_soilwatersummary_smnorm.csv'], soilwatersummary);