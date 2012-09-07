function m = test_files(interval)
% test_file.m
%
%
% arguments:
% interval = 'hourly' or 'daily' ... load hourly or daily data files
% ------------------------------------------------------------------------
fignum = 0;

% Load a list of daily files (all sensors) or hourly soil sensor data 
if strcmpi(interval, 'daily')
    datapath = '../rawdata/allsensors_daily/';
    

elseif strcmpi(interval, 'hourly')
    datapath = '../rawdata/soilsensors_hourly/';
    
end

filelist = textscan(ls([datapath '*.csv']), '%s');
filelist = sortrows(filelist{:});
sitelist = sortrows(csvread([datapath 'sitelist.txt'], 1, 0));
    
for i = 1:size(filelist)
    m = loadsnotel_oneyear(filelist{i}, interval);
    if strcmpi(interval, 'hourly')
        plothourly();
        %pause();
    elseif strcmpi(interval, 'daily')
        plotdaily();
        pause()
    end
    
end

    function null = plothourly()
        decday_h = datenum(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        Ts5 = m{7}; % column 7 is at -2 in (5cm depth)
        Ts20 = m{8}; % column 8 is at -8 in (20cm depth)
        Ts50 = m{9}; % column 9 is at -20 in (50cm depth)
        WC5 = m{4}; % column 4 is at -2 in (5cm depth)
        WC20 = m{5}; % column 5 is at -8 in (20cm depth)
        WC50 = m{6}; % column 6 is at -20 in (50cm depth)
        Ts5filt = filtertempseries(m{7}, 'shift', 2.5);
        Ts20filt = filtertempseries(m{8}, 'shift', 2.5);
        Ts50filt = filtertempseries(m{9}, 'shift', 2.5);
        WC5filt = filtertempseries(m{4}, 'shift', 2.5);
        WC20filt = filtertempseries(m{5}, 'shift', 2.5);
        WC50filt = filtertempseries(m{6}, 'shift', 2.5);
        
        fignum = fignum+1;
        h = figure(fignum);
        set(h, 'Name', [filelist{i} ' - Hourly Ts']);
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
        pause();
        
        fignum = fignum+1;
        h = figure(fignum);
        set(h, 'Name', [filelist{i} ' - Hourly VWC']);
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
        pause();
        close all;
        
    end
    
end