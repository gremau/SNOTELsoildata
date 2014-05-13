function testsensors()
% testsensors.m
%
% Open all SNOTEL data files in a directory (hourly or daily), creates
% plots of each sensor timeseries, and waits for user input. User input is 
% appended to a file (excludesensors.txt) that can be used by other
% scripts to remove bad sensor data during data analysis.
%
%  * Input: all SNOTEL .csv files in a directory
%  * Output: Sensor timeseries plots, and an excludesensors.txt file
%
% 9-20-2012, Greg Maurer
% ------------------------------------------------------------------------
interval = input('Hourly or daily files? ', 's');
type = 'sigma'; % type of filter to use
thresh = 2.5; % filtering threshold

fignum = 0;
% Set parameters for daily files (all sensors) or hourly soil sensor files
if strcmpi(interval, 'daily')
    windowsize = 11;
    datapath = '../processed_data/';
    outfile = 'excludesensors_daily.txt';
    addcols = 13;
    colheaders = ['site_id, year, wteq, prec, snwd, tobs, tmax, tmin, ' ...
        'tavg, ts5, ts20, ts50, vwc5, vwc20, vwc50'];
    filelist = sortrows(csvread([datapath 'filelist_daily.txt']));
elseif strcmpi(interval, 'hourly')
    windowsize = 25;
    datapath = '../processed_data/';
    outfile = 'excludesensors_hourly.txt';
    addcols = 6;
    colheaders = 'site_id, year, ts5, ts20, ts50, vwc5, vwc20, vwc50';
    filelist = sortrows(csvread([datapath 'filelist_hourly.txt']));
end

% If there is no excludesensors.txt file, it should be created.
% If it exists, check it against filelist to see if new files
% have been added and need to be checked.
try
    % Load the sites that have been checked from excludesensor.txt
    checked = csvread([datapath outfile], 2, 0);
    notchecked = ~ismember(filelist, checked(:,1:2), 'rows');
catch exception
    disp(['File: excludesensors.txt not found. '...
        'A new file will be created']);
    checked = [];
    notchecked = ~ismember(filelist, [0,0], 'rows');
end

% Mark unchecked sites in filelist and make a list
newfiles = filelist(notchecked, 1:end);
% Add columns to the newfiles list for excluding sensors
newfiles = [newfiles, zeros(length(newfiles), addcols)];

% Create the full checklist and sort
checklist = [checked; newfiles];
checklist = sortrows(checklist, [1,2]);

% Now, if there are new files, check them first
if ~isempty(newfiles)
    disp('There are new files to check!');
    plotfilelist(newfiles);
    writechecklist;
end
% If not, review the others
uinp = input('Continue reviewing other files? (y/n)', 's');
if strcmpi(uinp, 'y')
    plotfilelist(checklist);
    writechecklist;
end;

% Function to plot a list of files, get user input, and output an exclude
% file
    function plotfilelist(flist)
        % Loop through each file in flist, plot the sensor data 
        % and ask the user to mark sensors that should be excluded. 
        % The excluded sensors are written to the checklist.
        counter = 0;
        siteID_last = NaN;
        for i = 1:length(flist)
            counter = counter + 1;
            siteID = flist(i,1);
            year = flist(i,2);
            % Get current exclusion matrix for the file
            excluded = flist(i,3:end);
            checklisttest = checklist(:,1)==siteID & checklist(:,2)==year;
            m1 = loadsnotel_oneyear(siteID, year, interval);
            % Load raw and excluded multiyear data (when siteID changes)
            if siteID~=siteID_last
                m_raw = loadsnotel(siteID, interval);
                m_exc = loadsnotel(siteID, interval, 'exclude');
            end
            
            if strcmpi(interval, 'hourly')
                disp(['Hourly data, site = ' num2str(siteID) ...
                    ', year = ' num2str(year)]);
                % Create plots to check data, return modified (or not)
                % exclusion matrix
                excluded_r = plothourly(excluded, m1, m_raw, m_exc,...
                    siteID, year);
                % Add bad sensors to exclude matrix
                checklist(checklisttest,3:end) = excluded_r;
                
            elseif strcmpi(interval, 'daily')
                disp(['Daily data, site = ' num2str(siteID) ...
                    ', year = ' num2str(year)]);
                % Create plots to check data, return modified (or not)
                % exclusion matrix
                excluded_r = plotdaily(excluded, m1, m_raw, m_exc,...
                    siteID, year);
                % Add bad sensors to exclude matrix
                checklist(checklisttest,3:end) = excluded_r;
            end
            disp(['Writing ' num2str(excluded_r) ' to excludesensors matrix']);
            disp('OK!');
            if counter==5
                writechecklist;
                counter = 0;
            end
            siteID_last = siteID;
        end
        clear m m_raw m_exc excluded
    end

% Function to write the checklist to file
    function writechecklist()
        disp(['Writing to ' datapath outfile]);
        headertext = {['Last updated in testsensors.m script on ' date];...
            colheaders};
        dlmwrite([datapath outfile], char(headertext), 'delimiter', '');
        dlmwrite([datapath outfile], checklist, '-append');
        disp('OK!');
    end

% Function to filter the data
    function f_ser = filt(ser)
        f_ser = filterseries(ser, type, windowsize, thresh);
    end

% Function to ask the user if they want to exclude a sensor.
    function r_excl = getcheck(excl)
        r_excl = excl;
        for i=1:length(excl)
            if excl(i)==0
                current = 'OK';
            elseif excl(i)==1
                current = 'Excluded';
            end
            uinp = input(['Sensor ' num2str(i) ' currently ' current...
                '. Exclude? (y/n):'], 's');
            if strcmpi(uinp, 'y')
                r_excl(i) = 1;
            elseif strcmpi(uinp, 'n');
                r_excl(i) = 0;
            end
        end
    end

% Plot HOURLY sensors
    function excl_m = plothourly(excl, m, mraw, mexc, siteID, year)
        decday_h = datenum(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
        decday_raw = datenum(strcat(mraw{2}, mraw{3}), 'yyyy-mm-ddHH:MM');
        decday_exc = datenum(strcat(mexc{2}, mexc{3}), 'yyyy-mm-ddHH:MM');
        ts5 = 7; % column 7 is at -2 in (5cm depth)
        ts20 = 8; % column 8 is at -8 in (20cm depth)
        ts50 = 9; % column 9 is at -20 in (50cm depth)
        vwc5 = 4; % column 4 is at -2 in (5cm depth)
        vwc20 = 5; % column 5 is at -8 in (20cm depth)
        vwc50 = 6; % column 6 is at -20 in (50cm depth)
        %ts5F = filterseries(m{7}, type, windowsize, thresh);
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -10;
        ymax = 45;
        set(h, 'Name', ['Ts at site ' num2str(siteID) ' - ' num2str(year)]);
        sp1 = subplot(3, 1, 1);
        plot(decday_raw, mraw{ts5}, '.r');
        hold on;
        plot(decday_exc, mexc{ts5}, '.k', decday_h, filt(m{ts5}), '.b');
        title('1) Ts -5cm'); datetick(); ylim([ymin, ymax]);
        sp2 = subplot(3, 1, 2);
        plot(decday_raw, mraw{ts20}, '.r');
        hold on;
        plot(decday_exc, mexc{ts20}, '.k', decday_h, filt(m{ts20}), '.b');
        title('2) Ts -20cm'); datetick(); ylim([ymin, ymax]);
        sp3 = subplot(3, 1, 3);
        plot(decday_raw, mraw{ts50}, '.r');
        hold on;
        plot(decday_exc, mexc{ts50}, '.k', decday_h, filt(m{ts50}), '.b');
        title('3) Ts -50cm'); datetick(); ylim([ymin, ymax]);
        linkaxes([sp1, sp2, sp3]);
        
        badts = getcheck(excl(1:3));
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -5;
        ymax = 45;
        set(h, 'Name', ['VWC at site ' num2str(siteID) ' - ' num2str(year)]);
        sp4 = subplot(3, 1, 1);
        plot(decday_raw, mraw{vwc5}, '.r');
        hold on;
        plot(decday_exc, mexc{vwc5}, '.k', decday_h, filt(m{vwc5}), '.b');
        title('1) VWC -5cm'); datetick(); ylim([ymin, ymax]);
        sp5 = subplot(3, 1, 2);
        plot(decday_raw, mraw{vwc20}, '.r');
        hold on;
        plot(decday_exc, mexc{vwc20}, '.k', decday_h, filt(m{vwc20}), '.b');
        title('2) VWC -20cm'); datetick(); ylim([ymin, ymax]);
        sp6 = subplot(3, 1, 3);
        plot(decday_raw, mraw{vwc50}, '.r');
        hold on;
        plot(decday_exc, mexc{vwc50}, '.k', decday_h, filt(m{vwc50}), '.b');
        title('3) VWC -50cm'); datetick(); ylim([ymin, ymax]);
        linkaxes([sp4, sp5, sp6]);
        
        badvwc = getcheck(excl(4:6));

        close all;
        
        excl_m = [badts badvwc];
        clear excl m mraw mexc badts badvwc;
    end

% Plot DAILY sensors
    function excl_m = plotdaily(excl, m, mraw, mexc, siteID, year)
        decday_d = datenum(m{2}, 'yyyy-mm-dd');
        decday_raw = datenum(mraw{2}, 'yyyy-mm-dd');
        decday_exc = datenum(mexc{2}, 'yyyy-mm-dd');
        wteq = 4; % water equivalent
        prec = 5; % precip
        snwd = 10; % snow depth
        tobs = 6; % Air T
        tmax = 7; % Max AirT
        tmin = 8; % Min AirT
        tavg = 9; % Avg AirT
        vwc5 = 11; % 5cm VWC
        vwc20 = 12; % 20cm VWC
        vwc50 = 13; % 50cm VWC
        ts5 = 14; % 5cm Ts
        ts20 = 15; % 20cm Ts
        ts50 = 16; % 50cm Ts

        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = 0;
        ymax = 45;
        set(h, 'Name', ['Precip at site ' num2str(siteID) ' - ' num2str(year)]);
        sp1 = subplot(3, 1, 1)
        plot(decday_raw, mraw{wteq}, '.r');
        hold on
        plot(decday_exc, mexc{wteq}, '.k', decday_d, filt(m{wteq}), '.b');
        title('1) SWE (snow pillow)'); datetick(); %ylim([ymin, ymax]);
        sp2 = subplot(3, 1, 2)
        plot(decday_raw, mraw{prec}, '.r');
        hold on
        plot(decday_exc, mexc{prec}, '.k', decday_d, filt(m{prec}), '.b');
        title('2) Precip gauge'); datetick(); %ylim([ymin, 55]);
        sp3 = subplot(3, 1, 3)
        plot(decday_raw, mraw{snwd}, '.r');
        hold on
        plot(decday_exc, mexc{snwd}, '.k', decday_d, filt(m{snwd}), '.b');
        title('3) Snow Depth'); datetick(); %ylim([ymin, 100]);
        linkaxes([sp1, sp2, sp3]);
        
        badprec = getcheck(excl(1:3));
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -25;
        ymax = 50;
        set(h, 'Name', ['AirT at site ' num2str(siteID) ' - ' num2str(year)]);
        sp4 = subplot(4, 1, 1)
        plot(decday_raw, mraw{tobs}, '.r');
        hold on
        plot(decday_exc, mexc{tobs}, '.k', decday_d, filt(m{tobs}), '.b');
        title('1) Air T obs'); datetick(); %ylim([ymin, ymax]);
        sp5 = subplot(4, 1, 2)
        plot(decday_raw, mraw{tmax}, '.r');
        hold on
        plot(decday_exc, mexc{tmax}, '.k', decday_d, filt(m{tmax}), '.b');
        title('2) MAX Air T'); datetick(); %ylim([ymin, ymax]);
        sp6 = subplot(4, 1, 3)
        plot(decday_raw, mraw{tmin}, '.r');
        hold on
        plot(decday_exc, mexc{tmin}, '.k', decday_d, filt(m{tmin}), '.b');
        title('3) MIN Air T'); datetick(); %ylim([ymin, ymax]);
        sp7 = subplot(4, 1, 4)
        plot(decday_raw, mraw{tavg}, '.r');
        hold on
        plot(decday_exc, mexc{tavg}, '.k', decday_d, filt(m{tavg}), '.b');
        title('4) AVG Air T'); datetick(); %ylim([ymin, ymax]);
        linkaxes([sp4, sp5, sp6, sp7]);
        
        badairt = getcheck(excl(4:7));
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -10;
        ymax = 45;
        set(h, 'Name', ['Ts at site ' num2str(siteID) ' - ' num2str(year)]);
        sp8 = subplot(3, 1, 1)
        plot(decday_raw, mraw{ts5}, '.r');
        hold on
        plot(decday_exc, mexc{ts5}, '.k', decday_d, filt(m{ts5}), '.b');
        title('1) Ts -5cm'); datetick(); ylim([ymin, ymax]);
        sp9 = subplot(3, 1, 2)
        plot(decday_raw, mraw{ts20}, '.r');
        hold on
        plot(decday_exc, mexc{ts20}, '.k', decday_d, filt(m{ts20}), '.b');
        title('2) Ts -20cm'); datetick(); ylim([ymin, ymax]);
        sp10 = subplot(3, 1, 3)
        plot(decday_raw, mraw{ts50}, '.r');
        hold on
        plot(decday_exc, mexc{ts50}, '.k', decday_d, filt(m{ts50}), '.b');
        title('3) Ts -50cm'); datetick(); ylim([ymin, ymax]);
        linkaxes([sp8, sp9, sp10]);
        
        badts = getcheck(excl(8:10));
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -5;
        ymax = 45;
        set(h, 'Name', ['VWC at site ' num2str(siteID) ' - ' num2str(year)]);
        sp11 = subplot(3, 1, 1)
        plot(decday_raw, mraw{vwc5}, '.r');
        hold on
        plot(decday_exc, mexc{vwc5}, '.k', decday_d, filt(m{vwc5}), '.b');
        title('1) VWC -5cm'); datetick(); ylim([ymin, ymax]);
        sp12 = subplot(3, 1, 2)
        plot(decday_raw, mraw{vwc20}, '.r');
        hold on
        plot(decday_exc, mexc{vwc20}, '.k', decday_d, filt(m{vwc20}), '.b');
        title('2) VWC -20cm'); datetick(); ylim([ymin, ymax]);
        sp13 = subplot(3, 1, 3)
        plot(decday_raw, mraw{vwc50}, '.r');
        hold on
        plot(decday_exc, mexc{vwc50}, '.k', decday_d, filt(m{vwc50}), '.b');
        title('3) VWC -50cm'); datetick(); ylim([ymin, ymax]);
        linkaxes([sp11, sp12, sp13]);
        
        badvwc = getcheck(excl(11:13));
        
        close all;
        
        excl_m = [badprec badairt badts badvwc];
        clear excl m mraw mexc badprec badairt badts badvwc;
    end
end
