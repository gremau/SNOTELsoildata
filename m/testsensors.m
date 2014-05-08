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
    function excl_m = plotdaily(excl, m, mraw, mexc)
        decday_d = datenum(m{2});
        wteq = m{4}; % water equivalent
        prec = m{5}; % precip
        snwd = m{10}; % snow depth
        tobs = m{6}; % Air T
        tmax = m{7}; % Max AirT
        tmin = m{8}; % Min AirT
        tavg = m{9}; % Avg AirT
        vwc5 = m{11}; % 5cm VWC
        vwc20 = m{12}; % 20cm VWC
        vwc50 = m{13}; % 50cm VWC
        ts5 = m{14}; % 5cm Ts
        ts20 = m{15}; % 20cm Ts
        ts50 = m{16}; % 50cm Ts
        % Filtered
        wteqF = filterseries(m{4}, type, windowsize, thresh);
        precF = filterseries(m{5}, type, windowsize, thresh);
        snwdF = filterseries(m{10}, type, windowsize, thresh);
        tobsF = filterseries(m{6}, type, windowsize, thresh);
        tmaxF = filterseries(m{7}, type, windowsize, thresh);
        tminF = filterseries(m{8}, type, windowsize, thresh);
        tavgF = filterseries(m{9}, type, windowsize, thresh);
        vwc5F = filterseries(m{11}, type, windowsize, thresh);
        vwc20F = filterseries(m{12}, type, windowsize, thresh);
        vwc50F = filterseries(m{13}, type, windowsize, thresh);
        ts5F = filterseries(m{14}, type, windowsize, thresh);
        ts20F = filterseries(m{15}, type, windowsize, thresh);
        ts50F = filterseries(m{16}, type, windowsize, thresh);
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = 0;
        ymax = 45;
        set(h, 'Name', ['Precip at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(3, 1, 1)
        plot(decday_d, wteq, '.r', decday_d, wteqF, '.k');
        title('1) SWE (snow pillow)'); datetick(); %ylim([ymin, ymax]);
        subplot(3, 1, 2)
        plot(decday_d, prec, '.r', decday_d, precF, '.k');
        title('2) Precip gauge'); datetick(); %ylim([ymin, 55]);
        subplot(3, 1, 3)
        plot(decday_d, snwd, '.r', decday_d, snwdF, '.k');
        title('3) Snow Depth'); datetick(); %ylim([ymin, 100]);
        precinput = input('Bad precip sensors?  [1, 2, or 3]: ');
        badprec = zeros(1, 3);
        badprec(precinput) = 1;
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -25;
        ymax = 50;
        set(h, 'Name', ['AirT at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(4, 1, 1)
        plot(decday_d, tobs, '.r', decday_d, tobsF, '.k');
        title('1) Air T obs'); datetick(); %ylim([ymin, ymax]);
        subplot(4, 1, 2)
        plot(decday_d, tmax, '.r', decday_d, tmaxF, '.k');
        title('2) MAX Air T'); datetick(); %ylim([ymin, ymax]);
        subplot(4, 1, 3)
        plot(decday_d, tmin, '.r', decday_d, tminF, '.k');
        title('3) MIN Air T'); datetick(); %ylim([ymin, ymax]);
        subplot(4, 1, 4)
        plot(decday_d, tavg, '.r', decday_d, tavgF, '.k');
        title('4) AVG Air T'); datetick(); %ylim([ymin, ymax]);
        airtinput = input('Bad air T sensors?  [1, 2, 3, or 4]: ');
        badairt = zeros(1, 4);
        badairt(airtinput) = 1;
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -10;
        ymax = 45;
        set(h, 'Name', ['Ts at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(3, 1, 1)
        plot(decday_d, ts5, '.r', decday_d, ts5F, '.k');
        title('1) Ts -5cm'); datetick(); ylim([ymin, ymax]);
        subplot(3, 1, 2)
        plot(decday_d, ts20, '.r', decday_d, ts20F, '.k');
        title('2) Ts -20cm'); datetick(); ylim([ymin, ymax]);
        subplot(3, 1, 3)
        plot(decday_d, ts50, '.r', decday_d, ts50F, '.k');
        title('3) Ts -50cm'); datetick(); ylim([ymin, ymax]);
        tsinput = input('Bad Tsoil sensors?  [1, 2, or 3]: ');
        badts = zeros(1, 3);
        badts(tsinput) = 1;
        
        fignum = fignum+1;
        h = figure(fignum);
        ymin = -5;
        ymax = 45;
        set(h, 'Name', ['VWC at site ' num2str(siteID) ' - ' num2str(year)]);
        subplot(3, 1, 1)
        plot(decday_d, vwc5, '.r', decday_d, vwc5F, '.k');
        title('1) VWC -5cm'); datetick(); ylim([ymin, ymax]);
        subplot(3, 1, 2)
        plot(decday_d, vwc20, '.r', decday_d, vwc20F, '.k');
        title('2) VWC -20cm'); datetick(); ylim([ymin, ymax]);
        subplot(3, 1, 3)
        plot(decday_d, vwc50, '.r', decday_d, vwc50F, '.k');
        title('3) VWC -50cm'); datetick(); ylim([ymin, ymax]);
        vwcinput = input('Bad VWC sensors?  [1, 2, or 3]: ');
        badvwc = zeros(1, 3);
        badvwc(vwcinput) = 1;
        close all;
        
        excl_m = [badprec badairt badts badvwc];
        clear excl m mraw mexc badprec badairt badts badvwc;
    end
end
