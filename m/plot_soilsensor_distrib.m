% plot_soilsensor_distrib.m
%
% Plots distributions of SNOTEL soil sensor data for individual months 
% and years. Data for all years (available) is represented.
%
% started 120119 GM

close all; clear all;   % clear figures and variables in the workspace
fignum = 0;     % used to increment figure number for plots
addpath('/home/greg/data/programming_resources/m_common/');

% Ask user for site number and depth
siteID = str2num(input('Which SNOTEL station?: ', 's'));
sensoroutput = input('Which sensor output (vwc or temp)?: ', 's');
sensordepth = str2num(input(...
    'What sensor depth in histograms (1=5cm, 2=20cm, 3=50cm)?: ', 's'));

% Load hourly data from site  w/ loadsnotel:
siteHourly = loadsnotel(siteID, 'hourly', 'exclude');

% Filter the data with filterseries.m
for i=4:9
    siteHourly{i} = filterseries(siteHourly{i}, 'sigma', 25, 3);
end

% Parse out the date and soil moisture sensors and set normalizing/axes
datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
datenum_h = datenum(datevec_h);

% Select TEMP or VWC data and set distribution bins and plot axes
if strcmpi(sensoroutput, 'vwc');
    % Run RAW SENSOR DATA - no normalization
    % sensors = [siteHourly{4}, siteHourly{5}, siteHourly{6}];
    % xedges = 0:1:100; % raw sm data bins (0-100)
    % xmax = 75
    % Run NORMALIZED data with smnormalize
    sensors = [smnormalize(siteHourly{4}, 1), smnormalize(siteHourly{5}, 1), ...
        smnormalize(siteHourly{6}, 1)];
    xedges = 0:0.01:1; % normalized vwc bins (0-1)
    xmax = 1; % these axes are good for normalized data
    ymax = 0.2;
    disp('*** Running in normalized soil moisture data mode ***');
    xmin = 0;
elseif strcmpi(sensoroutput, 'temp');
    % Load soil temperature data(no normalization)
    sensors = [siteHourly{7}, siteHourly{8}, siteHourly{9}];
    xedges = -15:0.5:35; % soil temp bins -15-35 deg C
    xmax = 35; % corresponding x axis limits
    xmin = -15;
    ymax = 0.2;
end

% PLOT 1 - the entire timeseries - all depths
titles = {'5cm depth' '20cm depth' '50cm depth'};
ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Full sensor timeseries']);
for i = 1:size(sensors, 2)
    subplot (3, 1, i)
    plot(datenum_h, sensors(:, i), 'k');
    set(gca, 'XTick', ticklocations);
    set(gca, 'XTickLabel', ticklocations);
    datetick('x', 12, 'keepticks');
    ylabel(sensoroutput);
    title(titles{i});
end

% Create a water year vector that matches size of datevec_h year vector
wyearvec = datevec_h(:, 1);
wytest = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
wyearvec(wytest) = wyearvec(wytest) + 1;
% List of water years in the dataset
wyears = [unique(wyearvec); 0]; % 0 sets the all years option in loop

% Loop through each water year in dataset and select data, then plot
for i = 1:length(wyears)
    if wyears(i) > 0
        yeartest = wyearvec(:) == wyears(i); % Use water years to select data
    elseif wyears(i)==0
        yeartest = wyearvec(:) > 0;
        ymax = 0.1;
    end
    datevec_sel = datevec_h(yeartest, :);
    sensors_sel = sensors(yeartest, sensordepth);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_sel(:,2)==10 | datevec_sel(:,2)==11 | datevec_sel(:,2)==12);
    testJFM = (datevec_sel(:,2)==1 | datevec_sel(:,2)==2 | datevec_sel(:,2)==3);
    testAMJ = (datevec_sel(:,2)==4 | datevec_sel(:,2)==5 | datevec_sel(:,2)==6);
    testJAS = (datevec_sel(:,2)==7 | datevec_sel(:,2)==8 | datevec_sel(:,2)==9);
    sensors_OND = sensors_sel(testOND,:);
    sensors_JFM = sensors_sel(testJFM,:);
    sensors_AMJ = sensors_sel(testAMJ,:);
    sensors_JAS = sensors_sel(testJAS,:);
    
    % Calculate mean and standard deviations of this data
    sensors_ONDm = mean(sensors_OND(~isnan(sensors_OND)));
    sensors_JFMm = mean(sensors_JFM(~isnan(sensors_JFM)));
    sensors_AMJm = mean(sensors_AMJ(~isnan(sensors_AMJ)));
    sensors_JASm = mean(sensors_JAS(~isnan(sensors_JAS)));
    sensors_ONDsd = std(sensors_OND(~isnan(sensors_OND)));
    sensors_JFMsd = std(sensors_JFM(~isnan(sensors_JFM)));
    sensors_AMJsd = std(sensors_AMJ(~isnan(sensors_AMJ)));
    sensors_JASsd = std(sensors_JAS(~isnan(sensors_JAS)));
    
    % Logical tests, means and std deviations for months
    months = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
    for j=1:length(months)
        eval(['test' months{j} ' = datevec_sel(:,2)==' num2str(j) ';']);
        eval(['sensors_' months{j} ' = sensors_sel(test' months{j} ');']);
        eval(['sensors_' months{j} 'm = mean(sensors_' months{j}...
            '(~isnan(sensors_' months{j} ')));']);
        eval(['sensors_' months{j} 'sd = std(sensors_' months{j}...
            '(~isnan(sensors_' months{j} ')));']);
    end
    
    clear testOND testJFM testAMJ testJAS;
    
    %-------------------------------------------------------------
    % PLOTS 2 -> end. One per water year and one for all years
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' sensoroutput ...
        ' histograms - WY' num2str(wyears(i))]);
    
    % First four subplots on the left are for quarters
    subplot (4, 4, 1)
    n = histc(sensors_OND, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'y');
    hold on
    plot([sensors_ONDm sensors_ONDm], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title('Oct 1-Dec 31');
    
    subplot (4, 4, 5)
    n = histc(sensors_JFM, xedges);
    %normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'b');
    hold on
    plot([sensors_JFMm sensors_JFMm], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title('Jan 1 -Mar 31');
    
    subplot (4, 4, 9)
    n = histc(sensors_AMJ, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'y');
    hold on
    plot([sensors_AMJm sensors_AMJm], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title('April 1 - June 30');
    
    subplot (4, 4, 13)
    n = histc(sensors_JAS, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'r');
    hold on
    plot([sensors_JASm sensors_JASm], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title('July 1 - Sept 30');
    
    % Subplots of the monthly distributions in 3 rows 
    plotorder = [6 7 8 10 11 12 14 15 16 2 3 4];
    
    for j = 1:length(plotorder)
        subplot(4, 4, plotorder(j));
        eval(['n = histc(sensors_' months{j} ', xedges);']);
        % normalize
        nnorm = n./sum(n);
        bar (xedges, nnorm, 'g');
        axis([xmin xmax 0 0.25]);
        title(months{j});
    end
end

junk = 99;
