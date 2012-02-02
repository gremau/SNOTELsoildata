% plot_monthly_sm_distrib.m
%
% Plots distributions of soil moisture for individual months and years. 
% Data for all years (available) is represented.
%
% started 120119 GM

close all; clear all;   % clear figures and variables in the workspace
fignum = 0;     % used to increment figure number for plots
addpath('/home/greg/data/programming/m_common/');

% Ask user for site number and depth
siteID = str2num(input('Which SNOTEL station?: ', 's'));
sensordepth = str2num(input(...
    'What sensor for histograms (1 = 5cm,2 = 20cm, 3 = 50cm)?: ', 's'));

% Load hourly data from site  w/ loadsnotel:
siteHourly = loadsnotel('hourly', siteID);

% Parse out the date and soil moisture sensors and set normalizing/axes
datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
datenum_h = datenum(datevec_h);

% Run RAW SENSOR DATA
% sm = [siteHourly{4}, siteHourly{5}, siteHourly{6}]; % 5cm depth no norm
% xax = 75
% Run NORMALIZED data with smnormalize
sm = [smnormalize(siteHourly{4}, 1), smnormalize(siteHourly{5}, 1), ...
    smnormalize(siteHourly{6}, 1)].*100;
xax = 100; % these axes are good for normalized data
yax = 0.2;
disp('*** Running in normalized soil moisture data mode ***');

% PLOT 1 - the entire timeseries - all depths
titles = {'5cm depth' '20cm depth' '50cm depth'};
ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - SM timeseries - All data']);
for i = 1:size(sm, 2)
    subplot (3, 1, i)
    plot(datenum_h, sm(:, i), 'k');
    set(gca, 'XTick', ticklocations);
    set(gca, 'XTickLabel', ticklocations);
    datetick('x', 12, 'keepticks');
    ylabel('Normalized soil moisture');
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
        yax = 0.1;
    end
    datevec_sel = datevec_h(yeartest, :);
    sm_sel = sm(yeartest, sensordepth);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_sel(:,2)==10 | datevec_sel(:,2)==11 | datevec_sel(:,2)==12);
    testJFM = (datevec_sel(:,2)==1 | datevec_sel(:,2)==2 | datevec_sel(:,2)==3);
    testAMJ = (datevec_sel(:,2)==4 | datevec_sel(:,2)==5 | datevec_sel(:,2)==6);
    testJAS = (datevec_sel(:,2)==7 | datevec_sel(:,2)==8 | datevec_sel(:,2)==9);
    sm_OND = sm_sel(testOND,:);
    sm_JFM = sm_sel(testJFM,:);
    sm_AMJ = sm_sel(testAMJ,:);
    sm_JAS = sm_sel(testJAS,:);
    
    % Calculate mean and standard deviations of this data
    sm_ONDm = mean(sm_OND(~isnan(sm_OND)));
    sm_JFMm = mean(sm_JFM(~isnan(sm_JFM)));
    sm_AMJm = mean(sm_AMJ(~isnan(sm_AMJ)));
    sm_JASm = mean(sm_JAS(~isnan(sm_JAS)));
    sm_ONDsd = std(sm_OND(~isnan(sm_OND)));
    sm_JFMsd = std(sm_JFM(~isnan(sm_JFM)));
    sm_AMJsd = std(sm_AMJ(~isnan(sm_AMJ)));
    sm_JASsd = std(sm_JAS(~isnan(sm_JAS)));
    
    % Logical tests, means and std deviations for months
    months = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
    for j=1:length(months)
        eval(['test' months{j} ' = datevec_sel(:,2)==' num2str(j) ';']);
        eval(['sm_' months{j} ' = sm_sel(test' months{j} ');']);
        eval(['sm_' months{j} 'm = mean(sm_' months{j} '(~isnan(sm_' months{j} ')));']);
        eval(['sm_' months{j} 'sd = std(sm_' months{j} '(~isnan(sm_' months{j} ')));']);
    end
    
    clear testOND testJFM testAMJ testJAS;
    
    %-------------------------------------------------------------
    % PLOTS 2 -> end. One per water year and one for all years
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - SM histograms - WY'...
        num2str(wyears(i))]);
    
    % First four subplots on the left are for quarters
    subplot (4, 4, 1)
    xedges = 0:1:100;
    n = histc(sm_OND, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'y');
    hold on
    plot([sm_ONDm sm_ONDm], [0 1], ':k');
    axis([0 xax 0 yax]);
    ylabel('Frequency');
    title('Oct 1-Dec 31');
    
    subplot (4, 4, 5)
    xedges = 0:1:100;
    n = histc(sm_JFM, xedges);
    %normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'b');
    hold on
    plot([sm_JFMm sm_JFMm], [0 1], ':k');
    axis([0 xax 0 yax]);
    ylabel('Frequency');
    title('Jan 1 -Mar 31');
    
    subplot (4, 4, 9)
    xedges = 0:1:100;
    n = histc(sm_AMJ, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'y');
    hold on
    plot([sm_AMJm sm_AMJm], [0 1], ':k');
    axis([0 xax 0 yax]);
    ylabel('Frequency');
    title('April 1 - June 30');
    
   subplot (4, 4, 13)
    xedges = 0:1:100;
    n = histc(sm_JAS, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'r');
    hold on
    plot([sm_JASm sm_JASm], [0 1], ':k');
    axis([0 xax 0 yax]);
    ylabel('Frequency');
    title('July 1 - Sept 30');
    
    
    % Subplots of the monthly distributions in 3 rows 
    plotorder = [6 7 8 10 11 12 14 15 16 2 3 4];
    
    for j = 1:length(plotorder)
        subplot(4, 4, plotorder(j));
        xedges = 0:1:100;
        eval(['n = histc(sm_' months{j} ', xedges);']);
        % normalize
        nnorm = n./sum(n);
        bar (xedges, nnorm, 'g');
        axis([0 xax 0 0.25]);
        title(months{j});
    end
end

junk = 99;
