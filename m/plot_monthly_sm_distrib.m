% plot_monthly_sm_distrib.m
%
% Plots distributions of soil moisture for individual months and years. 
% Data for all years (available) is represented.
%
% started 120119 GM


close all; clear all;   % clear figures and variables in the workspace
fignum = 0;     % used to increment figure number for plots
%addpath('../m/');
addpath('/home/greg/data/programming/m_common/');

% Ask user for site number
siteID = str2num(input('Which SNOTEL station?: ', 's'));

% load hourly data from site  w/ loadsnotel:
[siteHourly, headers] = loadsnotel('hourly', siteID);

% parse out the date, soil moisture, and soil temperature sensors
datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
sm = siteHourly{4}; % 5cm depth
%sm = siteHourly{5};  % 20cm depth
%sm = siteHourly{6};  % 60cm depth

% create a water year vector that matches size of datevec_h year vector
wyearvec = datevec_h(:, 1);
wytest = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
wyearvec(wytest) = wyearvec(wytest) + 1;

% list of water years in the dataset
wyears = unique(wyearvec);

% Loop through each water year in dataset and select data, then plot
for i = 1:length(wyears)
    yeartest = wyearvec(:) == wyears(i); % Use water years to select data
    datevec_sel = datevec_h(yeartest, :);
    sm_sel = sm(yeartest, :);

     % this sorta works for wateryear
%     wyeartest = circshift(yeartest, -92*24);
%     datevec_sel = datevec_h(wyeartest, :);
%     sm_sel = sm(wyeartest, :);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_sel(:,2)==10 | datevec_sel(:,2)==11 | datevec_sel(:,2)==12);
    testJFM = (datevec_sel(:,2)==1 | datevec_sel(:,2)==2 | datevec_sel(:,2)==3);
    testAMJ = (datevec_sel(:,2)==4 | datevec_sel(:,2)==5 | datevec_sel(:,2)==6);
    testJAS = (datevec_sel(:,2)==7 | datevec_sel(:,2)==8 | datevec_sel(:,2)==9);
    sm_OND = sm_sel(testOND);
    sm_JFM = sm_sel(testJFM);
    sm_AMJ = sm_sel(testAMJ);
    sm_JAS = sm_sel(testJAS);
    
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
    % PLOT the data
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - SM histograms - WY'...
        num2str(wyears(i))]);
    
    % Four plots on the left are for quarters
    subplot (4, 4, 1)
    xedges = 0:1:100;
    n = histc(sm_OND, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'k');
    hold on
    plot([sm_ONDm sm_ONDm], [0 1], ':k');
    axis([0 75 0 0.35]);
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
    axis([0 75 0 0.35]);
    ylabel('Frequency');
    title('Jan 1 -Mar 31');
    
    subplot (4, 4, 9)
    xedges = 0:1:100;
    n = histc(sm_AMJ, xedges);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'r');
    hold on
    plot([sm_AMJm sm_AMJm], [0 1], ':k');
    axis([0 75 0 0.35]);
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
    axis([0 75 0 0.35]);
    ylabel('Frequency');
    title('July 1 - Sept 30');
    
    
    % Plot the month distributions in 3 rows 
    plotorder = [6 7 8 10 11 12 14 15 16 2 3 4];
    
    for i = 1:length(plotorder)
        subplot(4, 4, plotorder(i));
        xedges = 0:1:100;
        eval(['n = histc(sm_' months{i} ', xedges);']);
        % normalize
        nnorm = n./sum(n);
        bar (xedges, nnorm, 'g');
        axis([0 75 0 0.35]);
        title(months{i});
    end
end

junk = 99;
