% plot_monthly_sm_distrib.m
%
% Plots distributions of soil temperature for winter time periods and 
% individual months. Data for all years (available) is represented.
%
% ver 2: 111116 GM
% was soildistributions.m
% Changed to script with user input

close all;      % clear any figures
fignum = 0;     % used to increment figure number for plots
%addpath('../m/');
addpath('/home/greg/data/programming/m_common/');

% Ask user for site number
siteID = str2num(input('Which SNOTEL station?: ', 's'));

% load hourly data from 2 sites w/ loadsnotel:
[siteHourly, headers] = loadsnotel('hourly', siteID);

% parse out the date, soil moisture, and soil temperature sensors
datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
sm2= siteHourly{4};
sm8 = siteHourly{5};
sm20 = siteHourly{6};
st2 = siteHourly{7};
st8 = siteHourly{8};
st20 = siteHourly{9};

% load daily data
% [TRIALdaily, tr_headers_d] = loadsnotel('daily', 828);
% [LOUISdaily, lm_headers_d] = loadsnotel('daily', 972);
% tr_datenum_d = datenum(TRIALdaily{2});
% tr_wteq_d = TRIALdaily{4};
% tr_snwd_d = TRIALdaily{11};

% Create logical tests and pull desired multi-month timeperiods
testOND = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
testJFM = (datevec_h(:,2)==1 | datevec_h(:,2)==2 | datevec_h(:,2)==3);
testAM = (datevec_h(:,2)==4 | datevec_h(:,2)==5);
st2_OND = st2(testOND);
st2_JFM = st2(testJFM);
st2_AM = st2(testAM);

% Calculate mean and standard deviations of this data
st2_ONDm = mean(st2_OND(~isnan(st2_OND)));
st2_JFMm = mean(st2_JFM(~isnan(st2_JFM)));
st2_AMm = mean(st2_AM(~isnan(st2_AM)));
st2_ONDsd = std(st2_OND(~isnan(st2_OND)));
st2_JFMsd = std(st2_JFM(~isnan(st2_JFM)));
st2_AMsd = std(st2_AM(~isnan(st2_AM)));

% Logical tests, means and std deviations for months
months = {'JAN' 'FEB' 'MAR' 'APR' 'MAY' 'JUN' 'JUL' 'AUG' 'SEP' 'OCT' 'NOV' 'DEC'};
for i=1:length(months)
    eval(['test' months{i} ' = datevec_h(:,2)==' num2str(i) ';']);
    eval(['st2_' months{i} ' = st2(test' months{i} ');']);
    eval(['st2_' months{i} 'm = mean(st2_' months{i} '(~isnan(st2_' months{i} ')));']);
    eval(['st2_' months{i} 'sd = std(st2_' months{i} '(~isnan(st2_' months{i} ')));']);
end

clear testOND testJFM testAM;

%-------------------------------------------------------------
% PLOT the data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - 5cm temp histograms - All years']);

% Three plots on the left are for OND, JFM, and AM timeperiods
subplot (6, 3, [1 4])
xedges = -2:0.25:20;
n = histc(st2_OND, xedges);
% might want to normalize these
nnorm = n./sum(n);
bar (xedges, nnorm, 'k');
hold on
plot([st2_ONDm st2_ONDm], [0 1], ':k');
axis([-4 20 0 0.6]);
ylabel('Frequency');
title('Oct 1-Dec 31');

subplot (6, 3, [7 10])
xedges = -2:0.25:20;
n = histc(st2_JFM, xedges);
%normalize
nnorm = n./sum(n);
bar (xedges, nnorm, 'b');
hold on
plot([st2_JFMm st2_JFMm], [0 1], ':k');
axis([-4 20 0 0.6]);
ylabel('Frequency');
title('Jan 1 -Mar 31');

subplot (6, 3, [13 16])
xedges = -2:0.25:20;
n = histc(st2_AM, xedges);
% normalize
nnorm = n./sum(n);
bar (xedges, nnorm, 'r');
hold on
plot([st2_AMm st2_AMm], [0 1], ':k');
axis([-4 20 0 0.6]);
ylabel('Frequency');
title('April 1 - May 31');


% Plot the month distributions in 2 columns
plotorder = [2 5 8 11 14 17 3 6 9 12 15 18];

for i = 1:length(plotorder)
    subplot(6, 3, plotorder(i));
    xedges = -2:0.25:20;
    eval(['n = histc(st2_' months{i} ', xedges);']);
    % normalize
    nnorm = n./sum(n);
    bar (xedges, nnorm, 'g');
    axis([-4 20 0 0.5]);
    title(months{i});
end


junk = 99;
