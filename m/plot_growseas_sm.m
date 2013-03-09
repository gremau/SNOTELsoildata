function n = plot_growseas_sm()
% plot_growseas_sm.m
%
% Reads SNOTEL data files and makes plots characterizing growing season 
% variability in soil VWC (entire network or individual sites), and its 
% relationship to snowpack.
%
% Feb 20, 2012 - Greg Maurer

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Add any needed tools
addpath('/home/greg/data/code_resources/m_common/'); 
addpath('/home/greg/data/code_resources/m_common/nanstats/');
%addpath('/home/greg/data/code_resources/m_common/linreg/'); 
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

% Set data path and file name, read in file
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);
% Soil water content data
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt'], 1,0);
%     vwcData = csvread([processeddatapath ...
%         'wyear_soilwatersummary_daily_smnorm.txt'], 1,0);
vwcData = csvread([processeddatapath...
    'wyear_soilwatersummary_hourly.txt'], 1,0);
%     vwcData = csvread([processeddatapath...
%         'wyear_soilwatersummary_daily.txt'], 1,0);

% climData includes more than just soil sites, 
% Get a subset corresponding to the sites and years in tsData
matchsoil = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');

% Aggregation index for climData(matchsoil) data
[soilsites, ~, valindex] = unique(climData(matchsoil,1));
aggindex = [valindex ones(size(climData(matchsoil,1)))];

% Assign climData variables using the headers file USE MATCHSOIL
fid = fopen([processeddatapath 'headersClim.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:length(headers)
    eval([headers{i} ' = climData(matchsoil,i);']);
end
% Convert some units to mm
maxswe = maxswe.*25.4;
accumprecip = accumprecip*25.4;
JASprecip = JASprecip*25.4;

% Assign vwcData variables using the headers file
fid = fopen([processeddatapath 'headersVWC.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:length(headers)
    eval([headers{i} ' = vwcData(:,i);']);
end

% Assign NORMALIZED vwc (vwcDataN) variables using the headers file
fid = fopen([processeddatapath 'headersVWC.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:length(headers)
    eval([headers{i} 'N = vwcDataN(:,i);']);
end

%----------------------------------------------------
% Aggregate some values and select sets of sites
JASprecipAgg = accumarray(aggindex, JASprecip, [numel(soilsites) 1], @nanmean);
maxsweAgg = accumarray(aggindex, maxswe, [numel(soilsites) 1], @nanmean);
elevAgg = accumarray(aggindex, elev, [numel(soilsites) 1], @nanmean);

% Lists of sites with particular elev/swe characteristics - monsites get 
% high summer precip (> 200mm)
sites_hihi = soilsites(elevAgg>2750 & maxsweAgg>500);
sites_lohi = soilsites(elevAgg<2250 & maxsweAgg>500);
sites_hilo = soilsites(elevAgg>2750 & maxsweAgg<300);
sites_lolo = soilsites(elevAgg<2250 & maxsweAgg<300);
mtest = JASprecipAgg > 125;
monsites = soilsites(mtest);
monsites_hihi = monsites(elevAgg(mtest)>3000 & maxsweAgg(mtest)>500);
monsites_lohi = monsites(elevAgg(mtest)<2250 & maxsweAgg(mtest)>500);
monsites_hilo = monsites(elevAgg(mtest)>3000 & maxsweAgg(mtest)<300);
monsites_lolo = monsites(elevAgg(mtest)<2250 & maxsweAgg(mtest)<300);

%----------------------------------------------------
% FIG 1 - Plot soil moisture vs max swe snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe, amjVWC5mean,  '.g', ...
    maxswe, amjVWC20mean, '.b', ...
    maxswe, amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe, jasVWC5mean,  '.g', ...
    maxswe, jasVWC20mean, '.b', ...
    maxswe, jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy, amjVWC5mean,  '.g', ...
    meltdoy, amjVWC20mean, '.b', ...
    meltdoy, amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy, jasVWC5mean,  '.g', ...
    meltdoy, jasVWC20mean, '.b', ...
    meltdoy, jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');

%----------------------------------------------------
% FIG 2 - Plot NORMALIZED soil moisture vs max swe  and snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear NORMALIZED seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe, amjVWC5meanN,  '.g', ...
    maxswe, amjVWC20meanN, '.b', ...
    maxswe, amjVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe, jasVWC5meanN,  '.g', ...
    maxswe, jasVWC20meanN, '.b', ...
    maxswe, jasVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy, amjVWC5meanN,  '.g', ...
    meltdoy, amjVWC20meanN, '.b', ...
    meltdoy, amjVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy, jasVWC5meanN,  '.g', ...
    meltdoy, jasVWC20meanN, '.b', ...
    meltdoy, jasVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');


%----------------------------------------------------
% FIG 3 - Plot soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Jul, Aug, Sep VWC vs snowpack');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 145; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

x = meltdoy; %split into x and y
y = jasVWC5mean;
y2 = jasVWC5sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC');

y = jasVWC20mean;
y2 = jasVWC20sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
ylabel('20cm VWC');

% 50cm
y = jasVWC50mean;
y2 = jasVWC50sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(325); % 50cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
xlabel('Snowmelt day'); ylabel('50cm VWC');

% Right side - plot vs max SWE
topEdge = 1600; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+53;

x = maxswe; %split into x and y
y = jasVWC5mean;
y2 = jasVWC5sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
title('Summer VWC vs peak SWE');
ylabel('5cm VWC');

y = jasVWC20mean;
y2 = jasVWC20sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
ylabel('20cm VWC');

y = jasVWC50mean;
y2 = jasVWC50sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(326);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 4 - Plot NORMALIZED soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Normalized Jul, Aug, Sep VWC vs snowpack');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 145; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

x = meltdoy; %split into x and y
y = jasVWC5meanN;
y2 = jasVWC5sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jasVWC20meanN;
y2 = jasVWC20sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

% 50cm
y = jasVWC50meanN;
y2 = jasVWC50sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(325); % 50cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
xlabel('Snowmelt day'); ylabel('50cm VWC-norm');

% Right side - plot vs max SWE
topEdge = 1600; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+53;

x = maxswe; %split into x and y
y = jasVWC5meanN;
y2 = jasVWC5sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 
ylabel('5cm VWC');

y = jasVWC20meanN;
y2 = jasVWC20sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
ylabel('20cm VWC');

y = jasVWC50meanN;
y2 = jasVWC50sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(326);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 5 - Plot NORMALIZED soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Normalized Jul, Aug, Sep VWC vs snowpack');

polyorder = 1;

x = meltdoy; %split into x and y
y = jasVWC5meanN;
subplot(321); % 5cm  
errorbar(x, y, jasVWC5sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jasVWC20meanN;
subplot(323); % 20cm  
errorbar(x, y, jasVWC20sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

y = jasVWC50meanN;
subplot(325); % 50cm  
errorbar(x, y, jasVWC50sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('50cm VWC-norm');

% Right side - plot vs max SWE

x = maxswe; %split into x and y
y = jasVWC5meanN;
subplot(322); % 5cm  
errorbar(x, y, jasVWC5sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 

y = jasVWC20meanN;
subplot(324); % 20cm  
errorbar(x, y, jasVWC20sdN, 'ok', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);

y = jasVWC50meanN;
subplot(326); % 50cm  
errorbar(x, y, jasVWC50sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 6 - Plot normalized August soil moisture vs snowmelt date
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Aug soil moisture (normalized) vs snomelt date');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 175; % define limits
numBins = 12; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

% 5cm
x = meltdoy; %split into x and y
y = augVWC5meanN;
y2 = augVWC5sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]); %ylim([0 1]);
ticks = [175;240;300];
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
title('5cm VWC');

y = augVWC20meanN;
y2 = augVWC20sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]);%ylim([0 1]);
ticks = [180;240;300];
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
title('20cm VWC');

y = augVWC50meanN;
y2 = augVWC50sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]);%ylim([0 1]);
ticks = [180;240;300];
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
ylabel('Aug VWC (normalized)')
title('50cm VWC');

%--------------------------------------------------------------
% FIG 7 - Regress snowpack vs growing season vwc for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress vwc vs snowpack for 3 sites');
test = siteClim==828;
site1meltd = meltdoy(test);
site1pks = maxswe(test);
site1sm20 = jasVWC20meanN(siteVWC==828);
test = siteClim==330;
site2meltd = meltdoy(test);
site2pks = maxswe(test);
site2sm20 = jasVWC20meanN(siteVWC==330);
test = siteClim==582;
site3meltd = meltdoy(test);
site3pks = maxswe(test);
site3sm20 = jasVWC20meanN(siteVWC==582);

subplot (3, 2, 1)
plot(site1meltd, site1sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1meltd, site1sm20, 1, xrange);
plot(xfit, yfit,':k');
text(240, .8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Trial Lake');

subplot (3, 2, 3)
plot(site2meltd, site2sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2meltd, site2sm20, 1, xrange);
plot(xfit, yfit,':k');
text(210, .2, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Beaver Divide');

subplot (3, 2, 5)
plot(site3meltd, site3sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3meltd, site3sm20, 1, xrange);
plot(xfit, yfit,':k');
text(210, .45, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Little Bear');

subplot (3, 2, 2)
plot(site1pks, site1sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1pks, site1sm20, 1, xrange);
plot(xfit, yfit,':k');
text(200, .7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Trial Lake');

subplot (3, 2, 4)
plot(site2pks, site2sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2pks, site2sm20, 1, xrange);
plot(xfit, yfit,':k');
text(350, .2, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Beaver Divide');

subplot (3, 2, 6)
plot(site3pks, site3sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3pks, site3sm20, 1, xrange);
plot(xfit, yfit,':k');
text(250, .45, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Little Bear');

% -----------------------------------------------------------
% FIG 8 - Month soil water content vs elevation
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', ['Mean monthly 20cm VWC (normalized) vs elevation - '...
    'All SNOTEL/years']);

% Set some plotting parameters
plotorder = 1:3;
months = ['oct';'nov';'dec';'jan';'feb';'mar';'apr';'may';...
    'jun';'jul';'aug';'sep'] ;
polyorder = 1;

for i=plotorder
    subplot(4, 3, i);
    x = elev;
    eval(['y1 = ' months(i,:) 'VWC20meanN;']);
    %eval(['y2 = ' months(i,:) 'Ts5mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
    %plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+3)
    eval(['y1 = ' months(i+3,:) 'VWC20meanN;']);
    %eval(['y2 = ' months(i+3,:) 'Ts20mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+3,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+6);
    eval(['y1 = ' months(i+6,:) 'VWC20meanN;']);
    %eval(['y2 = ' months(i+6,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+6,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end

    subplot(4,3, i+9);
    eval(['y1 = ' months(i+9,:) 'VWC20meanN;']);
    %eval(['y2 = ' months(i+9,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+9,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    if i==1
        ylabel('^oC')
    elseif i==2
        xlabel('Elevation (m)');
        set(gca, 'YTickLabel', '');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
end

% -----------------------------------------------------------------------
% Examine the seasonal variability in temp or vwc at a set of 4 sites.
% These sites have been chosen to represent elevation/temp and SWE
% gradients
%
% 2 plots = Sensor timeseries for each site and then seasonal histograms

% 828 = TrialLk, 333 = BenLomTrail, 452=DonkeyRes, 573=Big BendNV
sensoroutput = 'vwc';
sensordepth = 2; %(1=5cm, 2=20cm, 3=50cm);
startwy = 2006;

% Select TEMP or VWC data and set distribution bins and plot axes
if strcmpi(sensoroutput, 'vwc');
    sensorcolumn = sensordepth + 3; % get proper column using sensordepth
    xmin = 0;
    % If running RAW SENSOR DATA (no normalization)
    % xedges = 0:1:100; % raw sm data bins (0-100)
    % xmax = 75
    % If running NORMALIZED data with smnormalize
    xedges = 0:0.02:1; % normalized vwc bins (0-1)
    xmax = 1; % these axes are good for normalized data
    ymax = 0.25;
    disp('*** Running in normalized soil moisture data mode ***');
elseif strcmpi(sensoroutput, 'temp');
    sensorcolumn = sensordepth + 6; % get proper column using sensordepth
    xedges = -15:0.5:35; % soil temp bins -15-35 deg C
    xmax = 35; % corresponding x axis limits
    xmin = -15;
    ymax = 0.2;
end

% Set up PLOT 1 - add each site's timeseries on iteration through following
% loop
figure3 = figure();

% Allocate for histogram and mean matrices - fill in following loop
histograms = zeros(length(xedges), 16);
means = zeros(16, 1);

for i = 1:length(siteIDs);
    % Load hourly data from site  w/ loadsnotel:
    siteHourly = loadsnotel(siteIDs(i), 'hourly', 'exclude');
    % Get rid of wateryears prior to startwy
    wyexclude = siteHourly{10}>startwy-1;
    for j = 1:10
        siteHourly{j} = siteHourly{j}(wyexclude);
    end
    % Parse out the desired sensor depth, normalize if plotting vwc
    if strcmpi(sensoroutput, 'vwc')
        sensordata = filterseries(siteHourly{sensorcolumn}, 'sigma', 25, 3);
        % SPECIAL Case for Taylor Cyn - There is some bad data that makes it
        % past filter and messes up normalization - remove it
        if siteIDs(i) == 336 %811n for TaylorCyn, 336 for BigBend
            test = sensordata<10; %18 for TaylorCyn, 10 for BigBend
            sensordata(test)=nan;
        end
        sensordata = smnormalize(sensordata, 1);
    else
        sensordata = siteHourly{sensorcolumn};
    end
    
    % Create date arrays
    datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}),...
        'yyyy-mm-ddHH:MM');
    datenum_h = datenum(datevec_h);
    
    % PLOT 1 - add the entire timeseries for site i
    ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
    set(figure3, 'Name', ['Site ' num2str(siteIDs(i))...
        ' - Full sensor timeseries']);
    subplot (4, 1, i)
    plot(datenum_h, sensordata, 'k');
    set(gca, 'XTick', ticklocations);
    set(gca, 'XTickLabel', ticklocations);
    datetick('x', 12, 'keepticks');
    ylabel(sensoroutput);
    title(['Site ' num2str(siteIDs(i))]);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
    testJFM = (datevec_h(:,2)==1 | datevec_h(:,2)==2 | datevec_h(:,2)==3);
    testAMJ = (datevec_h(:,2)==4 | datevec_h(:,2)==5 | datevec_h(:,2)==6);
    testJAS = (datevec_h(:,2)==7 | datevec_h(:,2)==8 | datevec_h(:,2)==9);
    sensordata_OND = sensordata(testOND,:);
    sensordata_JFM = sensordata(testJFM,:);
    sensordata_AMJ = sensordata(testAMJ,:);
    sensordata_JAS = sensordata(testJAS,:);
    
    % Generate histograms for each quarter's sensor data
    histOND = histc(sensordata_OND, xedges);
    histJFM = histc(sensordata_JFM, xedges);
    histAMJ = histc(sensordata_AMJ, xedges);
    histJAS = histc(sensordata_JAS, xedges);
    
    % Normalize and put histograms in the histograms matrix
    % (in plotting order)
    histograms(:, i) = histOND./sum(histOND);
    histograms(:, i+4) = histJFM./sum(histJFM);
    histograms(:, i+8) = histAMJ./sum(histAMJ);
    histograms(:, i+12) = histJAS./sum(histJAS);
    
    % Calculate mean and standard deviations of data from each quarter and
    % place in appropriate vector (in plotting order)
    means(i) = mean(sensordata_OND(~isnan(sensordata_OND)));
    means(i+4) = mean(sensordata_JFM(~isnan(sensordata_JFM)));
    means(i+8) = mean(sensordata_AMJ(~isnan(sensordata_AMJ)));
    means(i+12) = mean(sensordata_JAS(~isnan(sensordata_JAS)));
    % stddevOND = std(sensordata_OND(~isnan(sensordata_OND)));
    % stddevJFM = std(sensordata_JFM(~isnan(sensordata_JFM)));
    % stddevAMJ = std(sensordata_AMJ(~isnan(sensordata_AMJ)));
    % stddevJAS = std(sensordata_JAS(~isnan(sensordata_JAS)));

end
    clear testOND testJFM testAMJ testJAS;

% PLOT 2. Plot quarterly distributions for all sites
titlelabels = {'Hi SWE/Hi Elev' 'Hi SWE/Low Elev' 'Low SWE/Hi Elev'...
    'Low SWE/Low Elev'};
monthlabels = {'Oct-Dec' '' '' '' 'Jan-Mar' '' '' '' 'Apr-Jun' '' '' ''...
    'Jul-Sep'};

figure4 = figure('position',[100 0 1000 800],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');
set(figure4, 'Name', ['4 SNOTEL Sites - ' sensoroutput ...
    ' quarterly histograms - all years combined']);
set(figure4, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);

% Loop through 16 subplots and plot histograms and means
for i = 1:16;
    subplot (4, 4, i)
    bar (xedges, histograms(:, i), 'Facecolor', [0.7 0.7 0.7]);
    hold on
    plot([means(i) means(i)], [0 1], '--k', 'Linewidth', 1.5);
    axis([xmin xmax 0 ymax]);
    set(gca, 'position', [0.925 0.925 1.15 1.19] .* get(gca, 'position'),...
        'Xtick', [0.5]);
    if i==1 || i==5 || i==9 || i==13;
        %ylabel('Frequency');
        text(0.1, 0.8, monthlabels(i), 'Units', 'normalized');
        set(gca, 'XTick', [0;0.5])
        if i==5;
            text(-0.35, -1,'Normalized frequency of ocurrence',...
                'Units', 'normalized', 'Rotation', 90);
        end
    elseif i==16;
        set(gca, 'Xtick', [0.5;1]);
        set(gca, 'YtickLabel', '');
    else
        set(gca, 'YtickLabel', '');
    end
    if i < 5
        %title({['Site ' num2str(siteIDs(i))]; titlelabels{i}},...
        %    'Fontsize', 18, 'Fontangle', 'italic');
        title(titlelabels{i}, 'Fontsize', 18, 'Fontangle', 'italic');
        set(gca, 'XtickLabel', '');
    elseif i < 13
        set(gca, 'XtickLabel', '');
    elseif i 
    end
end

figpath = '../figures/';
print(figure4,'-depsc2','-painters',[figpath 'figJ.eps'])



