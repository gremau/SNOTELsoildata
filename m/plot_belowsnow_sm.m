function plot_belowsnow_sm
% plot_belowsnow_sm.m
%
% Reads SNOTEL data files and makes plots characterizing variability 
% in below-snow soil VWC (entire network or individual sites), and its 
% relationship to snowpack.
%
% Feb 14, 2013 - Greg Maurer
%

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Ask user whether use normalized soil moisture data
normalize = input('Use normalized soil moisture data?  (y/n) : ', 's');

% Access to nan situff, lines, etc
addpath('~/data/code_resources/m_common/');
addpath('~/data/code_resources/m_common/nanstuff/');
addpath('~/data/code_resources/m_common/hline_vline/');
addpath('~/data/code_resources/m_common/linreg/');
%addpath('/home/greg/data/code_resources/m_common/hist2/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));
allsites = unique(dailysites(:,1));
soilsites = unique(soilsites(:,1));

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);
% Soil water content data
if strcmpi(normalize, 'y')
    vwcData = csvread([processeddatapath ...
        'wyear_soilwatersummary_hourly_smnorm.txt'], 1,0);
%     vwcData = csvread([processeddatapath ...
%         'wyear_soilwatersummary_daily_smnorm.txt'], 1,0);
else
    vwcData = csvread([processeddatapath...
        'wyear_soilwatersummary_hourly.txt'], 1,0);
%     vwcData = csvread([processeddatapath...
%         'wyear_soilwatersummary_daily.txt'], 1,0);
end

% climData includes more than just soil sites, 
% Get a subset corresponding to the sites and years in tsData
matchsoil = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);

% Aggregation index for climData 
[vals, ~, valindex] = unique(climData(:,1));
climAggindex = [valindex ones(size(climData(:,1)))];

% Aggregation index for soilClim data
[vals, ~, valindex] = unique(soilClim(:,1));
soilAggindex = [valindex ones(size(soilClim(:,1)))];

% Unique elevations in climData
elevAllAgg = accumarray(climAggindex, climData(:,89), [numel(allsites) 1], @mean);
%Unique elevations in soilClim
elevSoilAgg = accumarray(soilAggindex, soilClim(:,89), [numel(soilsites) 1], @mean);

site_cl = climData(:, 1);
year_cl = climData(:, 2);
maxswe = climData(:, 3)*25.4;
maxsweday = climData(:, 4);
maxdepth = climData(:, 5);
onsetdoy = climData(:, 6);
meltdoy = climData(:, 7);
snowduration = climData(:, 8);
totaldaysSC = climData(:, 9); % Total days, may vary from duration above
maxcontinSC = climData(:, 10);% Length of longest continuos snowpack
numcontinSC = climData(:, 11);% # of continuous snowcovered periods
accumprecip = climData(:, 12)*25.4;
JASprecip = climData(:, 13)*25.4;

octSWEmean = climData(:, 19);
octSWEmed = climData(:, 20);
octSWEsd = climData(:, 21);
novSWEmean = climData(:, 22);
novSWEmed = climData(:, 23);
novSWEsd = climData(:, 24);
decSWEmean = climData(:, 25);
decSWEmed = climData(:, 26);
decSWEsd = climData(:, 27);
janSWEmean = climData(:, 28);
janSWEmed = climData(:, 29);
janSWEsd = climData(:, 30);
febSWEmean = climData(:, 31);
febSWEmed = climData(:, 32);
febSWEsd = climData(:, 33);
marSWEmean = climData(:, 34);
marSWEmed = climData(:, 35);
marSWEsd = climData(:, 36);
aprSWEmean = climData(:, 37);
aprSWEmed = climData(:, 38);
aprSWEsd = climData(:, 39);
maySWEmean = climData(:, 40);
maySWEmed = climData(:, 41);
maySWEsd = climData(:, 42);
junSWEmean = climData(:, 43);
junSWEmed = climData(:, 44);
junSWEsd = climData(:, 45);
julSWEmean = climData(:, 46);
julSWEmed = climData(:, 47);
julSWEsd = climData(:, 48);

octTairMean = climData(:, 55);
octTairSd = climData(:, 56);
novTairMean = climData(:, 57);
novTairSd = climData(:, 58);
decTairMean = climData(:, 59);
decTairSd = climData(:, 60);
janTairMean = climData(:, 61);
janTairSd = climData(:, 62);
febTairMean = climData(:, 63);
febTairSd = climData(:, 64);
marTairMean = climData(:, 65);
marTairSd = climData(:, 66);
aprTairMean = climData(:, 67);
aprTairSd = climData(:, 68);
mayTairMean = climData(:, 69);
mayTairSd = climData(:, 70);
junTairMean = climData(:, 71);
junTairSd = climData(:, 72);
julTairMean = climData(:, 73);
julTairSd = climData(:, 74);
augTairMean = climData(:, 75);
augTairSd = climData(:, 76);
sepTairMean = climData(:, 77);
sepTairSd = climData(:, 78);
freemat = climData(:, 79);
scovmat = climData(:, 80);
maat = climData(:, 81);
maat_sd = climData(:, 82);

preonsetTair = climData(:, 83);
preonsetTairSd = climData(:, 84);
premeltTair = climData(:, 85);
premeltTairSd = climData(:, 86);
postmeltTair = climData(:, 87);
postmeltTairSd = climData(:, 88);
elev = climData(:, 89);
lat = climData(:, 90);
lon = climData(:, 91);
ltMeanSWE = climData(:, 92);
ltMeanPrecip = climData(:, 93);

% Parse soilwatersummary
site_sw = vwcData(:, 1);
year_sw = vwcData(:, 2);
octVWC5mean = vwcData(:, 3);
octVWC5sd = vwcData(:, 4);
octVWC20mean = vwcData(:, 5);
octVWC20sd = vwcData(:, 6);
octVWC50mean = vwcData(:, 7);
octVWC50sd = vwcData(:, 8);
novVWC5mean = vwcData(:, 9);
novVWC5sd = vwcData(:, 10);
novVWC20mean = vwcData(:, 11);
novVWC20sd = vwcData(:, 12);
novVWC50mean = vwcData(:, 13);
novVWC50sd = vwcData(:, 14);
decVWC5mean = vwcData(:, 15);
decVWC5sd = vwcData(:, 16);
decVWC20mean = vwcData(:, 17);
decVWC20sd = vwcData(:, 18);
decVWC50mean = vwcData(:, 19);
decVWC50sd = vwcData(:, 20);
janVWC5mean = vwcData(:, 21);
janVWC5sd = vwcData(:, 22);
janVWC20mean = vwcData(:, 23);
janVWC20sd = vwcData(:, 24);
janVWC50mean = vwcData(:, 25);
janVWC50sd = vwcData(:, 26);
febVWC5mean = vwcData(:, 27);
febVWC5sd = vwcData(:, 28);
febVWC20mean = vwcData(:, 29);
febVWC20sd = vwcData(:, 30);
febVWC50mean = vwcData(:, 31);
febVWC50sd = vwcData(:, 32);
marVWC5mean = vwcData(:, 33);
marVWC5sd = vwcData(:, 34);
marVWC20mean = vwcData(:, 35);
marVWC20sd = vwcData(:, 36);
marVWC50mean = vwcData(:, 37);
marVWC50sd = vwcData(:, 38);
aprVWC5mean = vwcData(:, 39);
aprVWC5sd = vwcData(:, 40);
aprVWC20mean = vwcData(:, 41);
aprVWC20sd = vwcData(:, 42);
aprVWC50mean = vwcData(:, 43);
aprVWC50sd = vwcData(:, 44);
mayVWC5mean = vwcData(:, 45);
mayVWC5sd = vwcData(:, 46);
mayVWC20mean = vwcData(:, 47);
mayVWC20sd = vwcData(:, 48);
mayVWC50mean = vwcData(:, 49);
mayVWC50sd = vwcData(:, 50);
junVWC5mean = vwcData(:, 51);
junVWC5sd = vwcData(:, 52);
junVWC20mean = vwcData(:, 53);
junVWC20sd = vwcData(:, 54);
junVWC50mean = vwcData(:, 55);
junVWC50sd = vwcData(:, 56);
julVWC5mean = vwcData(:, 57);
julVWC5sd = vwcData(:, 58);
julVWC20mean = vwcData(:, 59);
julVWC20sd = vwcData(:, 60);
julVWC50mean = vwcData(:, 61);
julVWC50sd = vwcData(:, 62);
augVWC5mean = vwcData(:, 63);
augVWC5sd = vwcData(:, 64);
augVWC20mean = vwcData(:, 65);
augVWC20sd = vwcData(:, 66);
augVWC50mean = vwcData(:, 67);
augVWC50sd = vwcData(:, 68);

% These repeat through sept (end of wy)
ondVWC5mean = vwcData(:, 75);
ondVWC5sd = vwcData(:, 76);
ondVWC20mean = vwcData(:, 77);
ondVWC20sd = vwcData(:, 78);
ondVWC50mean = vwcData(:, 79);
ondVWC50sd = vwcData(:, 81);
jfmVWC5mean = vwcData(:, 81);
jfmVWC5sd = vwcData(:, 82);
jfmVWC20mean = vwcData(:, 83);
jfmVWC20sd = vwcData(:, 84);
jfmVWC50mean = vwcData(:, 85);
jfmVWC50sd = vwcData(:, 86);
amjVWC5mean = vwcData(:, 87);
amjVWC5sd = vwcData(:, 88);
amjVWC20mean = vwcData(:, 89);
amjVWC20sd = vwcData(:, 90);
amjVWC50mean = vwcData(:, 91);
amjVWC50sd = vwcData(:, 92);
jasVWC5mean = vwcData(:, 93);
jasVWC5sd = vwcData(:, 94);
jasVWC20mean = vwcData(:, 95);
jasVWC20sd = vwcData(:, 96);
jasVWC50mean = vwcData(:, 97);
jasVWC50sd = vwcData(:, 98);

preonsetVWC5 = vwcData(:, 99);
preonsetVWC5sd = vwcData(:, 100);
preonsetVWC20 = vwcData(:, 101);
preonsetVWC20sd = vwcData(:, 102);
preonsetVWC50 = vwcData(:, 103);
preonsetVWC50sd = vwcData(:, 104);

% Parse soiltempsummary
site_st = tsData(:, 1);
year_st = tsData(:, 2);
octTs5mean = tsData(:, 3);
octTs5sd = tsData(:, 4);
octTs20mean = tsData(:, 5);
octTs20sd = tsData(:, 6);
octTs50mean = tsData(:, 7);
octTs50sd = tsData(:, 8);
decTs5mean = tsData(:, 15);
decTs5sd = tsData(:, 16);
decTs20mean = tsData(:, 17);
decTs20sd = tsData(:, 18);
decTs50mean = tsData(:, 19);
decTs50sd = tsData(:, 20);
% These repeat through sept (end of wy)
ondTs5mean = tsData(:, 75);
ondTs5sd = tsData(:, 76);
ondTs20mean = tsData(:, 77);
ondTs20sd = tsData(:, 78);
ondTs50mean = tsData(:, 79);
ondTs50sd = tsData(:, 80);
jfmTs5mean = tsData(:, 81);
jfmTs5sd = tsData(:, 82);
jfmTs20mean = tsData(:, 83);
jfmTs20sd = tsData(:, 84);
jfmTs50mean = tsData(:, 85);
jfmTs50sd = tsData(:, 86);
amjTs5mean = tsData(:, 87);
amjTs5sd = tsData(:, 88);
amjTs20mean = tsData(:, 89);
amjTs20sd = tsData(:, 90);
amjTs50mean = tsData(:, 91);
amjTs50sd = tsData(:, 92);
jasTs5mean = tsData(:, 93);
jasTs5sd = tsData(:, 94);
jasTs20mean = tsData(:, 95);
jasTs20sd = tsData(:, 96);
jasTs50mean = tsData(:, 97);
jasTs50sd = tsData(:, 98);

% Seasonal/yearly soil temp means
mast5cm = tsData(:, 99);
sdast5cm = tsData(:, 100);
mast20cm = tsData(:, 101);
sdast20cm = tsData(:, 102);
mast50cm = tsData(:, 103);
sdast50cm = tsData(:, 104);

% Snowcovered soil temp means
snowcovTs5mean = tsData(:, 105);
snowcovTs5sd = tsData(:, 106);
snowcovTs20mean = tsData(:, 107);
snowcovTs20sd = tsData(:, 108);
snowcovTs50mean = tsData(:, 109);
snowcovTs50sd = tsData(:, 110);

% Snowfree soil temp means
snowfreeTs5mean = tsData(:, 111);
snowfreeTs5sd = tsData(:, 112);
snowfreeTs20mean = tsData(:, 113);
snowfreeTs20sd = tsData(:, 114);
snowfreeTs50mean = tsData(:, 115);
snowfreeTs50sd = tsData(:, 116);

preonsetTs5 = tsData(:, 117);
preonsetTs5sd = tsData(:, 118);
preonsetTs20 = tsData(:, 119);
preonsetTs20sd = tsData(:, 120);
preonsetTs50 = tsData(:, 121);
preonsetTs50sd = tsData(:, 122);
premeltTs5 = tsData(:, 123);
premeltTs5sd = tsData(:, 124);
postmeltTs5 = tsData(:, 125);
postmeltTs5sd = tsData(:, 126);

%-------------------------------------------------------------------------

% FIG 10 - Plot winter soil moisture vs onset soil moisture
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter VWC vs pre-snowpack VWC');

polyorder = 1;
plotorder = 1:4;
months = ['Oct';'Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May'];

x1 = preonsetVWC5;
x2 = preonsetVWC20;
x3 = preonsetVWC50;

for i=plotorder
    subplot(2, 4, i)
    eval(['y1 = ' lower(months(i,:)) 'VWC5mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC20mean;']);
    eval(['y3 = ' lower(months(i,:)) 'VWC50mean;']);
    
    plot(x1, y1, '.', 'Color', [0.7 0.7 0.7]);
    hold on;
    plot(x2, y2, '.', 'Color', [0.5 0.5 0.5]);
    plot(x3, y3, '.', 'Color', [0.3 0.3 0.3]);
    [b, rsq, xfit, yfit] = fitline(x2, y2, polyorder, [0, 1]);
    [b2,bint,resid,rint,stats] = regress2(y2, [x2 ones(size(x2))]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    % slope and r^2 values
    text(0.5, 0.08, ['Slope = ' num2str(b(1), 2) ...
        ', r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3))]);
    %xlim([0, 1500]); ylim([-10, 10]);
    title(months(i,:));
    set(gca, 'XTickLabel', '');
    if i==1
        legend('Location', 'Northwest', '5cm', '20cm', '50cm', '20cm polyfit');
        text(-0.35, -0.45, 'Mean monthly soil VWC', 'Units', 'normalized', 'Rotation', 90);
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(2, 4, i+4)
    eval(['y1 = ' lower(months(i+4,:)) 'VWC5mean;']);
    eval(['y2 = ' lower(months(i+4,:)) 'VWC20mean;']);
    eval(['y3 = ' lower(months(i+4,:)) 'VWC50mean;']);
    
    plot(x1, y1, '.', 'Color', [0.7 0.7 0.7]);
    hold on;
    plot(x2, y2, '.', 'Color', [0.5 0.5 0.5]);
    plot(x3, y3, '.', 'Color', [0.3 0.3 0.3]);
    [b, rsq, xfit, yfit] = fitline(x2, y2, polyorder, [0, 1]);
    [b2,bint,resid,rint,stats] = regress2(y2, [x2 ones(size(x2))]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    % slope and r^2 values
    text(0.5, 0.08, ['Slope = ' num2str(b(1), 2) ...
        ', r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3))]);
    %xlim([0, 1500]); ylim([-10, 10]);
    title(months(i+4,:));
    if i>1
        set(gca, 'YTickLabel', '');
    end
    if i==2
        text(1, -0.2, 'Pre-onset VWC', 'Units', 'normalized');
    end
end

%----------------------------------------------------
% FIG 11 - Plot winter soil moisture vs onset soil moisture
% Same as above, but binned, using only 20cm data
fignum = fignum+1;
fig= figure('position',[100 0 400 800],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Winter VWC vs pre-snowpack VWC');
set(fig, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);

% Set binning parameters
if strcmpi(normalize, 'n')
    topEdge = 45; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 45];
    yaxlim = [0 45];
elseif strcmpi(normalize, 'y')
    topEdge = 1; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 1];
    yaxlim = [0 1];
end
    
% And an xaxis to use
xax = (linspace(topEdge/numBins, topEdge, numBins) - (topEdge/numBins)/2);

% And months to plot
months = ['Nov'; 'Jan'; 'Mar'];

for i = 1:3;
    subplot (3,1, i)
    x = preonsetVWC5; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC5mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC5sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o', 'MarkerSize', 8,...
        'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6]);
    hold on;
    x = preonsetVWC20; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC20mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC20sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o', 'MarkerSize', 8,...
        'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3]);
    x = preonsetVWC50; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC50mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC50sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'ok',...
        'MarkerFaceColor', 'k', 'MarkerSize', 8);
    % Plot 1:1 line
    plot(0:45, 0:45, ':k');
    pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',pos.*[1 .90 1 1.22]); % change its height
    ylabel([months(i,:) ' VWC (normalized)']);
    xlim(xaxlim); ylim(yaxlim);
    if i==1
        legend('5cm', '20cm', '50cm', 'Location', 'NorthWest');
        title('Month mean VWC vs pre-snowpack VWC');
        set(gca, 'XTickLabels', '');
    elseif i < 3
        set(gca, 'XTickLabels', '');
    elseif i==3
        xlabel('Pre-snowpack VWC (normalized)');
    end
    
end

figpath = '../figures/';
print(fig,'-depsc2','-painters',[figpath 'figI.eps'])

%----------------------------------------------------
% FIG 12 - Plot winter soil moisture vs onset soil moisture
% Same as above, but binned, using only 20cm data
fignum = fignum+1;
fig= figure('position',[100 0 600 500],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Winter VWC vs pre-snowpack VWC');
set(fig, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);

% Set binning parameters
if strcmpi(normalize, 'n')
    topEdge = 45; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 45];
    yaxlim = [0 45];
elseif strcmpi(normalize, 'y')
    topEdge = 1; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 1];
    yaxlim = [0 1];
end
    
% And an xaxis to use
xax = (linspace(topEdge/numBins, topEdge, numBins) - (topEdge/numBins)/2);

% And months and colors to plot
months = ['Oct'; 'Dec'; 'Feb'; 'Apr'];
colors = {[0.9 0.9 0.9]; [0.6 0.6 0.6]; [0.3 0.3 0.3]; 'k'};

for i = 1:4;
    x = preonsetVWC20; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC20mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC20sd;']);
    [binMean, binSd] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMean, binSd, 'ok', 'MarkerSize', 10,...
        'MarkerFaceColor', colors{i});
    hold on;
end
% Plot 1:1 line
plot(0:45, 0:45, ':k', 'linewidth', 2);
pos = get(gca,'position'); % get subplot axis position
%set(gca,'position',pos.*[1 .90 1 1.22]); % change its height
ylabel('Mean monthly VWC (%)');
xlim(xaxlim); ylim(yaxlim);
legend('October','December','February','April','Location','Southeast');
xlabel('Pre-onset VWC (normalized)');


figpath = '../figures/';
print(fig,'-depsc2','-painters',[figpath 'figI2.eps'])

%----------------------------------------------------
% FIG 13 - Plot winter soil moisture vs onset soil moisture
% Same as above, with all months
fignum = fignum+1;
fig= figure('position',[100 0 600 500],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Winter VWC vs pre-snowpack VWC');
set(fig, 'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% Set binning parameters
if strcmpi(normalize, 'n')
    topEdge = 45; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 45];
    yaxlim = [0 45];
elseif strcmpi(normalize, 'y')
    topEdge = 1; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 1];
    yaxlim = [0 1];
end
    
% And an xaxis to use
xax = (linspace(topEdge/numBins, topEdge, numBins) - (topEdge/numBins)/2);

% And months and colors to plot
months = ['Oct';'Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May'];
colors = {'LightCyan';'LightSkyBlue';'DeepSkyBlue';'DodgerBlue';'MediumBlue';'DarkBlue';...
    'MidnightBlue';'Black';};
grays = {[0.9 0.9 0.9];[0.7 0.7 0.7]; [0.6 0.6 0.6]; [0.5 0.5 0.5]; ...
    [0.4 0.4 0.4]; [0.3 0.3 0.3];[0.2 0.2 0.2];'k'};

for i = 1:8;
    x = preonsetVWC20; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC20mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC20sd;']);
    [binMean, binSd] = binseries(x, y, y2, topEdge, botEdge, numBins);
    %errorbar(xax(1:numBins), binMean, binSd, 'ok', 'MarkerSize', 10,...
     %   'MarkerFaceColor', rgb(colors{i}));']);
    %eval(['h' int2str(i)  '= errorbar(xax(1:numBins), binMean, binSd,' ...
    %    char(39) 'ok' char(39) ', ' char(39) 'MarkerSize' char(39) ''...
    %    ', 10,' char(39) 'MarkerFaceColor' char(39) ', rgb(colors{i}));']);
    eval(['h' int2str(i)  '= errorbar(xax(1:numBins), binMean, binSd,' ...
        char(39) 'ok' char(39) ', ' char(39) 'MarkerSize' char(39) ''...
        ', 10,' char(39) 'MarkerFaceColor' char(39) ', grays{i});']);
    hold on;
end
% Plot 1:1 line
plot(0:45, 0:45, ':k', 'linewidth', 1.5);
pos = get(gca,'position'); % get subplot axis position
%set(gca,'position',pos.*[1 .90 1 1.22]); % change its height
ylabel('Mean monthly VWC (normalized)');
xlim(xaxlim); ylim(yaxlim);
legend([h1; h8], {'October','May'}, 'Location', 'Southeast');
xlabel('Pre-onset VWC (normalized)');


figpath = '../figures/';
print(fig,'-depsc2','-painters',[figpath 'figI3.eps'])
junk = 99;
end