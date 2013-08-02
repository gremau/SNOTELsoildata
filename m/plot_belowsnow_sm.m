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
addpath('~/data/code_resources/m_common/nanstats/');
addpath('~/data/code_resources/m_common/hline_vline/');
addpath('~/data/code_resources/m_common/linreg/');
%addpath('/home/greg/data/code_resources/m_common/hist2/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
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
matchsoil = ismember(climData(:, 1:2), vwcData(:, 1:2), 'rows');

% Assign climData variables using the headers file
fid = fopen([processeddatapath 'headersClim.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};

for i=1:11
    eval([headers{i} ' = climData(:,i);']);
end
% Load precip + SWE and convert to mm
for i=12:54
    eval([headers{i} ' = climData(:,i)*25.4;']);
end
% and the rest with no conversion
for i=55:length(headers)
    eval([headers{i} ' = climData(:,i);']);
end


% Assign vwcData variables using the headers file
fid = fopen([processeddatapath 'headersVWC.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};

for i=1:length(headers)
    eval([headers{i} ' = vwcData(:,i);']);
end


%-------------------------------------------------------------------------
% FIG 1 - Plot winter soil moisture vs onset soil moisture
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter monthly VWC vs pre-onset VWC');

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

%-------------------------------------------------------------------------
% FIG 2 - Plot winter soil moisture vs prior month's soil moisture
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter monthly VWC vs prior month VWC');

polyorder = 1;
plotorder = 1:4;
months = ['Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May';'Jun'];
priormonths = ['Oct';'Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May'];

for i=plotorder
    subplot(2, 4, i)
    eval(['x1 = ' lower(priormonths(i,:)) 'VWC5mean;']);
    eval(['x2 = ' lower(priormonths(i,:)) 'VWC20mean;']);
    eval(['x3 = ' lower(priormonths(i,:)) 'VWC50mean;']);
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
    xlim([0, 1]); ylim([0, 1]);
    title(months(i,:));
    set(gca, 'XTickLabel', '');
    if i==1
        legend('Location', 'Northwest', '5cm', '20cm', '50cm', '20cm polyfit');
        text(-0.35, -0.45, 'Mean monthly soil VWC', 'Units', 'normalized', 'Rotation', 90);
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(2, 4, i+4)
    eval(['x1 = ' lower(priormonths(i+4,:)) 'VWC5mean;']);
    eval(['x2 = ' lower(priormonths(i+4,:)) 'VWC20mean;']);
    eval(['x3 = ' lower(priormonths(i+4,:)) 'VWC50mean;']);
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
    xlim([0, 1]); ylim([0, 1]);
    title(months(i+4,:));
    if i>1
        set(gca, 'YTickLabel', '');
    end
    if i==2
        text(1, -0.2, 'Prior month VWC', 'Units', 'normalized');
    end
end

%-------------------------------------------------------------------------
% FIG 3 - Plot winter soil moisture vs onset soil moisture
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

% figpath = '../figures/';
% print(fig,'-depsc2','-painters',[figpath 'figI.eps'])

%-------------------------------------------------------------------------
% FIG 4 - Plot winter soil moisture vs onset soil moisture
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


% figpath = '../figures/';
% print(fig,'-depsc2','-painters',[figpath 'figI2.eps'])

%-------------------------------------------------------------------------
% FIG 5 - Plot winter soil moisture vs onset soil moisture
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


% figpath = '../figures/';
% print(fig,'-depsc2','-painters',[figpath 'figI3.eps'])

%-------------------------------------------------------------------------
% FIG 6 - Plot change in monthly soil moisture from pre-onset VWC
fignum = fignum+1;
fig= figure('position',[100 0 600 500],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Monthly VWC - Pre-onset VWC (pre-onset delta VWC)');
set(fig, 'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% Change months a bit
months = ['Nov';'Dec';'Jan';'Feb';'Mar';'Apr';'May'];

for i = 1:length(months);
    eval(['delta5cm' months(i,:) ' = ' lower(months(i,:)) 'VWC5mean' ...
         '- preonsetVWC5;']);
    eval(['delta20cm' months(i,:) ' = ' lower(months(i,:)) 'VWC20mean' ...
         '- preonsetVWC20;']);
     eval(['delta50cm' months(i,:) ' = ' lower(months(i,:)) 'VWC50mean' ...
         '- preonsetVWC50;']);
    n = sum(~isnan(['delta5cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta5cm' months(i,:)...
        '),  nanstd(delta5cm' months(i,:) ')/sqrt(n),' char(39) 'or' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
    hold on;
    n = sum(~isnan(['delta20cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta20cm' months(i,:)...
        '),  nanstd(delta20cm' months(i,:) ')/sqrt(n),' char(39) 'og' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
    n = sum(~isnan(['delta50cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta50cm' months(i,:)...
        '),  nanstd(delta50cm' months(i,:) ')/sqrt(n),' char(39) 'ob' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
end
ylabel('Change in VWC');
xlabel('Month');

%-------------------------------------------------------------------------
% FIG 7 - Plot change in monthly soil moisture from prior month VWC
% Same as above, with all months
fignum = fignum+1;
fig= figure('position',[100 0 600 500],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Monthly VWC - Prior month VWC (1-month delta VWC)');
set(fig, 'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% And prior months
priormonths = ['Oct';'Nov';'Dec';'Jan';'Feb';'Mar';'Apr'];

for i = 1:length(months);
    eval(['delta5cm' months(i,:) ' = ' lower(months(i,:)) 'VWC5mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC5mean;']);
    eval(['delta20cm' months(i,:) ' = ' lower(months(i,:)) 'VWC20mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC20mean;']);
    eval(['delta50cm' months(i,:) ' = ' lower(months(i,:)) 'VWC50mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC50mean;']);
    n = sum(~isnan(['delta5cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta5cm' months(i,:)...
        '),  nanstd(delta5cm' months(i,:) ')/sqrt(n),' char(39) 'or' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
    hold on;
    n = sum(~isnan(['delta20cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta20cm' months(i,:)...
        '),  nanstd(delta20cm' months(i,:) ')/sqrt(n),' char(39) 'og' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
    n = sum(~isnan(['delta50cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar(i, nanmean(delta50cm' months(i,:)...
        '),  nanstd(delta50cm' months(i,:) ')/sqrt(n),' char(39) 'ob' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 10);']);
end
ylabel('Change in VWC');
xlabel('Month');

% -------------------------------------------------------------
% FIG 8 - Dual panel with Fig 6 and 7
fignum = fignum+1;    
h = figure('position',[100 0 1100 500],'paperpositionmode',...
    'auto','color', 'white','InvertHardcopy','off');
set(h, 'Name','Change in below-snow VWC',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

subplot(1,2,1);
% Plot change in VWC from prior month - WATCH OUT for wierd points (error
% bar fix)
for i = 1:length(months);
    eval(['delta5cm' months(i,:) ' = ' lower(months(i,:)) 'VWC5mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC5mean;']);
    eval(['delta20cm' months(i,:) ' = ' lower(months(i,:)) 'VWC20mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC20mean;']);
    eval(['delta50cm' months(i,:) ' = ' lower(months(i,:)) 'VWC50mean ' ...
         '- ' lower(priormonths(i,:)) 'VWC50mean;']);
    n = sum(~isnan(['delta50cm' months(i,:)]));
    eval(['h' int2str(i+2)  ' = errorbar([i;10;1], [nanmean(delta50cm' months(i,:)...
        ');999;999],  [nanstd(delta50cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{7});']);
    hold on;
    n = sum(~isnan(['delta20cm' months(i,:)]));
    eval(['h' int2str(i+1)  ' = errorbar([i;10;1], [nanmean(delta20cm' months(i,:)...
        ');999;999],  [nanstd(delta20cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{4});']);
    n = sum(~isnan(['delta5cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar([i;10;1], [nanmean(delta5cm' months(i,:)...
        ');999;999],  [nanstd(delta5cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{1});']);

end
set(gca, 'position', [0.90 1 1.15 1] .* get(gca, 'position'));
ylim([-0.2 0.6]);
hline(0, ':k');
ylabel('Change in WC (normalized)');
set(gca, 'xtick', [1;2;3;4;5;6;7], 'xticklabel', months);
legend([h7;h8;h9], {'5cm', '20cm', '50cm'}, 'location', 'northeast');
text(0.05, 0.9, 'a. Month-to-month', 'units', 'normalized', 'fontangle', 'italic',...
    'Fontsize', 20);
%text(0.95, -0.11, 'Month', 'units', 'normalized');


subplot(1,2,2);
% Plot Change in VWC from pre-onset - WATCH OUT for wierd points (error
% bar fix)
for i = 1:length(months);
    eval(['delta5cm' months(i,:) ' = ' lower(months(i,:)) 'VWC5mean' ...
         '- preonsetVWC5;']);
    eval(['delta20cm' months(i,:) ' = ' lower(months(i,:)) 'VWC20mean' ...
         '- preonsetVWC20;']);
     eval(['delta50cm' months(i,:) ' = ' lower(months(i,:)) 'VWC50mean' ...
         '- preonsetVWC50;']);
    n = sum(~isnan(['delta50cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar([i;10;1], [nanmean(delta50cm' months(i,:)...
        ');999;999],  [nanstd(delta50cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{7});']);
    hold on;
    n = sum(~isnan(['delta20cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar([i;10;1], [nanmean(delta20cm' months(i,:)...
        ');999;999],  [nanstd(delta20cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{4});']);
    n = sum(~isnan(['delta5cm' months(i,:)]));
    eval(['h' int2str(i)  ' = errorbar([i;10;1], [nanmean(delta5cm' months(i,:)...
        ');999;999],  [nanstd(delta5cm' months(i,:) ')/sqrt(n);100;100],' char(39) 'ok' char(39)...
        ', ' char(39) 'MarkerSize' char(39) ', 12,' char(39) ...
        'MarkerFaceColor' char(39) ', grays{1});']);

end
set(gca, 'position', [0.90 1 1.15 1] .* get(gca, 'position'));
ylim([-0.2 0.6]);
hline(0, ':k');
set(gca, 'yticklabel', [], 'xtick', [1;2;3;4;5;6;7], 'xticklabel', months);
text(0.05, 0.9, 'b. Cumulative', 'Units', 'normalized', 'Fontangle', 'italic',...
    'Fontsize',20);

figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figF.eps']) 

junk = 99;
end