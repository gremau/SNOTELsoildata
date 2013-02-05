%soilseason_lm.m

close all;
clear all;
fignum = 1;

% Add any needed tools
addpath('/home/greg/data/code_resources/m_common/'); 
addpath('/home/greg/data/code_resources/m_common/nanstuff/');
addpath('/home/greg/data/code_resources/m_common/linear/'); 
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

% Use this site as an example:
% examplesite = 828;
% examplelabel = 'Trial Lake';
% examplesite = 393;
% examplelabel = 'Chalk Creek';
% examplesite = 332;
% examplelabel = 'Ben Lomond Peak';
% examplesite = 333;
% examplelabel = 'Ben Lomond Trail';
% examplesite = 390;
% examplelabel = 'Castle Valley';
%examplesite = 766;
%examplelabel = 'Snowbird';
examplesite = 839;
examplelabel = 'Upper Rio Grande';

% Set processed data path
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'],1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'],1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt'],1,0);
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt'],1,0);

% Get a subset of climData that corresponds with available soildata
[matchsoil, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);
% matchsoil2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

% Climate data
site_cl = soilClim(:, 1);
year_cl = soilClim(:, 2);
maxswe = soilClim(:, 3)*25.4;
maxsweday = soilClim(:, 4);
maxdepth = soilClim(:, 5);
onsetdoy = soilClim(:, 6);
meltdoy = soilClim(:, 7);
snowduration = soilClim(:, 8);
totaldaysSC = soilClim(:, 9); % Total days, may vary from duration above
maxcontinSC = soilClim(:, 10);% Length of longest continuos snowpack
numcontinSC = soilClim(:, 11);% # of continuous snowcovered periods
accumprecip = soilClim(:, 12)*25.4;
JASprecip = soilClim(:, 13)*25.4;
maat = soilClim(:, 80);
maat_sd = soilClim(:, 81);

preonsetTair = soilClim(:, 82);
preonsetTairSd = soilClim(:, 83);
premeltTair = soilClim(:, 84);
premeltTairSd = soilClim(:, 85);
postmeltTair = soilClim(:, 86);
postmeltTairSd = soilClim(:, 87);
elev = soilClim(:, 88);
lat = soilClim(:, 89);
lon = soilClim(:, 90);
ltMeanSWE = soilClim(:, 91);
ltMeanPrecip = soilClim(:, 92);

% Seasonal/yearly soil temp means
site_ts = tsData(:, 1);
mast5cm = tsData(:, 99);
sdast5cm = tsData(:, 100);
mast20cm = tsData(:, 101);
sdast20cm = tsData(:, 102);
mast50cm = tsData(:, 103);
sdast50cm = tsData(:, 104);

% Seasonal soil moisture
amjVWC5mean = vwcDataN(:, 87);
amjVWC5sd = vwcDataN(:, 88);
amjVWC20mean = vwcDataN(:, 89);
amjVWC20sd = vwcDataN(:, 90);
amjVWC50mean = vwcDataN(:, 91);
amjVWC50sd = vwcDataN(:, 92);
jasVWC5mean = vwcDataN(:, 93);
jasVWC5sd = vwcDataN(:, 94);
jasVWC20mean = vwcDataN(:, 95);
jasVWC20sd = vwcDataN(:, 96);
jasVWC50mean = vwcDataN(:, 97);
jasVWC50sd = vwcDataN(:, 98);

preonsetVWC5 = vwcDataN(:, 99);
preonsetVWC5sd = vwcDataN(:, 100);
preonsetVWC20 = vwcDataN(:, 101);
preonsetVWC20sd = vwcDataN(:, 102);
preonsetVWC50 = vwcDataN(:, 103);
preonsetVWC50sd = vwcDataN(:, 104);

% Monthly soil moisture
novVWC5mean = vwcDataN(:, 9);
novVWC5sd = vwcDataN(:, 10);
novVWC20mean = vwcDataN(:, 11);
novVWC20sd = vwcDataN(:, 12);
novVWC50mean = vwcDataN(:, 13);
novVWC50sd = vwcDataN(:, 14);
decVWC5mean = vwcDataN(:, 15);
decVWC5sd = vwcDataN(:, 16);
decVWC20mean = vwcDataN(:, 17);
decVWC20sd = vwcDataN(:, 18);
decVWC50mean = vwcDataN(:, 19);
decVWC50sd = vwcDataN(:, 20);
janVWC5mean = vwcDataN(:, 21);
janVWC5sd = vwcDataN(:, 22);
janVWC20mean = vwcDataN(:, 23);
janVWC20sd = vwcDataN(:, 24);
janVWC50mean = vwcDataN(:, 25);
janVWC50sd = vwcDataN(:, 26);
febVWC5mean = vwcDataN(:, 27);
febVWC5sd = vwcDataN(:, 28);
febVWC20mean = vwcDataN(:, 29);
febVWC20sd = vwcDataN(:, 30);
febVWC50mean = vwcDataN(:, 31);
febVWC50sd = vwcDataN(:, 32);
marVWC5mean = vwcDataN(:, 33);
marVWC5sd = vwcDataN(:, 34);
marVWC20mean = vwcDataN(:, 35);
marVWC20sd = vwcDataN(:, 36);
marVWC50mean = vwcDataN(:, 37);
marVWC50sd = vwcDataN(:, 38);
julVWC5mean = vwcDataN(:, 57);
julVWC5sd = vwcDataN(:, 58);
julVWC20mean = vwcDataN(:, 59);
julVWC20sd = vwcDataN(:, 60);
julVWC50mean = vwcDataN(:, 61);
julVWC50sd = vwcDataN(:, 62);
augVWC5mean = vwcDataN(:, 63);
augVWC5sd = vwcDataN(:, 64);
augVWC20mean = vwcDataN(:, 65);
augVWC20sd = vwcDataN(:, 66);
augVWC50mean = vwcDataN(:, 67);
augVWC50sd = vwcDataN(:, 68);

sites = unique(site_cl);

% PUBLICATION PLOTS   ------------------------------------------------
diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
Yvars = {'mast20cm' 'diff20cm' 'jasVWC20mean','jasVWC20mean','jasVWC20mean'};
Xvars = {'totaldaysSC' 'totaldaysSC' 'meltdoy', 'maxswe', 'JASprecip'};

% EXPLORATORY ANALYSIS  ------------------------------------------------
% Set regression variables:

% Yvars = {'mast50cm','mast50cm','mast50cm'};
% Xvars = {'maat', 'totaldaysSC', 'maxswe'};

% Yvars = {'julVWC20mean','augVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'meltdoy', 'meltdoy'};

% Yvars = {'jasVWC20mean','jasVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'maxswe', 'JASprecip'};

% Yvars = {'novVWC20mean','febVWC20mean','marVWC20mean'};
% Xvars = {'preonsetVWC20', 'preonsetVWC20', 'preonsetVWC20'};

% diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
% Yvars = {'diff5cm','diff20cm','diff50cm'};
% Xvars = {'onsetdoy', 'onsetdoy', 'onsetdoy'};

% diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
% int = totaldaysSC.*maat;
% Yvars = {'mast50cm','mast50cm','mast50cm'};
% Xvars = {'maat', 'totaldaysSC', 'int'};

for i = 1:length(Yvars)
    % Gather regression results for all sites
    slopes = [];
    yints = [];
    rsqs = [];
    pvals = [];
    xmeans = [];
    elevs = [];
    siteIDs = [];
    for j = 1:length(sites)
        getsite = site_cl==sites(j);
        eval(['y = ' Yvars{i} '(getsite);']);
        eval(['x = ' Xvars{i} '(getsite);']);
        x(:,2) = ones(numel(x), 1); % Regress needs 2nd constant x column
        % Only use sites with at least 3 datapoints after nan removal
        nantest = isnan(x(:,1)) | isnan(y);
        if (sum(getsite) - sum(nantest)) > 4 % Do sites w/ this many years
            %[coefficients, rsq, ~, ~] = fitline(x(:,1), y, 1, [0, 1]);
            [b,bint,resid,rint,stats] = shregress(y, x);
            slopes = [slopes; b(1)];%coefficients(1)];
            %slopes(j, i) = b(1);
            yints = [yints; b(2)];%coefficients(2)];
            rsqs = [rsqs; stats(1)];%rsq];
            pvals = [pvals; stats(3)];
            xmeans = [xmeans; nanmean(x(:,1))];
            elevs = [elevs; mean(elev(getsite))];
            siteIDs = [siteIDs; sites(j)];
        end
    end
    
    % Spit out the results
    regressTemp = [siteIDs elevs xmeans slopes yints rsqs pvals];
    eval(['regress' num2str(i) 'dat = regressTemp;']);
    
    % Stouffers combined p value
    pStouff = inline('(1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2');
    
    % Fig 1, 3, 5 - Plot regression results for all sites
    h = figure(fignum);
    set(h, 'Name', ['Regressions: Y = ' Yvars{i} ', X = ' Xvars{i}]);
    
    subplot(3, 3, 1);
    plot(xmeans, slopes, 'ok');
    title('Slope');
    hline(mean(slopes), ':k')
    xlabel('Site mean of x');
    
    subplot(3, 3, 2);
    plot(xmeans, pvals, 'ob');
    hline(mean(pvals), ':b')
    title('P-values');
    text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
        'Units', 'normalized');
    xlabel('Site mean of x');
    
    subplot(3, 3, 3);
    plot(xmeans, rsqs, 'or');
    hline(mean(rsqs), ':r')
    title('R squared');
    xlabel('Site mean of x');

    subplot (3, 3, 4)
    plot(elevs, slopes, '.k');
    xlabel('Elev (m)');
    
    subplot(3, 3, 5);
    plot(elevs, pvals, '.b');
    xlabel('Elev (m)');
    
    subplot(3, 3, 6);
    plot(elevs, rsqs, '.r');
    xlabel('Elev (m)');
    
    subplot(3, 3, 7);
    xedges = linspace(min(slopes), max(slopes), 20);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    %xlim([-0.1 0.1]); 
    %ylim([0 35]);
    vline(mean(slopes), '-r'); vline(0, ':k');
    title('Frequency of Slope values');
    xlabel('Slope'); ylabel('n');
        
    subplot(3, 3, 8);
    xedges = linspace(min(pvals), max(pvals), 20);
    pHist = histc(pvals, xedges);
    bar(xedges, pHist, 'k');
    xlim([-0.1 1.1]); 
    %ylim([0 35]);
    vline(mean(pvals), '-r'); vline(0.5, ':k');
    title('Frequency of P values');
    xlabel('p'); ylabel('n');
    
    subplot(3, 3, 9);
    xedges = linspace(0, 1, 20);
    rsqHist = histc(rsqs, xedges);
    bar(xedges, rsqHist, 'k');
    xlim([-0.1 1.1]);% ylim([0 35]);
    vline(mean(rsqs), '-r'); vline(0.5, ':k');
    title('Frequency of R-sq values');
    xlabel('r^2'); ylabel('n');
    
    fignum = fignum + 1;
    
    %--------------------------------------------------------------
    % FIG 2, 4, 6 - Plot the results for one example site
    h = figure(fignum);
    set(h, 'Name', ['Regression: Y = ' Yvars{i} ', X = ' Xvars{i} ...
        ', Site = ' num2str(examplesite)]);
    test = site_cl==examplesite;
    % Be sure that the variables here are used for all sites above
    eval(['ysite = ' Yvars{i} '(test);']);
    eval(['xsite = ' Xvars{i} '(test);']);
   
    subplot(2, 2, [1 2])
    plot(xsite, ysite, '.b');
    hold on;
    xrange = xlim(gca);
    [coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
    [b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
    plot(xfit, yfit,'--k');
    text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
        ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
    xlabel(Xvars{i}); ylabel(Yvars{i});
    title(examplelabel);
    
    subplot(2, 2, 3);
    xedges = linspace(min(slopes), max(slopes), 50);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    %xlim([-0.1 0.1]);
    %ylim([0 35]);
    vline(mean(slopes), '-r'); vline(0, ':k');
    title('Frequency of Slope values (all sites)');
    xlabel('Slope'); ylabel('n');
    
    subplot(2, 2, 4);
    xedges = linspace(min(pvals), max(pvals), 50);
    pvalHist = histc(pvals, xedges);
    bar(xedges, pvalHist, 'k');
    xlim([-0.1 1.1]);
    %ylim([0 35]);
    vline(mean(pvals), '-r'); vline(0.5, ':k');
    title('Frequency of P values (all sites)');
    xlabel('p value'); ylabel('n');
    
    fignum = fignum + 1;
end

% PUBLICATION PLOTS   ------------------------------------------------
% First plot regressions 1 and 2
h = figure(fignum);
examplesite = 393;%393
examplelabel = 'Chalk Creek';%'Chalk Creek';
set(h, 'Name', ['Regressions: MAST and Offset on # snowcovered days.'...
    ' Site = ' num2str(examplesite)]);
test = site_cl==examplesite;

% Example regression 1 - Set the x and y variables
eval(['ysite = ' Yvars{1} '(test);']);
eval(['xsite = ' Xvars{1} '(test);']);

subplot(2, 4, [1 2])
plot(xsite, ysite, 'ok', 'MarkerFaceColor', 'b');
hold on;
xrange = get(gca, 'Xlim');
%xrange = [190 270];% Works best for trial lake
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'units', 'normalized');
xlabel('No. snowcovered days'); ylabel('MAST (50cm)');
xlim(xrange);
title(examplelabel);

% Regression 1 - all sites
slopes = regress1dat(:, 4);
pvals = regress1dat(:, 7);
% Slope histogram
subplot(2, 4, 5);
sigtest = pvals<0.05;
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistSig = histc(slopes(sigtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistSig, 'b');
xlim([-0.14 0.14]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(2, 4, 6);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '-r'); vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
xlabel('p values'); ylabel('n');

% Example regression 2 - Set the x and y variables
eval(['ysite = ' Yvars{2} '(test);']);
eval(['xsite = ' Xvars{2} '(test);']);
subplot(2, 4, [3 4])
plot(xsite, ysite, 'ok', 'MarkerFaceColor', 'b');
hold on;
xrange = get(gca, 'Xlim');
%xrange = [190 270];% Works best for trial lake
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'units', 'normalized');
xlim(xrange);
xlabel('No. snowcovered days'); ylabel('MAST - MAT');

% Regression 2 - all sites
slopes = regress2dat(:, 4);
pvals = regress2dat(:, 7);
% Slope histogram
subplot(2, 4, 7);
sigtest = pvals<0.05;
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistSig = histc(slopes(sigtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistSig, 'b');
xlim([-0.14 0.14]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(2, 4, 8);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '-r'); vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized'); 
xlabel('p values'); ylabel('n');


% Plot regressions 3, 4, 5 ---------------------------------------
fignum = fignum + 1;
h = figure(fignum);
examplesite = 333;
examplelabel = 'Ben Lomond Trail';
set(h, 'Name', ['Regressions: Growing season vwc on SWE, melt day, and '...
    'summer rain. Site = ' num2str(examplesite)]);
test = site_cl==examplesite;

% Example regression 3 - Set the x and y variables
eval(['ysite = ' Yvars{3} '(test);']);
eval(['xsite = ' Xvars{3} '(test);']);

subplot(2, 6, [1 2])
plot(xsite, ysite, 'ok', 'MarkerFaceColor', 'b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'units', 'normalized');
xlabel('Day of snowmelt'); ylabel('50cm VWC');
ylim([0.3, 0.65]);

% Regression 3 - all sites
slopes = regress3dat(:, 4);
pvals = regress3dat(:, 7);
% Slope histogram
subplot(2, 6, 7);
sigtest = pvals<0.05;
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistSig = histc(slopes(sigtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistSig, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
%title('Frequency of slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(2, 6, 8);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '-r'); vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p values');% ylabel('n');

% Example regression 4 - Set the x and y variables
eval(['ysite = ' Yvars{4} '(test);']);
eval(['xsite = ' Xvars{4} '(test);']);
subplot(2, 6, [3 4])
plot(xsite, ysite, 'ok', 'MarkerFaceColor', 'b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'units', 'normalized');
set(gca, 'YTickLabel', '');
xlabel('Peak SWE (mm)');
ylim([0.3, 0.65]);
title(examplelabel);

% Regression 4 - all sites
slopes = regress4dat(:, 4);
pvals = regress4dat(:, 7);
% Slope histogram
subplot(2, 6, 9);
sigtest = pvals<0.05;
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistSig = histc(slopes(sigtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistSig, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
%title('Frequency of Slope values (all sites)');
xlabel('Slope');% ylabel('n');

% Pval histogram
subplot(2, 6, 10);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '-r'); vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p');% ylabel('n');

% Example regression 5 - Set the x and y variables
eval(['ysite = ' Yvars{5} '(test);']);
eval(['xsite = ' Xvars{5} '(test);']);
subplot(2, 6, [5 6])
plot(xsite, ysite, 'ok', 'MarkerFaceColor', 'b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'units', 'normalized');
set(gca, 'YTickLabel', '');
xlabel('Summer precip (mm)');
ylim([0.3, 0.65]);

% Regression 5 - all sites
slopes = regress5dat(:, 4);
pvals = regress5dat(:, 7);
% Slope histogram
subplot(2, 6, 11);
sigtest = pvals<0.05;
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistSig = histc(slopes(sigtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistSig, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
%title('Frequency of Slope values (all sites)');
xlabel('Slope');% ylabel('n');

% Pval histogram
subplot(2, 6, 12);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '-r'); vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p');% ylabel('n');

%##########################################################################
% NEW fig - same as above but only low JAS Precip sites
% Plot regressions 3, 4, 5 ---------------------------------------
fignum = fignum + 1;
h = figure(fignum);

set(h, 'Name', ['Regressions: Growing season vwc on SWE, melt day, and '...
    'summer rain. Slope/p of low precip sites']);

% Regression 3 - all sites
lowsites = regress5dat(regress5dat(:,3) < 275, 1);
lowtest = ismember(regress3dat(:,1), lowsites)
slopes = regress3dat(:, 4);
pvals = regress3dat(:, 7);
% Slope histogram
subplot(2, 3, 1);
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistLow = histc(slopes(lowtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistLow, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
title('JAS soil moisture vs melt day');
xlabel('Slope');% ylabel('n');

% Pval histogram
subplot(2, 3, 4);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
pvalHistLow = histc(pvals(lowtest), xedges);
bar(xedges, pvalHist, 'k');
hold on;
bar(xedges, pvalHistLow, 'b');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '--k'); vline(mean(pvals(lowtest)), '--b');
vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p');% ylabel('n');

% Regression 4 - all sites
lowtest = ismember(regress4dat(:,1), lowsites)
slopes = regress4dat(:, 4);
pvals = regress4dat(:, 7);
% Slope histogram
subplot(2, 3, 2);
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistLow = histc(slopes(lowtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistLow, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
title('JAS soil moisture vs peak SWE');
xlabel('Slope');% ylabel('n');

% Pval histogram
subplot(2, 3, 5);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
pvalHistLow = histc(pvals(lowtest), xedges);
bar(xedges, pvalHist, 'k');
hold on;
bar(xedges, pvalHistLow, 'b');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '--k'); vline(mean(pvals(lowtest)), '--b');
vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p');% ylabel('n');

% Regression 5 - all sites
lowtest = ismember(regress5dat(:,1), lowsites)
slopes = regress5dat(:, 4);
pvals = regress5dat(:, 7);
% Slope histogram
subplot(2, 3, 3);
xedges = linspace(min(slopes), max(slopes), 20);
slopeHist = histc(slopes, xedges);
slopeHistLow = histc(slopes(lowtest), xedges);
bar(xedges, slopeHist, 'k');
hold on;
bar(xedges, slopeHistLow, 'b');
xlim([-0.06 0.06]);
%ylim([0 35]);
vline(mean(slopes),'--k');vline(mean(slopes(sigtest)),'--b');vline(0,':k');
title('JAS soil moisture vs JAS Precip');
xlabel('Slope');% ylabel('n');

% Pval histogram
subplot(2, 3, 6);
xedges = linspace(min(pvals), max(pvals), 20);
pvalHist = histc(pvals, xedges);
pvalHistLow = histc(pvals(lowtest), xedges);
bar(xedges, pvalHist, 'k');
hold on;
bar(xedges, pvalHistLow, 'b');
xlim([-0.1 1.1]);
%ylim([0 35]);
vline(mean(pvals), '--k'); vline(mean(pvals(lowtest)), '--b');
vline(0.5, ':k');
text(0.45, 0.7, ['Combined p = ' num2str(pStouff(pvals),2)],...
    'Units', 'normalized');
%title('Frequency of P values (all sites)');
xlabel('p');% ylabel('n');