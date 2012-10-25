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
examplesite = 828;
examplelabel = 'Trial Lake';
% examplesite = 393;
% examplelabel = 'Chalk Creek';
% examplesite = 332;
% examplelabel = 'Ben Lomond Peak';
% examplesite = 333;
% examplelabel = 'Ben Lomond Trail';

% Set processed data path
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt']);

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
JASprecip = soilClim(:, 13);
maat = soilClim(:, 74);
maat_sd = soilClim(:, 75);

preonsetTair = soilClim(:, 76);
preonsetTairSd = soilClim(:, 77);
premeltTair = soilClim(:, 78);
premeltTairSd = soilClim(:, 79);
postmeltTair = soilClim(:, 80);
postmeltTairSd = soilClim(:, 81);
elev = soilClim(:, 82);
lat = soilClim(:, 83);
lon = soilClim(:, 84);
ltMeanSWE = soilClim(:, 85);
ltMeanPrecip = soilClim(:, 86);

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
diff20cm = mast20cm-maat;
Yvars = {'mast20cm' 'diff20cm' 'jasVWC50mean','jasVWC50mean','jasVWC50mean'};
Xvars = {'totaldaysSC' 'totaldaysSC' 'meltdoy', 'maxswe', 'JASprecip'};

% EXPLORATORY ANALYSIS  ------------------------------------------------
% Set regression variables:

% Yvars = {'mast20cm','mast20cm','mast20cm'};
% Xvars = {'maat', 'totaldaysSC', 'snowduration'};

% Yvars = {'julVWC20mean','augVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'meltdoy', 'meltdoy'};

% Yvars = {'jasVWC20mean','jasVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'maxswe', 'JASprecip'};

% Yvars = {'janVWC20mean','febVWC20mean','marVWC20mean'};
% Xvars = {'preonsetVWC20', 'preonsetVWC20', 'preonsetVWC20'};

% diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
% Yvars = {'diff5cm','diff20cm','diff50cm'};
% Xvars = {'onsetdoy', 'onsetdoy', 'onsetdoy'};

% diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
% Yvars = {'diff20cm','diff20cm','diff20cm'};
% Xvars = {'maat', 'totaldaysSC', 'snowduration'};

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
        if (sum(getsite) - sum(nantest)) > 2
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
    xedges = linspace(min(slopes), max(slopes), 50);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    %xlim([-0.1 0.1]); 
    %ylim([0 35]);
    vline(mean(slopes), '-b');
    title('Frequency of Slope values');
    xlabel('Slope'); ylabel('n');
        
    subplot(3, 3, 8);
    xedges = linspace(min(pvals), max(pvals), 50);
    pHist = histc(pvals, xedges);
    bar(xedges, pHist, 'k');
    %xlim([-0.1 0.1]); 
    %ylim([0 35]);
    vline(mean(pvals), '-b');
    title('Frequency of P values');
    xlabel('p'); ylabel('n');
    
    subplot(3, 3, 9);
    xedges = linspace(0, 1, 50);
    rsqHist = histc(rsqs, xedges);
    bar(xedges, rsqHist, 'k');
    xlim([-0.1 1.1]);% ylim([0 35]);
    vline(mean(rsqs), '-b');
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
   
    subplot(2, 1, 1)
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
    
    subplot(2, 1, 2);
    xedges = linspace(min(slopes), max(slopes), 50);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    %xlim([-0.1 0.1]);
    %ylim([0 35]);
    vline(mean(slopes), '-b');
    title('Frequency of Slope values (all sites)');
    xlabel('Slope'); ylabel('n');
    
    fignum = fignum + 1;
end

% PUBLICATION PLOTS   ------------------------------------------------
% First plot regressions 1 and 2
h = figure(fignum);
examplesite = 828;
examplelabel = 'Trial Lake';
set(h, 'Name', ['Regressions: MAST and Offset on # snowcovered days.'...
    ' Site = ' num2str(examplesite)]);
test = site_cl==examplesite;

% Example regression 1 - Set the x and y variables
eval(['ysite = ' Yvars{1} '(test);']);
eval(['xsite = ' Xvars{1} '(test);']);

subplot(3, 4, [1 2 5 6])
plot(xsite, ysite, '.b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
xlabel(Xvars{1}); ylabel(Yvars{1});
title(examplelabel);

% Regression 1 - all sites
% Slope histogram
subplot(3, 4, 9);
slopes = regress1dat(:, 4);
xedges = linspace(min(slopes), max(slopes), 50);
slopeHist = histc(slopes, xedges);
bar(xedges, slopeHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(slopes), '-b');
title('Frequency of slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(3, 4, 10);
pvals = regress1dat(:, 7);
xedges = linspace(min(pvals), max(pvals), 50);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(pvals), '-b');
title('Frequency of P values (all sites)');
xlabel('p'); ylabel('n');

% Example regression 2 - Set the x and y variables
eval(['ysite = ' Yvars{2} '(test);']);
eval(['xsite = ' Xvars{2} '(test);']);
subplot(3, 4, [3 4 7 8])
plot(xsite, ysite, '.b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
xlabel(Xvars{2}); ylabel(Yvars{2});
title(examplelabel);

% Regression 2 - all sites
% Slope histogram
subplot(3, 4, 11);
slopes = regress2dat(:, 4);
xedges = linspace(min(slopes), max(slopes), 50);
slopeHist = histc(slopes, xedges);
bar(xedges, slopeHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(slopes), '-b');
title('Frequency of Slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(3, 4, 12);
pvals = regress2dat(:, 7);
xedges = linspace(min(pvals), max(pvals), 50);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(pvals), '-b');
title('Frequency of P values (all sites)');
xlabel('p'); ylabel('n');

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

subplot(3, 6, [1 2 7 8])
plot(xsite, ysite, '.b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
xlabel(Xvars{3}); ylabel(Yvars{3});
title(examplelabel);

% Regression 3 - all sites
% Slope histogram
subplot(3, 6, 13);
slopes = regress3dat(:, 4);
xedges = linspace(min(slopes), max(slopes), 50);
slopeHist = histc(slopes, xedges);
bar(xedges, slopeHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(slopes), '-b');
title('Frequency of slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(3, 6, 14);
pvals = regress3dat(:, 7);
xedges = linspace(min(pvals), max(pvals), 50);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(pvals), '-b');
title('Frequency of P values (all sites)');
xlabel('p'); ylabel('n');

% Example regression 4 - Set the x and y variables
eval(['ysite = ' Yvars{4} '(test);']);
eval(['xsite = ' Xvars{4} '(test);']);
subplot(3, 6, [3 4 9 10])
plot(xsite, ysite, '.b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
xlabel(Xvars{4}); ylabel(Yvars{4});
title(examplelabel);

% Regression 4 - all sites
% Slope histogram
subplot(3, 6, 15);
slopes = regress4dat(:, 4);
xedges = linspace(min(slopes), max(slopes), 50);
slopeHist = histc(slopes, xedges);
bar(xedges, slopeHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(slopes), '-b');
title('Frequency of Slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(3, 6, 16);
pvals = regress4dat(:, 7);
xedges = linspace(min(pvals), max(pvals), 50);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(pvals), '-b');
title('Frequency of P values (all sites)');
xlabel('p'); ylabel('n');

% Example regression 5 - Set the x and y variables
eval(['ysite = ' Yvars{5} '(test);']);
eval(['xsite = ' Xvars{5} '(test);']);
subplot(3, 6, [5 6 11 12])
plot(xsite, ysite, '.b');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = shregress(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)]); % r^2 & p
xlabel(Xvars{5}); ylabel(Yvars{5});
title(examplelabel);

% Regression 5 - all sites
% Slope histogram
subplot(3, 6, 17);
slopes = regress5dat(:, 4);
xedges = linspace(min(slopes), max(slopes), 50);
slopeHist = histc(slopes, xedges);
bar(xedges, slopeHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(slopes), '-b');
title('Frequency of Slope values (all sites)');
xlabel('Slope'); ylabel('n');

% Pval histogram
subplot(3, 6, 18);
pvals = regress5dat(:, 7);
xedges = linspace(min(pvals), max(pvals), 50);
pvalHist = histc(pvals, xedges);
bar(xedges, pvalHist, 'k');
%xlim([-0.1 0.1]);
%ylim([0 35]);
vline(mean(pvals), '-b');
title('Frequency of P values (all sites)');
xlabel('p'); ylabel('n');