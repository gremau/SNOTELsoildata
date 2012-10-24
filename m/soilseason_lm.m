%soilseason_lm.m

close all;
clear all;
fignum = 1;

% Add any needed tools
addpath('/home/greg/data/code_resources/m_common/'); 
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

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


% Set regression variables:

% Yvars = {'mast20cm','mast20cm','mast20cm'};
% Xvars = {'maat', 'totaldaysSC', 'maxswe'};

% Yvars = {'julVWC20mean','augVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'meltdoy', 'meltdoy'};

% Yvars = {'jasVWC20mean','jasVWC20mean','jasVWC20mean'};
% Xvars = {'meltdoy', 'maxswe', 'JASprecip'};

% Yvars = {'janVWC20mean','febVWC20mean','marVWC20mean'};
% Xvars = {'preonsetVWC20', 'preonsetVWC20', 'preonsetVWC20'};

% diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
% Yvars = {'diff5cm','diff20cm','diff50cm'};
% Xvars = {'onsetdoy', 'onsetdoy', 'onsetdoy'};

diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
Yvars = {'diff5cm','diff5cm','diff5cm'};
Xvars = {'maat', 'totaldaysSC', 'maxswe'};


for i = 1:length(Yvars)
    % Do a regression of MAST on snowcovered days at each site
    slopes = [];
    yints = [];
    rsqs = [];
    xmeans = [];
    elevs = [];
    for j = 1:length(sites)
        getsite = site_cl==sites(j);
        eval(['y = ' Yvars{i} '(getsite);']);
        eval(['x = ' Xvars{i} '(getsite);']);
        % Only use sites with at least 3 datapoints after nan removal
        nantest = isnan(x) | isnan(y);
        if (sum(getsite) - sum(nantest)) > 2
            [coefficients, rsq, ~, ~] = fitline(x, y, 1, [0, 1]);
            slopes = [slopes; coefficients(1)];
            yints = [yints; coefficients(2)];
            rsqs = [rsqs; rsq];
            xmeans = [xmeans; nanmean(x)];
            elevs = [elevs; mean(elev(getsite))];
        end
    end
    
    
    % Plot the results of the regression analysis
    h = figure(fignum);
    set(h, 'Name', ['Regressions: Y = ' Yvars{i} ', X = ' Xvars{i}]);
    
    subplot(3, 3, 1);
    plot(xmeans, slopes, 'ok');
    title('Slope');
    hline(mean(slopes), ':k')
    xlabel('Site mean of x');
    
    subplot(3, 3, 2);
    plot(xmeans, yints, 'ob');
    hline(mean(yints), ':b')
    title('Y - int');
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
    plot(elevs, yints, '.b');
    xlabel('Elev (m)');
    
    subplot(3, 3, 6);
    plot(elevs, rsqs, '.r');
    xlabel('Elev (m)');
    
    subplot (3, 3, 7)
    xedges = linspace(min(slopes), max(slopes), 50);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    vline(mean(slopes), ':k');
    %xlim([-0.1 0.1]); 
    ylim([0 35]);
    title('Frequency of Slope values');
    xlabel('Slope'); ylabel('n');
        
    subplot (3, 3, 8)
    xedges = linspace(min(slopes), max(slopes), 50);
    slopeHist = histc(slopes, xedges);
    bar(xedges, slopeHist, 'k');
    vline(mean(slopes), ':k');
    %xlim([-0.1 0.1]); 
    ylim([0 35]);
    title('Frequency of Slope values');
    xlabel('Slope'); ylabel('n');
    
    subplot (3, 3, 9)
    xedges = linspace(0, 1, 50);
    rsqHist = histc(rsqs, xedges);
    bar(xedges, rsqHist, 'k');
    vline(mean(rsqs), ':k');
    xlim([-0.1 1.1]);% ylim([0 35]);
    title('Frequency of R-sq values');
    xlabel('r^2'); ylabel('n');
    
    fignum = fignum + 1;
end