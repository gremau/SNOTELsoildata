% plot_ts_snowpack.m
%
% Plots Tsoil and Tsoil-Tair vs SWE for SNOTEL sites (and subsets 
% of SNOTELS) and fits some lines to the data.
%
% Greg Maurer - was plot_ts_vs_snowpack.m at one point
%

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

addpath('~/data/code_resources/m_common/');
addpath('~/data/code_resources/m_common/nonlinear/');
addpath('~/data/code_resources/m_common/nanstuff/');
addpath('~/data/code_resources/m_common/hline_vline/');


% Set processed data path
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));

% Import list of wasatch + uinta sites
formatstr = '%s%f%s%s';
fid = fopen([processeddatapath 'SNOTELrangelist.csv']);
wasatchUintaCell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

% Creat list of wasatch and uinta sites
wasatchTest = strcmpi(wasatchUintaCell{4}, 'WASATCH');
uintaTest =  strcmpi(wasatchUintaCell{4},'UINTA');
wasatch = wasatchUintaCell{2}(wasatchTest);
uintas = wasatchUintaCell{2}(uintaTest);
clear wasatchTest uintaTest wasatchUintaCell;

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% Get a subset of climData that corresponds with available soildata
[matchsoil, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);
% matchsoil2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

% Now assign variables
janTairMean = soilClim(:, 56);

octSWEmean = soilClim(:, 19)*25.4;
octSWEmed = soilClim(:, 20)*25.4;
octSWEsd = soilClim(:, 21)*25.4;
novSWEmean = soilClim(:, 22)*25.4;
novSWEmed = soilClim(:, 23)*25.4;
novSWEsd = soilClim(:, 24)*25.4;
decSWEmean = soilClim(:, 25)*25.4;
decSWEmed = soilClim(:, 26)*25.4;
decSWEsd = soilClim(:, 27)*25.4;
janSWEmean = soilClim(:, 28)*25.4;
janSWEmed = soilClim(:, 29)*25.4;
janSWEsd = soilClim(:, 30)*25.4;
febSWEmean = soilClim(:, 31)*25.4;
febSWEmed = soilClim(:, 32)*25.4;
febSWEsd = soilClim(:, 33)*25.4;
marSWEmean = soilClim(:, 34)*25.4;
marSWEmed = soilClim(:, 35)*25.4;
marSWEsd = soilClim(:, 36)*25.4;
aprSWEmean = soilClim(:, 37)*25.4;
aprSWEmed = soilClim(:, 38)*25.4;
aprSWEsd = soilClim(:, 39)*25.4;
maySWEmean = soilClim(:, 40)*25.4;
maySWEmed = soilClim(:, 41)*25.4;
maySWEsd = soilClim(:, 42)*25.4;
junSWEmean = soilClim(:, 43)*25.4;
junSWEmed = soilClim(:, 44)*25.4;
junSWEsd = soilClim(:, 45)*25.4;
julSWEmean = soilClim(:, 46)*25.4;
julSWEmed = soilClim(:, 47)*25.4;
julSWEsd = soilClim(:, 48)*25.4;

octTs5mean = tsData(:, 3);
octTs5sd = tsData(:, 4);
octTs20mean = tsData(:, 5);
octTs20sd = tsData(:, 6);
octTs50mean = tsData(:, 7);
octTs50sd = tsData(:, 8);
novTs5mean = tsData(:, 9);
novTs5sd = tsData(:, 10);
novTs20mean = tsData(:, 11);
novTs20sd = tsData(:, 12);
novTs50mean = tsData(:, 13);
novTs50sd = tsData(:, 14);
decTs5mean = tsData(:, 15);
decTs5sd = tsData(:, 16);
decTs20mean = tsData(:, 17);
decTs20sd = tsData(:, 18);
decTs50mean = tsData(:, 19);
decTs50sd = tsData(:, 20);
janTs5mean = tsData(:, 21);
janTs5sd = tsData(:, 22);
janTs20mean = tsData(:, 23);
janTs20sd = tsData(:, 24);
janTs50mean = tsData(:, 25);
janTs50sd = tsData(:, 26);
febTs5mean = tsData(:, 27);
febTs5sd = tsData(:, 28);
febTs20mean = tsData(:, 29);
febTs20sd = tsData(:, 30);
febTs50mean = tsData(:, 31);
febTs50sd = tsData(:, 32);
marTs5mean = tsData(:, 33);
marTs5sd = tsData(:, 34);
marTs20mean = tsData(:, 35);
marTs20sd = tsData(:, 36);
marTs50mean = tsData(:, 37);
marTs50sd = tsData(:, 38);
aprTs5mean = tsData(:, 39);
aprTs5sd = tsData(:, 40);
aprTs20mean = tsData(:, 41);
aprTs20sd = tsData(:, 42);
aprTs50mean = tsData(:, 43);
aprTs50sd = tsData(:, 44);
mayTs5mean = tsData(:, 45);
mayTs5sd = tsData(:, 46);
mayTs20mean = tsData(:, 47);
mayTs20sd = tsData(:, 48);
mayTs50mean = tsData(:, 49);
mayTs50sd = tsData(:, 50);
junTs5mean = tsData(:, 51);
junTs5sd = tsData(:, 52);
junTs20mean = tsData(:, 53);
junTs20sd = tsData(:, 54);
junTs50mean = tsData(:, 55);
junTs50sd = tsData(:, 56);
% These repeat through sept (end of wy)


% PLOTS
%----------------------------------------------------------------------
% FIG 1 - Month soil temps vs mean SWE
% Note that this can easily be set to use median SWE
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly temp vs Mean SWE - all sites/wateryears');

% Set some plotting parameters
plotorder = 1:5;
months = ['Dec';'Jan';'Feb';'Mar';'Apr'] ;
polyorder = 2;

for i=plotorder
    subplot(3, 5, i);
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    % eval(['x = ' lower(months(i,:)) 'SWEmed;']);
    eval(['y = ' lower(months(i,:)) 'Ts5mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1000]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    title(months(i,:));
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('5cm ^oC');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+5)
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    % eval(['x = ' lower(months(i,:)) 'SWEmed;']);
    eval(['y = ' lower(months(i,:)) 'Ts20mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1000]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('20cm ^oC');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+10);
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    % eval(['x = ' lower(months(i,:)) 'SWEmed;']);
    eval(['y = ' lower(months(i,:)) 'Ts50mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    if i==1
        ylabel('50cm ^oC');
    elseif i==3
        xlabel('Mean SWE');
        % xlabel('Median SWE');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
end


%----------------------------------------------------------------------
% FIG 2 - Month soil temps vs mean SWE - use bounded exponential fit
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly temp vs Median SWE - all sites/wateryears');

% Set some plotting parameters
plotorder = 1:5;
months = ['Dec';'Jan';'Feb';'Mar';'Apr'];
% Set up non-linear fits - Bounded exponential
nlfunc = inline('b(1)*(1 - b(2)*exp(-b(3).*x))', 'b', 'x');
beta_init = [1.0 .2 0.01];

for i=plotorder
    subplot(3, 5, i)
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts5mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    
    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k');
    % Mean Squared error and r-squared values
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 1.5);
    
    xlim([0, 1500]); ylim([-10, 10]);
    title(months(i,:));
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('5cm ^oC');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+5)
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts20mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;

    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k');
    % Mean Squared error and r-squared values
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 1.5);
    
    xlim([0, 1500]); ylim([-10, 10]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('20cm ^oC');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+10)
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts50mean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    
    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k');
    % Mean Squared error and r-squared values
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 1.5);
    

    xlim([0, 1500]); ylim([-10, 10]);
    if i==1
        ylabel('50cm ^oC');
    elseif i==3
        xlabel('Mean SWE');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end

end

%--------------------------------------------------------
% FIG 3 - Subset of above - this was in AGU 2011 poster
fignum = fignum+1;    
h = figure(fignum);

subplot 121;
plot(janSWEmean, janTs5mean, 'ob');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
title('December soil temperatures');
legend('5cm one-month mean');

% Mean month soil temp by SWE - 20cm
subplot 122;
plot(janSWEmean, janTs20mean, 'ok');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
legend('20cm one-month mean');
%title([monthlabel ' 20cm temp vs ' monthlabel ' SWE']);

%--------------------------------------------------------
% FIG 4 - Like above - but try some new curve fitting
% This was where I figured out the best curve, several are possible here
fignum = fignum+1;    
h = figure(fignum);

%Maybe this is the right line to use
% y = arrayfun(@(x) tmax*(1-exp(-x/tmax)) * exp(-50/tmax), 1:100)
% or func = inline('b1*(1-exp(-x/b1)) * exp(-50/b1)')

subplot 121;
plot(janSWEmean, janTs5mean, 'ob');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
title('December soil temperatures');
legend('5cm one-month mean');
hold on

% Bounded exponential
boundedexp = inline('b(1)*(1 - b(2)*exp(-b(3).*x))', 'b', 'x');
beta_init = [1.5 2 0.01];
[coeffs, r, ~] = nlinfit(janSWEmean, janTs5mean, boundedexp, beta_init);
rmse = sqrt(nansum(r.^2)/(length(janSWEmean)-1));
disp(['Bounded exponential: ' num2str(coeffs) ' RMSE = ' num2str(rmse)]);
plot(0:1500, boundedexp(coeffs, 0:1500), '--r');

% Logistic fit -> K./(1+exp(-r*(t-t0)));
logist = inline('(b(1)./(1 + exp(-b(2).*x))) - b(3)', 'b', 'x');
beta_init = [1.5 0.01 1];
[coeffs, r, ~] = nlinfit(janSWEmean, janTs5mean, logist, beta_init);
rmse = sqrt(nansum(r.^2)/(length(janSWEmean)-1));
disp(['Logistic function: ' num2str(coeffs) ' RMSE = ' num2str(rmse)]);
plot(0:1500, logist(coeffs, 0:1500), '-r');

% Mean month soil temp by SWE - 20cm
subplot 122;
plot(janSWEmean, janTs20mean, 'ok');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
legend('20cm one-month mean');
%title([monthlabel ' 20cm temp vs ' monthlabel ' SWE']);
hold on

% Asymptotic fit
% asymp = inline('a + b/x'); % asymptotic line
%plot(1:1500, arrayfun(@(x) asymp(1, -11, x), 1:1500), ':k');

% Bounded exponential
boundedexp = inline('b(1)*(1 - b(2)*exp(-b(3).*x))', 'b', 'x');
beta_init = [1.5 2 0.01];
[coeffs, r, ~] = nlinfit(janSWEmean, janTs20mean, boundedexp, beta_init);
rmse = sqrt(nansum(r.^2)/(length(janSWEmean)-1));
disp(['Bounded exponential: ' num2str(coeffs) ' RMSE = ' num2str(rmse)]);
plot(0:1500, boundedexp(coeffs, 0:1500), '--r');

% Logistic fit -> K./(1+exp(-r*(t-t0)));
logist = inline('(b(1)./(1 + exp(-b(2).*x))) - b(3)', 'b', 'x');
beta_init = [1.5 0.01 1];
[coeffs, r, ~] = nlinfit(janSWEmean, janTs20mean, logist, beta_init);
rmse = sqrt(nansum(r.^2)/(length(janSWEmean)-1));
disp(['Logistic function: ' num2str(coeffs) ' RMSE = ' num2str(rmse)]);
plot(0:1500, logist(coeffs, 0:1500), '-r');

%----------------------------------------------------------------------
% FIG 5 - Plot a small subset of months/depths and use better fitlines
fignum = fignum+1;    
h = figure('position',[100 0 1000 800],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(h, 'Name', 'Mean monthly temp vs Mean SWE - all sites/wateryears');
set(h, 'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% Set some plotting parameters
plotorder = [1, 2];
months = ['Dec';'Jan'] ;
polyorder = 1;

% Set up non-linear fits
% Bounded exponential
nlfunc = inline('b(1)*(1 - b(2)*exp(-b(3).*x))', 'b', 'x');
beta_init = [1.5 2 0.01];
% Logistic fit -> K./(1+exp(-r*(t-t0)));
% nlfunc = inline('(b(1)./(1 + exp(-b(2).*x))) - b(3)', 'b', 'x');
% beta_init = [1.5 0.01 1];

for i=plotorder
    subplot(2, 2, i);
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts5mean;']);
    plot(x, y, '.', 'Color', [0.5,0.5,0.5], 'markersize', 10);
    hold on;
    
    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k');
    % Mean Squared error and r-squared values
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 1.5);
    
    set(gca, 'position', [0.90 0.91 1.15 1.2] .* get(gca, 'position'));
    xlim([0, 700]); ylim([-6, 6]);
    title(months(i,:), 'Fontsize', 20, 'Fontangle', 'italic');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('5cm T_{soil} (^oC)');
        set(gca,'Ytick', [-4,-2,0,2,4,6]);
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(2, 2, i+2)
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts20mean;']);
    plot(x, y, '.', 'Color', [0.5,0.5,0.5], 'markersize', 10);
    hold on;
    
    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k');
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 1.5);
    
    set(gca, 'position', [0.90 0.95 1.15 1.2] .* get(gca, 'position'));
    xlim([0, 700]); ylim([-6, 6]);
    xlabel('Mean SWE (mm)');
    if i==1
        ylabel('20cm T_{soil} (^oC)');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
end
figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figG.eps']) 

% FIG 6 - Plot a small subset of months/depths and use better fitlines
fignum = fignum+1;    
h = figure('position',[100 0 1100 450],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(h, 'Name', 'Mean monthly temp vs Mean SWE - all sites/wateryears');
set(h, 'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% Set some plotting parameters
plotorder = [1, 2];
months = ['Dec';'Jan'] ;
polyorder = 1;

% Set up non-linear fits
% Bounded exponential
nlfunc = inline('b(1)*(1 - b(2)*exp(-b(3).*x))', 'b', 'x');
beta_init = [1.5 2 0.01];
% Logistic fit -> K./(1+exp(-r*(t-t0)));
% nlfunc = inline('(b(1)./(1 + exp(-b(2).*x))) - b(3)', 'b', 'x');
% beta_init = [1.5 0.01 1];

for i=plotorder
    subplot(1, 2, i);
    eval(['x = ' lower(months(i,:)) 'SWEmean;']);
    eval(['y = ' lower(months(i,:)) 'Ts5mean;']);
    plot(x, y, '.', 'Color', [0.5 0.5 0.5], 'markersize', 10);
    hold on;
    
    % Non-linear fit
    [coeffs, r, ~] = nlinfit(x, y, nlfunc, beta_init);
    rmse = sqrt(nansum(r.^2)/(length(x)-1));
    rsq = 1 - (nansum(r.^2)/((length(y)-1) * nanvar(y)));
    plot(0:1500, nlfunc(coeffs, 0:1500), '-k', 'Linewidth', 1.5);
    % Mean Squared error and r-squared values
    text(0.6, 0.8, ['RMSE = ' num2str(rmse, 2)], 'units', 'normalized');
    text(0.6, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized');
    % Plot horizontal zero line
    l = hline(0, ':k');
    set(l, 'linewidth', 2);
    
    set(gca, 'position', [0.90 1 1.15 1] .* get(gca, 'position'));
    xlim([0, 700]); ylim([-6, 6]);
    %title(months(i,:), 'Fontsize', 20, 'Fontangle', 'italic');
    %set(gca, 'XTickLabel', '');
    if i==1
        ylabel('Mean monthly T_{soil} (^oC)');
        set(gca,'Ytick', [-4,-2,0,2,4,6]);
        xlabel('Mean December SWE (mm)');
    elseif i>1
        set(gca, 'YTickLabel', '');
        xlabel('Mean January SWE (mm)');
    end
end
figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figG2.eps']) 

%------------------------------------------------------
% FIG 7 - Soil temp and offset - plot for wasatch and uinta sites 
uintaTest = ismember(soilClim(:, 1), uintas);
wasatchTest = ismember(soilClim(:, 1), wasatch);

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean January Ts vs snowpack - Wasatch & Uintas');
% Mean month soil temp by SWE - 5cm
subplot 321;
plot(janSWEmean(wasatchTest), janTs5mean(wasatchTest), 'om');
hold on;
plot(janSWEmean(uintaTest), janTs5mean(uintaTest), 'ob');
bothTest = wasatchTest | uintaTest;
x = janSWEmean(bothTest);
y = janTs5mean(bothTest);
yNoNan = ~isnan(y);
coefficients = polyfit(x(yNoNan), y(yNoNan), 2);
pnFit = polyval(coefficients, (0:1:700));
plot((0:1:700), pnFit,':k');
xlabel('Mean SWE');
ylabel('5cm Ts (^oC)');
ylim([-4, 4]);
legend('Wasatch mtns', 'Uinta mtns');
title('January Ts vs January SWE');

% Mean month soil temp by SWE - 20cm
subplot 323;
plot(janSWEmean(wasatchTest), janTs20mean(wasatchTest), 'om');
hold on;
plot(janSWEmean(uintaTest), janTs20mean(uintaTest), 'ob');
bothTest = wasatchTest | uintaTest;
x = janSWEmean(bothTest);
y = janTs20mean(bothTest);
yNoNan = ~isnan(y);
coefficients = polyfit(x(yNoNan), y(yNoNan), 2);
pnFit = polyval(coefficients, (0:1:700));
plot((0:1:700), pnFit,':k');
ylim([-4, 4]);
ylabel('20cm Ts (^oC)');

% Mean month soil temp by snow depth - 20cm
subplot 325;
plot(janSWEmean(wasatchTest), janTs50mean(wasatchTest), 'om');
hold on;
plot(janSWEmean(uintaTest), janTs50mean(uintaTest), 'ob');
bothTest = wasatchTest | uintaTest;
x = janSWEmean(bothTest);
y = janTs50mean(bothTest);
yNoNan = ~isnan(y);
coefficients = polyfit(x(yNoNan), y(yNoNan), 2);
pnFit = polyval(coefficients, (0:1:700));
plot((0:1:700), pnFit,':k');
ylim([-4, 4]);
xlabel('Mean SWE');
ylabel('50cm Ts (^oC)');

% Mean month soil temps vs air temps
offset5 = janTs5mean - janTairMean;
offset20 = janTs20mean - janTairMean;
offset50 = janTs50mean - janTairMean;

subplot 322;
plot(janSWEmean(wasatchTest), offset5(wasatchTest), 'om');
hold on
plot(janSWEmean(uintaTest), offset5(uintaTest), 'ob');
title('Tsoil - Tair vs snowpack');

subplot 324;
plot(janSWEmean(wasatchTest), offset20(wasatchTest), 'om');
hold on
plot(janSWEmean(uintaTest), offset20(uintaTest), 'ob');

subplot 326;
plot(janSWEmean(wasatchTest), offset50(wasatchTest), 'om');
hold on
plot(janSWEmean(uintaTest), offset50(uintaTest), 'ob');
xlabel('Mean SWE');


junk = 99;
