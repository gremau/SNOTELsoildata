function plot_snoteltests(interval, siteID, mraw, mtest)
%
% Takes two matrices from a given site and interval. The columns in these 2
% matrices are then plotted on top of one another. This allows a view of
% what values change as a dataset is put through cleaning/filtering
% functions.
%
% args:  interval - 'daily' or 'hourly'
%        siteID - site identifier integer (only used to title figures)
%        mraw - first matrix (plotted first - should be "raw" data)
%        mtest - second matrix (will contain a subset or transformation of
%                 data in mraw

fignum=0; close all;

% Plot using the record sequence number - DEPRECATED
% seq = 1:length(mraw{1});
% seq_test = 1:length(mtest{1});

% Create datetime vector to use in plotting
if strcmp(interval, 'hourly')
    tvec = datenum(strcat(mraw{2}, mraw{3}), 'yyyy-mm-ddHH:MM');
    tvec_test = datenum(strcat(mtest{2}, mtest{3}), 'yyyy-mm-ddHH:MM');
elseif strcmp(interval, 'daily')
    tvec = datenum(mraw{2}, 'yyyy-mm-dd');
    tvec_test = datenum(mtest{2}, 'yyyy-mm-dd');
end

% PLOT the first and second dataset on top of one another.
% m is in red, mtest is in black
%
% Different plots for different file types
if strcmp(interval, 'daily')
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - SWE, Precip & Snow Depth']);
    
    subplot(3, 1, 1)
    plot(tvec, mraw{4}, '.r', tvec_test, mtest{4}, '.k');
    title('SWE'); datetick();
    subplot(3, 1, 2)
    plot(tvec, mraw{5}, '.r', tvec_test, mtest{5}, '.k');
    title('Precip'); datetick();
    subplot(3, 1, 3)
    plot(tvec, mraw{10}, '.r', tvec_test, mtest{10}, '.k');
    title('Snow depth'); datetick();
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Tobs, Tmax, Tmin, Tavg']);
    
    subplot(4, 1, 1)
    plot(tvec, mraw{6}, '.r', tvec_test, mtest{6}, '.k');
    title('Tobs'); datetick();
    subplot(4, 1, 2)
    plot(tvec, mraw{7}, '.r', tvec_test, mtest{7}, '.k');
    title('Tmax'); datetick();
    subplot(4, 1, 3)
    plot(tvec, mraw{8}, '.r', tvec_test, mtest{8}, '.k');
    title('Tmin'); datetick();
    subplot(4, 1, 4)
    plot(tvec, mraw{9}, '.r', tvec_test, mtest{9}, '.k');
    title('Tavg'); datetick();
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil VWC @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(tvec, mraw{11}, '.r', tvec_test, mtest{11}, '.k');
    title('VWC -2in'); datetick();
    subplot(3, 1, 2)
    plot(tvec, mraw{12}, '.r', tvec_test, mtest{12}, '.k');
    title('VWC -8in'); datetick();
    subplot(3, 1, 3)
    plot(tvec, mraw{13}, '.r', tvec_test, mtest{13}, '.k');
    title('VWC -20in'); datetick();
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil temp @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(tvec, mraw{14}, '.r', tvec_test, mtest{14}, '.k');
    title('Ts -2in'); datetick();
    subplot(3, 1, 2)
    plot(tvec, mraw{15}, '.r', tvec_test, mtest{15}, '.k');
    title('Ts -8in'); datetick();
    subplot(3, 1, 3)
    plot(tvec, mraw{16}, '.r', tvec_test, mtest{16}, '.k');
    title('Ts -20in'); datetick();
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ...
        ' - Dielectric const @ 3 depths and battery voltage']);
    
    subplot(4, 1, 1)
    plot(tvec, mraw{17}, '.r', tvec_test, mtest{17}, '.k');
    title('RDC -2in'); datetick();
    subplot(4, 1, 2)
    plot(tvec, mraw{18}, '.r', tvec_test, mtest{18}, '.k');
    title('RDC -8in'); datetick();
    subplot(4, 1, 3)
    plot(tvec, mraw{19}, '.r', tvec_test, mtest{19}, '.k');
    title('RDC -20in'); datetick();
    subplot(4, 1, 4)
    plot(tvec, mraw{20}, '.r', tvec_test, mtest{20}, '.k');
    title('BattVolt'); datetick();
    
    
    
elseif strcmp(interval, 'hourly')
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil VWC @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(tvec, mraw{4}, '.r', tvec_test, mtest{4}, '.k');
    title('VWC -2in'); datetick();
    subplot(3, 1, 2)
    plot(tvec, mraw{5}, '.r', tvec_test, mtest{5}, '.k');
    title('VWC -8in'); datetick();
    subplot(3, 1, 3);
    plot(tvec, mraw{6}, '.r', tvec_test, mtest{6}, '.k');
    title('VWC -20in'); datetick();

    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil temp @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(tvec, mraw{7}, '.r', tvec_test, mtest{7}, '.k');
    title('Ts -2in'); datetick();
    subplot(3, 1, 2)
    plot(tvec, mraw{8}, '.r', tvec_test, mtest{8}, '.k');
    title('Ts -8in'); datetick();
    subplot(3, 1, 3)
    plot(tvec, mraw{9}, '.r', tvec_test, mtest{9}, '.k');
    title('Ts -20in'); datetick();
    
end