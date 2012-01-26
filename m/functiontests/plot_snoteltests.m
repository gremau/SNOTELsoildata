function plot_snoteltests(interval, siteID, m, m_test)
%
% Takes two matrices from a given site and interval. The columns in these 2
% matrices are then plotted on top of one another. This allows a view of
% what values change as a dataset is put through cleaning/filtering
% functions.
%
% args:  interval - 'daily' or 'hourly'
%        siteID - site identifier integer (only used to title figures)
%        m - first matrix (plotted first - should be "raw" data)
%        m_test - second matrix (will contain a subset or transformation of
%                 data in m

fignum=0; close all;

seq = 1:length(m{1});
seq_test = 1:length(m_test{1});

% PLOT the first and second dataset on top of one another.
% m is in red, m_test is in black
%
% Different plots for different file types
if strcmp(interval, 'daily')
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - SWE, Precip & Snow Depth']);
    
    subplot(3, 1, 1)
    plot(seq, m{4}, '.r', seq_test, m_test{4}, '.k');
    title('SWE');
    subplot(3, 1, 2)
    plot(seq, m{5}, '.r', seq_test, m_test{5}, '.k');
    title('Precip');
    subplot(3, 1, 3)
    plot(seq, m{10}, '.r', seq_test, m_test{10}, '.k');
    title('Snow depth');
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Tobs, Tmax, Tmin, Tavg']);
    
    subplot(4, 1, 1)
    plot(seq, m{6}, '.r', seq_test, m_test{6}, '.k');
    title('Tobs');
    subplot(4, 1, 2)
    plot(seq, m{7}, '.r', seq_test, m_test{7}, '.k');
    title('Tmax');
    subplot(4, 1, 3)
    plot(seq, m{8}, '.r', seq_test, m_test{8}, '.k');
    title('Tmin');
    subplot(4, 1, 4)
    plot(seq, m{9}, '.r', seq_test, m_test{9}, '.k');
    title('Tavg');
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil VWC @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(seq, m{11}, '.r', seq_test, m_test{11}, '.k');
    title('VWC -2in');
    subplot(3, 1, 2)
    plot(seq, m{12}, '.r', seq_test, m_test{12}, '.k');
    title('VWC -8in');
    subplot(3, 1, 3)
    plot(seq, m{13}, '.r', seq_test, m_test{13}, '.k');
    title('VWC -20in');
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil temp @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(seq, m{14}, '.r', seq_test, m_test{14}, '.k');
    title('Ts -2in');
    subplot(3, 1, 2)
    plot(seq, m{15}, '.r', seq_test, m_test{15}, '.k');
    title('Ts -8in');
    subplot(3, 1, 3)
    plot(seq, m{16}, '.r', seq_test, m_test{16}, '.k');
    title('Ts -20in');
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ...
        ' - Dielectric const @ 3 depths and battery voltage']);
    
    subplot(4, 1, 1)
    plot(seq, m{17}, '.r', seq_test, m_test{17}, '.k');
    title('RDC -2in');
    subplot(4, 1, 2)
    plot(seq, m{18}, '.r', seq_test, m_test{18}, '.k');
    title('RDC -8in');
    subplot(4, 1, 3)
    plot(seq, m{19}, '.r', seq_test, m_test{19}, '.k');
    title('RDC -20in');
    subplot(4, 1, 4)
    plot(seq, m{20}, '.r', seq_test, m_test{20}, '.k');
    title('BattVolt');
    
    
    
elseif strcmp(interval, 'hourly')
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil VWC @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(seq, m{4}, '.r', seq_test, m_test{4}, '.k');
    title('VWC -2in');
    subplot(3, 1, 2)
    plot(seq, m{5}, '.r', seq_test, m_test{5}, '.k');
    title('VWC -8in');
    subplot(3, 1, 3)
    plot(seq, m{6}, '.r', seq_test, m_test{6}, '.k');
    title('VWC -20in');
    
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Soil temp @ 3 depths']);
    
    subplot(3, 1, 1)
    plot(seq, m{7}, '.r', seq_test, m_test{7}, '.k');
    title('Ts -2in');
    subplot(3, 1, 2)
    plot(seq, m{8}, '.r', seq_test, m_test{8}, '.k');
    title('Ts -8in');
    subplot(3, 1, 3)
    plot(seq, m{9}, '.r', seq_test, m_test{9}, '.k');
    title('Ts -20in');
    
end