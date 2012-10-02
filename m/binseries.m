% binseries.m
% Binning function
% Bin the vwc datapoints and then generate a mean for each bin
% Data:
% x, y1, y2 = x data, ydata and ydata2 (usually std of y1)
% Binning parameters:
% topEdge = 300; % define upper x range of bins
% botEdge = 90; % define lower x range of bins
% numBins = 15; % define number of bins


function [binMean1, binMean2] = binseries(x, y1, y2, topEdge, botEdge, numBins)
    addpath('/home/greg/data/code_resources/m_common/nanstuff/');
    binEdges = linspace(botEdge, topEdge, numBins+1);
    [h,whichBin] = histc(x, binEdges);
    for i = 1:numBins
        flagBinMembers = (whichBin == i);
        binMembers1     = y1(flagBinMembers);
        binMembers2     = y2(flagBinMembers);
        binMean1(i)     = nanmean(binMembers1);
        %binStd1(i)      = nanstd(binMembers1);
        binMean2(i)     = nanmean(binMembers2);
        %binStd2(i)      = nanstd(binMembers2);
    end
end