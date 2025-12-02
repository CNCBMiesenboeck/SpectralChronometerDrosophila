%% To reproduce the MBON timing and stats shown in fig 6
%
% % % Copyright (C) 2025 Clifford Talbot
% % % This program is free software: you can redistribute it and/or modify
% % % it under the terms of the GNU General Public License as published by
% % % the Free Software Foundation, either version 3 of the License, or
% % % (at your option) any later version.
% % % 
% % % This program is distributed in the hope that it will be useful,
% % % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% % % GNU General Public License for more details.
%
% This script performs the stats and response timing analysis on the MBON
% e-phys recordings presented in figure 6. It produces a figure matching
% figure 6B and outputs the statistics to the command window.

% clear the workspace, screen and load the data
clear
clc
load('data\MBON_spikingData.mat');
% MBON_spikeTimes is a 4-by-20 cell array with the action potential times
% of the 20 MBONs in response to the four odour presentations:
%   row 1: early control
%   row 2: early test
%   row 3: late control
%   row 4: late test
% "early" refers to the odour reinforced between 3-4 seconds after onset. 
% "late" refers to the odour reinforced between 7-8 seconds after onset.
% "control" refers to the recordings performed before reinforcement.
% "test" refers to the recordings performed after reinforcement.

% the psth binning width:
binWidth = 0.5; % In the publication, we used 0.5.
% time range for the baseline subtraction:
baselineRange = [-5,-0.5]; % publication range: [-5,-0.5]
% time range for the divergence calculations (the odour presentation time)
divergenceTimeRange = [0,10]; % publication: [0,10]
% Options for the noise stabilisation. See Pouzat as referenced in the
% publication for details. This makes virtually no difference for the count
% rates in the MBON measurements.
stabilisationInd = 1; % 1, 2, or 3. In the publication, we use 1.
% the time range to plot the figures:
plotRange = [-3,13];
% the reinforcement times (over which the stats are performed):
reinforcementTimes = [3,4; 7,8];

% obtain the first and last spike time of all cells and then determine the
% histogram bin edges and centres (for the PSTHs, and divergence
% computation):
tempMin = cellfun(@min, MBON_spikeTimes, 'UniformOutput', false);
tempMax = cellfun(@max, MBON_spikeTimes, 'UniformOutput', false);
minSpikeTime = min(cell2mat(tempMin(~cellfun(@isempty,tempMin))))/1000;
maxSpikeTime = max(cell2mat(tempMax(~cellfun(@isempty,tempMax))))/1000;
timeEdges = binWidth*(floor(minSpikeTime/binWidth):ceil(maxSpikeTime/binWidth));
timeCentres = timeEdges(1:(end-1))+binWidth/2;
nBins = numel(timeEdges)-1;
divergenceBins = find(timeCentres>divergenceTimeRange(1)&timeCentres<divergenceTimeRange(2));
nBinsDivergence = numel(divergenceBins);
% which PSTH bins to use for plotting:
plotBins = find(timeCentres>plotRange(1)&timeCentres<plotRange(2));

% some simple helper functions
% to count the number of spikes in a time range:
hCountTime = @(x,t1,t2) sum(x>(t1*1000) & x<=(t2*1000));
% to count the number of spikes in a psth bin (pp indexes odour 
% presentations and cc indexes cell number - see description of
% MBON_spikeTimes. bb indexes the PSTH time bin):
hPSTHHelper = @(pp,cc,bb) hCountTime(MBON_spikeTimes{pp,cc},timeEdges(bb),timeEdges(bb+1));
% to compute the baseline spike rate:
hBaseline = @(x) round(hCountTime(x,baselineRange(1),baselineRange(2)) * ...
    (binWidth/(baselineRange(2)-baselineRange(1))));
% to count the number of spikes in the reinforcement intervals (for stats):
hStatsCountsHelper = @(ii,kk) cellfun(@(x) ...
    hCountTime(x,reinforcementTimes(ii,1),reinforcementTimes(ii,2)),MBON_spikeTimes(kk,:));
% to perform the noise stabilisation used in the random walk (divergence)
% calculation:
switch stabilisationInd
    case 1
        hStabilised = @(psthIn) 2*abs(sqrt(psthIn+0.25));
    case 2
        hStabilised = @(psthIn) abs(sqrt(psthIn))+abs(sqrt(psthIn+1));
    case 3
        hStabilised = @(psthIn) 2*abs(sqrt(psthIn+0.375));
end
rt2=sqrt(2);
hStabilisedDifference = @(test,ctrl) sum(hStabilised(ctrl)-hStabilised(test))/rt2;

% compute the baseline subtracted PSTHs:
nCells = size(MBON_spikeTimes,2);
[earlyCtrlPSTH, earlyTestPSTH, lateCtrlPSTH, lateTestPSTH] = ...
    deal(zeros(nCells,nBins));
[earlyCtrlBL, earlyTestBL, lateCtrlBL, lateTestBL] = deal(zeros(nCells, 1));
for cc = 1:nCells
    earlyCtrlBL(cc) = hBaseline(MBON_spikeTimes{1,cc});
    earlyTestBL(cc) = hBaseline(MBON_spikeTimes{2,cc});
    lateCtrlBL(cc) = hBaseline(MBON_spikeTimes{3,cc});
    lateTestBL(cc) = hBaseline(MBON_spikeTimes{4,cc});
    for bb = 1:nBins
        earlyCtrlPSTH(cc,bb) = hPSTHHelper(1,cc,bb)-earlyCtrlBL(cc);
        earlyTestPSTH(cc,bb) = hPSTHHelper(2,cc,bb)-earlyTestBL(cc);
        lateCtrlPSTH(cc,bb) = hPSTHHelper(3,cc,bb)-lateCtrlBL(cc);
        lateTestPSTH(cc,bb) = hPSTHHelper(4,cc,bb)-lateTestBL(cc);
    end
end

% preallocate some arrays for the divergence calculation:
[cumdiff, normdiff] = deal(zeros(nBinsDivergence));

% for output:
hFig = figure();
fprintf('\n');

% loop over odours to perform the random walk divergence calcuation:
for ee = 1:2
    % to extract the data from the correct odour presentations:
    switch ee
        case 1 % for 'early'
            odourTitle = 'Early';
            % for the calculations:
            ctrl = sum(earlyCtrlPSTH(:, divergenceBins));
            test = sum(earlyTestPSTH(:, divergenceBins));
            % for the plotting:
            plotCtrl = mean(earlyCtrlPSTH(:, plotBins))/binWidth;
            plotTest = mean(earlyTestPSTH(:, plotBins))/binWidth;
        case 2 % for 'late'
            odourTitle = 'Late';
            ctrl = sum(lateCtrlPSTH(:, divergenceBins));
            test = sum(lateTestPSTH(:, divergenceBins));
            plotCtrl = mean(lateCtrlPSTH(:, plotBins))/binWidth;
            plotTest = mean(lateTestPSTH(:, plotBins))/binWidth;
    end
    
    % plot the PSTHs:
    hAxData = subplot(2,2,1+(ee-1),'NextPlot','add');
    bar(hAxData,timeCentres(plotBins),plotCtrl,'FaceAlpha',0.75);
    bar(hAxData,timeCentres(plotBins),plotTest,'FaceAlpha',0.75);
    title([odourTitle,' odour']);
    xlabel('time (s)');
    ylabel('firing rate (/s)');
    legend('Control', 'Test');
    
    % perform the random walk calculations (see methods for details):
    cumdiff(:) = 0;
    normdiff(:) = 0;
    for kk = 1:nBinsDivergence
        for jj = kk:nBinsDivergence
            cumdiff(kk,jj) = hStabilisedDifference(test(kk:jj),ctrl(kk:jj));
            normdiff(kk,jj) = cumdiff(kk,jj)/sqrt(jj+1-kk);
        end
    end
    % normdiff is now an upper triangular matrix, with each row
    % representing a noise adjusted random walk starting at each time bin.
    % Obtain the maximum along each row to obtain the maximum divergence
    % for each time bin (symbol delta_k in the publication), and then the
    % time bin at which the largest of the random walks started:
    maxDiff = max(normdiff, [], 2);
    [~,ind] = max(maxDiff);
    
    % plot the maximum random walk distance versus start time of the random
    % walks:
    hAxSig = subplot(2,2,3+(ee-1));
    scatter(hAxSig, timeCentres(divergenceBins),maxDiff);
    hAxSig.XLim = hAxData.XLim;
    xlabel('accumulation start time (s)');
    ylabel('maximum divergence');
    % print the start time of the largest random walk to command window:
    fprintf('%s max deviation when start clock at bin %d = %.1f s\n', ...
        odourTitle,ind,(ind-1)*binWidth);
end

% to perform the stats showing that there is a difference between the
% response to odours in the corresponding reinforcement windows, but "not"
% (i.e. null hypothesis not rejected) in the non-reinforced windows:
normSigTh = 0.05; % for the normality test
odourLabels = {'early','late'}; % for labelling/ouput to command window
passfail = {'pass','fail'}; % for labelling/ouput to command window
fprintf('\n\n');
for rr = 1:2
    % rr loops over reinforcement windows
    for od = 1:2
        % od loops over odours
        % get the spike counts control/test odour presentation:
        ctrlData = hStatsCountsHelper(rr,(od-1)*2+1);
        testData = hStatsCountsHelper(rr,(od-1)*2+2);
        % perform t-test:
        [~,p_tTest,~,stats_tTest] = ttest(ctrlData,testData);
        % output to command window:
        % odour identity and reinforcement window:
        fprintf('%s reinforced, %d to %d seconds\n',odourLabels{od},reinforcementTimes(rr,:));
        % normality test output:
        fprintf('\tcontrol %s normality\n',passfail{jbtest(ctrlData,normSigTh)+1});
        fprintf('\ttest %s normality\n',passfail{jbtest(testData,normSigTh)+1});
        % t-test output:
        fprintf('\tt-test:\t\tp=%.4f\tdof=%d,\tt=%.4f\n',p_tTest,stats_tTest.df,stats_tTest.tstat);
        % wilcoxon output (this wasn't used in the publication because the
        % normality tests are all okay. But it shows that we still reach
        % the p<0.05 significance level without assuming normality.):
        fprintf('\tWilcoxon:\tp=%.4f\n\n',signrank(ctrlData,testData));
    end
end