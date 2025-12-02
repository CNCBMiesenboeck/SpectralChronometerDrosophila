%% To preprocess the behaviour data for Bayesian analysis and to produce the
% heatmaps in figure 5A and extended 3B.
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
% This script takes the the magnitude of the ball rotation during the odour
% presentation and (calls a function that) extracts the movement bouts and
% returns the start time of the most vigorous movement bout. These are the
% time points used for the Bayesian analysis. 
%
% This script generates the behaviour heatmaps and probability of movement
% initiation histograms shown in figure 5A and extended 3B.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To clear the workspace, screen, and load the data
clear
clc
loaded = load('data\behaviourResponses.mat');
% The loaded struct contains the responses to the testing odour pulse of 
% three experimental groups. Each data set is a struct with the following
% members:
% times             the (approx) time points of the camera frames relative 
%                   to the odour presentation. 
% movtMagnitude     the magnitude of the ball rotation in radians, as 
%                   described in the manuscript. The first array index
%                   iterates over the odours (paired with early aversive 
%                   stimulus, late punishment or not paired). The second
%                   index iterates over the flies and the final index
%                   corresponds to time points.
% sortOrder         The order of the flies once sorted as described in the
%                   manuscript and as we will reproduce here.
% sortedTimes       The start times of the maximum bouts as described in
%                   the manuscript and reproduced here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To extract the start times of the maximum movt bouts and sort the data.
% after running this section, "sortedStartTimes" should be identical to
% those used to produce the figures in the manuscript, which are also in
% the "sortedTimes" member of the corresponding data struct in loaded (e.g.
% sortedloaded.DANData.sortedTimes). These are the time points used for the
% Bayesian analysis.

% To select the data set:
dataSetInd = 1;         % 1=DAN, 2=chrimson control

threshold = exp(-4.5); % as determined by inspection of the histogram of
%                        the magnitude of the ball movement (with a log
%                        x-axis). This corresponds to approx = 1.75 mm/s,
%                        as stated in the manuscript.
minBoutLen = 8;        % 8 frames = ~160 ms
timeRange = [0,10];
epochNames = {'Early', 'Late', 'Neutral'};
dataSetOpts = {'DAN stim','Chr Ctrl'};
switch dataSetInd
    case 1
        times = loaded.DANData.times;
        movt = loaded.DANData.movtMagnitude;
    case 2
        times = loaded.chrCtrlData.times;
        movt = loaded.chrCtrlData.movtMagnitude;
end
frameRate = 1/(times(2)-times(1));
[nEpoch, nFlies, nTimepoints] = size(movt);
sortedBehaviour = zeros(size(movt));
[sortOrder, sortedStartTimes] = deal(zeros(nEpoch,nFlies));
for ee = 1:nEpoch
    % The actual processing is performed in the getMaxBoutStarts function
    % please read the comments (in the file "getMaxBoutStarts.m").
    % The function returns the sorted start index of the maximum movement
    % bout and the sorting order.
    [sortedTimeIndsTemp, sortOrder(ee,:)] = ...
        BehaviourFns.getMaxBoutStarts(...
        movt(ee,:,:), threshold, minBoutLen);
    % flies with no movement bouts above threshold have a start of max bout
    % at time index -1. We also exclude flies that have their most vigorous
    % movement bout starting before the odour onset (which is detected as
    % time index 1)
    moverInds = sortedTimeIndsTemp>0;
    nValidMoversTemp = sum(moverInds);
    % these are assigned a dummy value. It doesn't matter exactly what it
    % is. These flies aren't used in the histograms, or Bayesian analysis.
    sortedStartTimes(ee,~moverInds) = -abs(times(1));
    % convert indicies to time points
    sortedStartTimes(ee,moverInds) = ...
        times(sortedTimeIndsTemp(moverInds)+1);
    % to reproduce the behaviour heatmaps in the manuscript
    sortedBehaviour(ee,:,:) = movt(ee,sortOrder(ee,:),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To compute the probability of movement initiation histograms.
% The probability of movement initiation is a smoothed histogram. It is
% produced by convolution of the rasterised start times with a 1 second top
% hat function and normalised by the number of valid flies and area of the
% smoothing window to obtain expected counts per second
displaySubSample = 10; % this is just for display (matching the figure)
kernelWindow = 1; % 1 second window
dt = kernelWindow/displaySubSample;
pMovtTimeScale = timeRange(1):dt:timeRange(2);
nPMovtTime = numel(pMovtTimeScale);
[movtHistInit,pMovtInit] = deal(zeros(nEpoch,nPMovtTime));
kernelRelativeLowerInd = -ceil((displaySubSample-1)/2);
kernelRelativeUpperInd = floor((displaySubSample-1)/2);
hKernelCentreInd = @(tt) 1+round((tt-timeRange(1))/dt);
hKernelIndsHelper = @(tt,cc) max(1,cc+kernelRelativeLowerInd):min(nPMovtTime,cc+kernelRelativeUpperInd);
hKernelInds = @(tt) hKernelIndsHelper(tt,hKernelCentreInd(tt));
for ee = 1:nEpoch
    % exclude flies that are moving at the begining of the odour
    % presentation. NB times(1) is the dummy value encoding no movement
    validStartInds = sortedStartTimes(ee,:)>times(2);
    for tt = sortedStartTimes(ee,validStartInds)
        movtHistInit(ee,hKernelInds(tt))=movtHistInit(ee,hKernelInds(tt))+1;
    end
    % normalise to obtain the probability of movement initiation per
    % second. NB the normalisation is not by the number of flies because
    % the kernel is truncated when an initiation time is near the edge of
    % the time range.
    pMovtInit(ee,:) = movtHistInit(ee,:)/(dt*sum(movtHistInit(ee,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To display the sorted heatmaps and movement initiation histograms in
% figures 5A and extended 3B
figure('OuterPosition',[100,50,1200,1000]);
zScale = 1.5; % upper limit of the movement scale in ball rotations per second
plotOrder = [3,1,2];
for ee=1:nEpoch
    pp = plotOrder(ee);
    subplot(3,nEpoch,ee+[0,nEpoch]);
    imagesc(-squeeze(sortedBehaviour(pp,:,:)));
    colormap('gray');
    caxis([-zScale*2*pi/frameRate,0]);
    title(epochNames{pp});
    subplot(3,nEpoch,ee+2*nEpoch);
    plot(pMovtTimeScale,pMovtInit(pp,:),'LineWidth',3,'Color','k');
    set(gca,'YLim',[0,0.4]);
    title('probability of movement initiation');
end