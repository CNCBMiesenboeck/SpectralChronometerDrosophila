%% To reproduce the temporal classification shown in figure 2D and 3C.
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
% This script performs the temporal decoding of the e-phys recordings of
% Kenyon cells in response to odour pulses. The mathematical details are
% presented in the publication by Kropf, Talbot and Miesenboeck.
%
% This script will generate a multi-panel figure. The top row corresponds
% to the predictor matricies in publication figure 2d and 3c labelled 
% "Observed". The second row are the matricies labelled "Shuffled". The
% third row are not shown in the publication. These are the average over
% all of the permutations of the shuffling. A set of tables containing the
% coefficient of determination for the five predictor matricies for each
% individual odour presentation, and for the average of these predictor
% matricies (NB this is not the same as the average of the R^2 values for
% each matrix). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the data set:
dataSetID = 0;  % enter 0 or 1 for the 4-second, or the 10-second recordings

% clear the workspace, screen, and load the data
clearvars -except dataSetID
clc
if dataSetID==0
    load('data\KC_data.mat', 'stimComplete', 'on_significant', 'off_significant', ...
        'anatomicalNames');
    % please see the readme file for a description of the data.
    % Each variable is a 7-element cell array, one element for each cell type
    % in the order indicated in anatomicalNames. The remaining arrays have size
    % nCells-by-9-by-5. The first index corresponds to cell identity (so nCells
    % is the number of cells measured of that type). The second index is over
    % the presentations of clean air, mineral oil, followed by the seven 
    % odours. The final index is over the five repetitions. Further details:
    % stimComplete
    %   an array of the times of each spike in that recording
    % on_significant & off_significant
    %   if that response exceeded the significance threshold in the odour on
    %   phase, or odour off phase.
    
    % compute a logical array for each of the seven KC types indicating if each
    % cell odour combination is significant. The source arrays (on_significant
    % and off_significant) are both 1-by-7 cell arrays of logical arrays with
    % three dimensions indexing cell, odour, and repetition identities. The
    % logical value is identical for the five repetitions of any cell-odour
    % combination. The following line reorganises the source arrays into a
    % 1-by-7 cell array, each with one dimension indexing over cell-odour
    % combination (having removed the redundant repetition entries).
    significant = cellfun(@(x) squeeze(x(:,:,1)),...
        arrayfun(@(x) on_significant{x}|off_significant{x}, ...  
        1:7,'UniformOutput',false), 'UniformOutput',false);

    % Don't change these. These define the (combinations of) cell types as used 
    % in the publication, but in a different order. I maintain the original
    % order so that it is possible to reproduce the shuffling. publicationOrder
    % decodes the default cellTypes order into that shown in the publication.
    % To use different combinations of cell types, change the variable
    % "cellTypes" in the following section. 
    defaultCellTypes = {1,3,4,5,6,7,[1,3:7]};
    publicationOrder = [5,2,3,4,1,6,7];
else
    loaded = load('data\10s_KC_data.mat');
    stimComplete = loaded.odourOn;
    significant = loaded.on_significant;
    cellTypes = {1};
    publicationOrder = 1;
    anatomicalNames = {'a''/b'' middle (10 second responses)'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters:
if dataSetID==0
    % for the 4-second recordings
    
    % the time range of the classification. 
    timeRange = [0,4]; % [0,4] used in publication
    % the resolution of the classification. 
    temporalWindow = 0.2; % 0.2 in publication
    temporalStepSize = temporalWindow;% 0.2 in publication
    % the (combinations of) cell types to use. use defaultCellTypes to
    % reproduce the results in the publication. Alternatively, perhaps try
    % other anatomically appropriate combinations eg {[1,7],[3,6],[4,5]}.
    cellTypes = defaultCellTypes;
else
    % for the 10-second recordings
    timeRange = [0,10]; % [0,10] in publication 
    temporalWindow = 0.5; % 0.5 in publication
    temporalStepSize = 0.5; % 0.5 in publication
end
% the repetitions to use. You may want to see how well the temporal coding
% is preserved in eg repetitions 2:5. Obviously, with fewer repetitions,
% the noise will be increase.
repOpts = 1:5;

% whether to shuffle the data or not.
shuffleBetweenRepsFlag = true;
rndSeed = 1; % The random number generator seed. -1 for no seeding. 
nPermutations = 10000; % when shuffling, the number of permutations to average
% rndSeed=1 and nPermutations=10000 were used for the publication

% define some useful constants and arrays.
nReps = numel(repOpts); nRepComb = nReps + 1;
nLineCombs = numel(cellTypes);
temp = timeRange(1):temporalStepSize:(timeRange(2)-temporalWindow);
timeBins = [temp; temp+temporalWindow]';
nTimeBins = numel(temp);
timeCentres = temp+0.5*temporalWindow;

% this will hold the number of significant cell-odour responses
nResponsive = zeros(1,nLineCombs);
% this will hold the psth's for each cell type (combination).
[binnedResponses,logFactBinnedResp] = deal(cell(nLineCombs,1));
% the computation requires taking the log of the factorial of a large
% array of spike counts. But the spike counts will be integers up to about
% 20. Its much faster to make a lookup table. The lookup table will be
% expanded if necessary at the end of the following loop.
currentMaxBinCount = 20;
logFactorialLookup=log(factorial(0:currentMaxBinCount));
for cc = 1:nLineCombs
    % extract and print to screen the number of responses
    nResponsive(cc) = sum(cellfun(@(x) sum(x,'all'),significant(cellTypes{cc})));
    fprintf('%s %d pulse responses\n', ...
        sprintf('%s, ',anatomicalNames{cellTypes{cc}}),nResponsive(cc));
    
    % to temporarily hold the spike times of all significant cell-odour
    % responses before binning in the next loop.
    allResponsesTemp = cell(nResponsive(cc),nReps);
    tempCounter = 1;
    for ii = cellTypes{cc}
        stimCompleteTemp = stimComplete{ii};
        for jj = 1:size(significant{ii},1)
            for kk = 1:size(significant{ii},2)
                if significant{ii}(jj,kk)
                    for ll = 1:5
                        allResponsesTemp(tempCounter,ll) = ...
                            stimCompleteTemp{jj,kk,ll};
                    end
                    tempCounter = tempCounter+1;
                end
            end
        end
    end
    
    % to bin the spike times into psth's
    binnedResponses{cc} = zeros(nTimeBins,nResponsive(cc),nReps);
    for tt = 1:nTimeBins
        for rr = 1:nReps
            rInd = repOpts(rr);
            binnedResponses{cc}(tt,:,rr) = ...
                cellfun(@(x) sum(x>timeBins(tt,1)&x<=timeBins(tt,2)), ...
                allResponsesTemp(:,rInd));
        end
    end
    % check and update the log-factorial lookup table
    tempMax = max(binnedResponses{cc},[],'all');
    if tempMax>currentMaxBinCount
        currentMaxBinCount=tempMax;
        logFactorialLookup = log(factorial(0:currentMaxBinCount));
    end
    % the log-factorial can be calculated and stored once rather than at
    % each iteration of the shuffling.
    logFactBinnedResp{cc} = logFactorialLookup(binnedResponses{cc}+1);
end
if shuffleBetweenRepsFlag 
    if rndSeed>=0
        % seed the random number generator if necessary
        rng(rndSeed,'threefry');
    end
else
    nPermutations = 1;
end
permutationsFlag = nPermutations>1;

% pmfs are probability that temporal pattern of individual repetitions arise
% from the poisson distribution measured from the combination of the other 
% repetitions. 
% When shuffling, the pmfs at each permutation are computed independently as 
% currentPMFs and averaged as avePmfs. The currentPMFs overwritten at each 
% permutation. When performing the analysis on the observed data, currentPMFs is
% identical to avePMFs. 
% Index 1: cell type (or combination of types)
%       2: actual time bin
%       3: test time bin
%       4: repetition number. Plus last entry = average.
[currentPMFs,obsPMFs,shufPMFs] = deal(zeros(nLineCombs,nTimeBins,nTimeBins,nRepComb));
% To store the coefficient of determination of each permutation to build
% the distribution:
rSqDistr = zeros(nPermutations,nLineCombs,nRepComb);
onesTemp = ones(nTimeBins,1); % for convenience
% for indexing in a leave one out methodology:
rrTilde = cell2mat(arrayfun(@(x) sort(setxor((1:nReps),x)),1:nReps,'UniformOutput',false)');
% the calculation will be performed on a copy of the data. The following 
% is easier to write than preallocating:
extractedData = binnedResponses;
logFactData = logFactBinnedResp;
% used to reduce the number of times indexing into currentPMFs:
pmfTemp = zeros(nTimeBins,nTimeBins);
for bb = 0:nPermutations
    % First iteration of the loop (when bb==0), performs the calculation on the 
    % observed data. Then for bb>0, the calculation is performed on
    % nPermutations of shuffled data.
    if shuffleBetweenRepsFlag && bb>0
        % Shuffle if necessary
        for cc = 1:nLineCombs
            for rr = 1:nReps
                tempPerm = randperm(nResponsive(cc));
                extractedData{cc}(:,:,rr) = binnedResponses{cc}(...
                    :,tempPerm,rr);
                logFactData{cc}(:,:,rr) = logFactBinnedResp{cc}(:,tempPerm,rr);
            end
        end
        if mod(bb,100)==0
            % print a progress update every 100 iterations
            fprintf('shuffling progress: %d/%d complete\r',bb,nPermutations);
        end
    end
    for cc = 1:nLineCombs
        for rr = 1:nReps
            % loop over the five odour pulses
            
            % compute the average spike rate of the other four repetitions:
            aveTrainRate = (sum(extractedData{cc}(:,:,rrTilde(rr,:)),3))'/...
                (nReps-1);
            % to deal with zeros (see the publication for details):
            aveZeroInds = aveTrainRate==0;
            aveTrainRate(aveZeroInds) = -log(0.5)/(nReps-1);
            logAveTrainRate = log(aveTrainRate);
            sumAveTrainRate = sum(aveTrainRate,1);
            pmfTemp(:)=0;
            for tt = 1:nTimeBins
                % compute the Poisson log-probability that the data along
                % the time axis of the current odour presentation (rr)
                % originates from that of the other presentations rrTilde.
                % The matrix multiplication sums the log-probabilities
                % across cells (ie "AND"s the probabilities).
                indLogP = extractedData{cc}(tt,:,rr)*logAveTrainRate - ...
                    sum(logFactData{cc}(tt,:,rr)) - sumAveTrainRate;
                % "de-log":
                nonNormPMF = exp(indLogP);
                % normalise:
                pmfTemp(tt,:) =  nonNormPMF/sum(nonNormPMF,'all');
            end
            currentPMFs(cc,:,:,rr) = pmfTemp;
        end
    end
    % compute the average over the five odour presentations for the current
    % permutation
    currentPMFs(:,:,:,nRepComb) = mean(currentPMFs(:,:,:,1:nReps),4);
    if permutationsFlag && bb>0
        rSqDistr(bb,:,:) = coefficientOfDetermination(currentPMFs,timeCentres);
    end
    if bb==0
        obsPMFs=currentPMFs;
    else
        shufPMFs=shufPMFs+currentPMFs;
    end
end
shufPMFs=shufPMFs/nPermutations;

% compute the average over the five odour presentations for the observed and 
% then average of the shuffled data:
obsPMFs(:,:,:,nRepComb) = mean(obsPMFs(:,:,:,1:nReps),4);
shufPMFs(:, :, :, nRepComb) = mean(shufPMFs(:, :, :, 1:nReps), 4);
% and computed the coefficients of determination:
[rSqObs,rSqObsMean] = coefficientOfDetermination(obsPMFs,timeCentres);
[rSqShufAve, rSqShufAveRepMean] = coefficientOfDetermination(shufPMFs,timeCentres);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To produce the prediction matrix images. To reproduce the figs, set (below)
% repInd=nRepComb and displayOrder=publicationOrder (assuming default value
% for cellTypes). 

% which repetition index to display (set to nRepComb for the average
% prediction matrix):
repInd = nRepComb;

if nLineCombs==numel(publicationOrder)
    displayOrder = publicationOrder;
else
    displayOrder = 1:nLineCombs;
end
nCTI = numel(displayOrder);
nRowsTemp = 1+shuffleBetweenRepsFlag+2*(permutationsFlag);
figure('OuterPosition',[50,30,225*nCTI,250*nRowsTemp]);
[predMatObs, predMatShufSingle, predMatShufAve] = deal(zeros(nCTI,nTimeBins,nTimeBins));
histBinWidth = 0.1;
rSqHistEdges = [-inf,-3:histBinWidth:1];
rSqHistCounts = zeros(nCTI,numel(rSqHistEdges)-1);
rSqObsRep = squeeze(rSqObs(displayOrder,repInd))';
rSqShufAveRep = squeeze(rSqShufAve(displayOrder,repInd))';
rSqShufExample = squeeze(rSqDistr(end,displayOrder,repInd));
rSqDistribution = squeeze(rSqDistr(:,displayOrder,repInd));
hSubplot = @(ii) subplot(nRowsTemp,nCTI,ii);
hImShow = @(im) imshow(1-imrotate(squeeze(im),90),'InitialMagnification',700);
hRSqTitle = @(rSq) title(sprintf('R^2=%.2f',rSq));
for cc = 1:nCTI
    predMatObs(cc,:,:) = squeeze(obsPMFs(displayOrder(cc),:,:,repInd));
    hSubplot(cc);
    hImShow(predMatObs(cc,:,:));
    hRSqTitle(rSqObsRep(cc));
    if shuffleBetweenRepsFlag
        predMatShufSingle(cc,:,:) = squeeze(currentPMFs(displayOrder(cc),:,:,repInd));
        hSubplot(cc+nCTI);
        hImShow(predMatShufSingle(cc,:,:));
        hRSqTitle(rSqShufExample(cc));
        if permutationsFlag
            predMatShufAve(cc,:,:) = squeeze(shufPMFs(displayOrder(cc),:,:,repInd));
            hSubplot(cc+2*nCTI);
            hImShow(predMatShufAve(cc,:,:));
            hRSqTitle(rSqShufAveRep(cc));
            hSubplot(cc+3*nCTI);
            rSqHistCounts(cc,:) = histcounts(rSqDistribution(:,displayOrder(cc)),rSqHistEdges);
            bar(rSqHistEdges(2:(end-1))+histBinWidth/2,rSqHistCounts(cc,2:end),1);
            set(gca,'XLim',[-2,1]);
            set(gca,'YLim',[0,2500]);
            title(sprintf('R^2=%.2f +/- %.2f', ...
                mean(rSqDistribution(:,displayOrder(cc))), ...
                std(rSqDistribution(:,displayOrder(cc)))));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To tabulate the coefficient of determination
if nLineCombs==numel(publicationOrder)
    displayOrder = publicationOrder;
else
    displayOrder = 1:nLineCombs;
end

hRsqTableHelper = @(str,rSqArr,offset) uitable(figure(...
    'Name',['Coefficient of determination ',str],'OuterPosition',[350,50+offset,1050,250]),...
    'Data',rSqArr(displayOrder,:)','ColumnName',cellfun(@(x) x(1:(end-1)),...
    cellfun(@(x) sprintf('%s,',anatomicalNames{x}),cellTypes(displayOrder),...
    'UniformOutput',false),'UniformOutput',false),...
    'RowName',[cellstr("repetition "+string(repOpts)),'R^2 of mean pmf'], ...
    'Position',[20,20,1000,132]);
hRsqTableHelper('observed',rSqObs,500);
if shuffleBetweenRepsFlag
    hRsqTableHelper('shuffle average',rSqShufAve,0);
    if permutationsFlag
        hRsqTableHelper('shuffle example',reshape(rSqDistr(end,:,:),nLineCombs,nRepComb),250);
    end
end

%%
function [rSqOut, rSqRepMean] = coefficientOfDetermination(pmfsIn,timeCentres)
% compute the coefficient of determination for prediction matricies for
% each repetition and cell-type combination (ie R^2 over the two inner
% indicies of pmfs). The comparison is to perfect prediction. The common
% terminology is:
%   R^2 = 1 - SS_res/SS_tot
%     SS_res = sum of squared residuals i.e. (pmfs-perfect prediction)^2
%     SS_tot = variance of marginal pmfs along predicted time axis

nLineCombs = size(pmfsIn,1);
nRepComb = size(pmfsIn,4);
nReps = nRepComb-1;

% Compute the matrix representing the deviation from perfect prediction
predTimes2D = timeCentres - timeCentres';
% Replicate and permute into the 4-D format of pmfs:
predTSq4D = permute(repmat(predTimes2D.*predTimes2D, ...
    1,1,nLineCombs,nRepComb),[3,1,2,4]);
% compute the deviation from perfect prediction (= SS_res)
SS_res = mean(sum(pmfsIn.*predTSq4D,3),2);

% compute the marginal pmfs along predicted time axis
marginalPmfs = mean(pmfsIn,2);
% replicate and permute the time axis into dimensionality of marginalPmfs
times4D = permute(repmat(timeCentres,1,1,nLineCombs,nRepComb),...
    [3,1,2,4]);
% compute SS_tot using var(x) = E(x^2) - E(x)^2:
mean_predTimeSquared = sum(marginalPmfs.*times4D.*times4D,3);
mean_predictedTime = sum(marginalPmfs.*times4D,3);
SS_tot = mean_predTimeSquared - mean_predictedTime.*mean_predictedTime;
% compute R^2. First index is over cell-type combination. Second index is
% over the repetition, with the final entry being for the average pmf.
rSqOut = reshape(1-SS_res./SS_tot,nLineCombs,nRepComb);
rSqRepMean = mean(rSqOut(:,1:nReps),2);
end