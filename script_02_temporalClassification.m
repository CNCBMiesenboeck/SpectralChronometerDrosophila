%% To reproduce the temporal classification and resolution analysis
% as shown in figures 2D, 2E, 3C, 3D, 4C, and S2C
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
% The script will generate (versions of) the plots shown in figures 2D, 2E,
% 3C, 3D, 4C, and S2C. Alternatively, a customised analysis can be
% performed and plots equivalent (or similar) to those in the figures can
% be produced.
%
% Execute section 1 to perform the analysis, changing the options as
% indicated to select the data and set the parameters. Figures will be
% automatically generated for preset options. Execute sections 2 to 5,
% modifying options, generate alternative figures.

%% Section 1
% There are various settings in the script. Either choose from a preset 
% corresponding to one of the figures, or set them in the customOptions 
% struct below.
optionsIndex = 1; % 1-9 = figures (see figureID_Opts below), 10 = custom
figureID_Opts = {...
    '2D', ...           1
    '2E', ...           2
    '3C', ...           3
    '3D', ...           4
    '4C top', ...       5
    '4C middle', ...    6
    '4C bottom', ...    7
    'S2C grey', ...     8
    'S2C black', ...    9
    'custom'};        % 10 to use customOptions. See below

% set to 0 or 1 for the 4-second or 10-second recordings:
customOptions.dataSetID = 0;
% Include responses that pass significance test (true) or to include all
% responses (false). Figures equivalent to 3C, and 4D that include all 
% responses can be produced by setting this to false:
customOptions.significantFlag = true;
% array of odours to include (1=sulcatone, 2=heptanone, 3=ethyl acetate, 
% 4=benzaldehyde, 5=pentanol, 6=octanol,7=isopentyl acetate):
customOptions.odourSelection = 1:7;
% time window to analyse:
customOptions.timeRange = [0,4];
% the binning window, and step sizes (to loop over). In the publication, 
% the step sizes were identical to the window sizes, but there there are
% possible combinations to improve the temporal resolution:
customOptions.temporalWindowOpts = 0.2;
customOptions.temporalWindowStepSizeOpts = 0.2;
% the cell type (combinations) to analyse. Alternatively, to the figure 
% options, perhaps try other anatomically appropriate combinations eg 
% {[1,7],[3,6],[4,5]}. For the 10-s recordings, this is limited to {1}.
% 1=g_main, 2=a/b_posterior, 3=a/b_surface, 4=a'/b'_middle, 5=a'p'_ap,
% 6=a/b_core, 7=g_dorsal:
customOptions.cellTypes = {1,3,4,5,6,7,[1,3:7]};
% the repetitions to use. You may want to see how well the temporal coding
% is preserved in eg repetitions 2:5. Obviously, with fewer repetitions,
% the noise will be increase:
customOptions.repOpts = 1:5;
% whether to perform the shuffling or not, which seed to use (-1 for no
% seeding), and the number of permutations to use:
customOptions.shuffleBetweenRepsFlag = true;
customOptions.rndSeed = 1;
customOptions.nPermutations = 10000;

switch optionsIndex
    case num2cell(1:9)
        optionsStruct = getPresetOptions(figureID_Opts{optionsIndex});
    case 10
        optionsStruct = customOptions;
    otherwise
        disp('Invalid optionsIndex. See line 31.')
        return;
end

% clear the workspace, screen, and load the data
clearvars -except optionsStruct optionsIndex figureID_Opts
clc
if optionsStruct.dataSetID==0
    load('data\KC_data.mat', 'stimComplete', 'on_significant', 'off_significant', ...
        'anatomicalNames');
    % please see the readme file for a description of the data.
    % Each variable is a 7-element cell array, one element for each cell
    % type in the order indicated in anatomicalNames. The remaining arrays 
    % have size nCells-by-9-by-5. The first index corresponds to cell 
    % identity (so nCells is the number of cells measured of that type). 
    % The second index is over the presentations of clean air, mineral oil,
    % followed by the seven odours. The final index is over the five
    % repetitions. Further details:
    % stimComplete
    %   an array of the times of each spike in that recording
    % on_significant & off_significant
    %   if that response exceeded the significance threshold in the odour on
    %   phase, or odour off phase.
    
    % compute a logical array for each of the seven KC types indicating if 
    % each cell odour combination is significant. The source arrays 
    % (on_significant and off_significant) are both 1-by-7 cell arrays of 
    % logical arrays with three dimensions indexing cell, odour, and 
    % repetition identities. The logical value is identical for the five 
    % repetitions of any cell-odour combination. The following line
    % reorganises the source arrays into a 1-by-7 cell array, each with one
    % dimension indexing over cell-odour combination (having removed the 
    % redundant repetition entries).
    significant = cellfun(@(x) squeeze(x(:,:,1)),...
        arrayfun(@(x) on_significant{x}|off_significant{x}, ...  
        1:7,'UniformOutput',false), 'UniformOutput',false);

    % Don't change these. These define the (combinations of) cell types as 
    % used in fig 2 of the publication, but in a different order. I have 
    % maintained the original order so that it is possible to reproduce the
    % shuffling. publicationOrder decodes the default cellTypes order into 
    % that shown in fig 2. To use different combinations of cell types, 
    % change the "cellTypes" field of the, and use the, customOptions
    % struct near the top of this script.
    defaultCellTypes = {1,3,4,5,6,7,[1,3:7]};
    publicationOrder = [5,2,3,4,1,6,7];
else
    loaded = load('data\10s_KC_data.mat');
    stimComplete = loaded.odourOn;
    significant = loaded.on_significant;
    optionsStruct.cellTypes = {1};
    defaultCellTypes = {1};
    publicationOrder = 1;
    anatomicalNames = {'a''/b'' middle (10 second responses)'};
end

odourNames = {...
    'sulcatone','2-heptanone','ethyl acetate','benzaldehyde',... 1, 2, 3, 4
    '2-pentanol','1-octanol','isopentyl acetate'};             % 5, 6, 7
if numel(optionsStruct.odourSelection) == 7
    odourString = 'all odours';
else
    odourString = sprintf('%s, ',odourNames{optionsStruct.odourSelection});
    odourString = ['odours: ', odourString(1:(end-2))];
end
fprintf('Analysing %s\n', odourString);

if ~optionsStruct.shuffleBetweenRepsFlag
    optionsStruct.nPermutations = 0;
end

% The analysis will be performed on the recordings indicated in the logical
% cell array selectedRecordings. In the first case (all significant odours,
% corresponding to e.g. figure 2D), it is a copy "significant"
selectedRecordings = significant;
if ~optionsStruct.significantFlag
    for cc = 1:numel(selectedRecordings)
        selectedRecordings{cc}(:,3:9) = true;
    end
end
if ~isempty(optionsStruct.odourSelection)
    odourSelectionBool = false(1,9);
    odourSelectionBool(optionsStruct.odourSelection+2) = true;
    for cc = 1:numel(selectedRecordings)
        selectedRecordings{cc}(:,:) = selectedRecordings{cc}(:,:) & ...
            repmat(odourSelectionBool,size(selectedRecordings{cc},1),1);
    end
end

% define some useful constants and arrays.
nReps = numel(optionsStruct.repOpts); nRepComb = nReps + 1;
nLineCombs = numel(optionsStruct.cellTypes);

nWindowOpts = numel(optionsStruct.temporalWindowOpts);
[currentPMFs,obsPMFs,shufPMFs,RMS_t_Obs,RMS_t_distr] = ...
    deal(cell(1,nWindowOpts));
% To store the coefficient of determination of each permutation to build
% the distribution:
rSqDistr = zeros(nWindowOpts,optionsStruct.nPermutations,nLineCombs,nRepComb);
[rSqObs,rSqShufAve,RMS_Obs] = deal(zeros(nWindowOpts,nLineCombs,nRepComb));

for ww = 1:nWindowOpts
    % loop through the temporal binning window widths
    temporalWindow = optionsStruct.temporalWindowOpts(ww);
    temporalStepSize = optionsStruct.temporalWindowStepSizeOpts(ww);
    
    temp = optionsStruct.timeRange(1):temporalStepSize:...
        (optionsStruct.timeRange(2)-temporalWindow);
    timeBins = [temp; temp+temporalWindow]';
    nTimeBins = numel(temp);
    timeCentres = temp+0.5*temporalWindow;

    % this will hold the number of selected cell-odour responses
    nResponsive = zeros(1,nLineCombs);
    % this will hold the psth's for each cell type (combination).
    [binnedResponses,logFactBinnedResp] = deal(cell(nLineCombs,1));
    % the computation requires taking the log of the factorial of a large
    % array of spike counts. But the spike counts will be integers up to
    % about 20. Its much faster to make a lookup table. The lookup table 
    % will be expanded if necessary at the end of the following loop.
    currentMaxBinCount = 20;
    logFactorialLookup=log(factorial(0:currentMaxBinCount));
    for cc = 1:nLineCombs
        % extract and print to screen the number of responses
        nResponsive(cc) = sum(cellfun(@(x) sum(x,'all'),...
            selectedRecordings(optionsStruct.cellTypes{cc})));
        if ww==1
            fprintf('%s %d pulse responses\n', ...
                sprintf('%s, ',...
                anatomicalNames{optionsStruct.cellTypes{cc}}), ...
                nResponsive(cc));
        end

        % to temporarily hold the spike times of all selected cell-odour
        % responses before binning in the next loop.
        allResponsesTemp = cell(nResponsive(cc),nReps);
        tempCounter = 1;
        for ii = optionsStruct.cellTypes{cc}
            stimCompleteTemp = stimComplete{ii};
            for jj = 1:size(selectedRecordings{ii},1)
                for kk = 1:size(selectedRecordings{ii},2)
                    if selectedRecordings{ii}(jj,kk)
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
                rInd = optionsStruct.repOpts(rr);
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
        % the log-factorial can be calculated and stored once rather than
        % at each iteration of the shuffling.
        logFactBinnedResp{cc} = logFactorialLookup(binnedResponses{cc}+1);
    end
    if optionsStruct.shuffleBetweenRepsFlag 
        if optionsStruct.rndSeed>=0
            % seed the random number generator if necessary
            rng(optionsStruct.rndSeed,'threefry');
        end
    end
    permutationsFlag = optionsStruct.nPermutations>0;

    % pmfs are probability that temporal pattern of individual repetitions 
    % arisefrom the poisson distribution measured from the combination of
    % the other repetitions. 
    % When shuffling, the pmfs at each permutation are computed 
    % independently as currentPMFs and averaged as avePmfs. The currentPMFs
    % overwritten at each permutation. When performing the analysis on the
    % observed data, currentPMFs is identical to avePMFs. 
    % Index 1: cell type (or combination of types)
    %       2: actual time bin
    %       3: test time bin
    %       4: repetition number. Plus last entry = average.
    [currentPMFs{ww},obsPMFs{ww},shufPMFs{ww}] = deal(...
        zeros(nLineCombs,nTimeBins,nTimeBins,nRepComb));
    RMS_t_Obs{ww} = zeros(nTimeBins,nLineCombs,nRepComb);

    onesTemp = ones(nTimeBins,1); % for convenience
    % for indexing in a leave one out methodology:
    rrTilde = cell2mat(arrayfun(@(x) sort(setxor((1:nReps),x)),1:nReps,...
        'UniformOutput',false)');
    % the calculation will be performed on a copy of the data. The
    % following is easier to write than preallocating:
    extractedData = binnedResponses;
    logFactData = logFactBinnedResp;
    % used to reduce the number of times indexing into currentPMFs:
    pmfTemp = zeros(nTimeBins,nTimeBins);
    for bb = 0:optionsStruct.nPermutations
        % First iteration of the loop (when bb==0), performs the
        % calculation on the observed data. Then for bb>0, the calculation
        % is performed on nPermutations of shuffled data.
        if optionsStruct.shuffleBetweenRepsFlag && bb>0
            % Shuffle if necessary
            for cc = 1:nLineCombs
                for rr = 1:nReps
                    tempPerm = randperm(nResponsive(cc));
                    extractedData{cc}(:,:,rr) = binnedResponses{cc}(...
                        :,tempPerm,rr);
                    logFactData{cc}(:,:,rr) = logFactBinnedResp{cc}(...
                        :,tempPerm,rr);
                end
            end
            if mod(bb,100)==0
                % print a progress update every 100 iterations
                if nWindowOpts > 1
                    fprintf('temporal binning step: %d/%d (%.2fs)\t',...
                        ww,nWindowOpts,temporalWindow);
                end
                fprintf('shuffling progress: %d/%d complete\r',...
                    bb,optionsStruct.nPermutations);
            end
        elseif ~optionsStruct.shuffleBetweenRepsFlag
            if nWindowOpts > 1
                fprintf('temporal binning step: %d/%d (%.2fs)\r',...
                    ww,nWindowOpts,temporalWindow);
            end        
        end
        for cc = 1:nLineCombs
            for rr = 1:nReps
                % loop over the five odour pulses

                % compute the average spike rate of the other four 
                % repetitions:
                aveTrainRate = ...
                    sum(extractedData{cc}(:,:,rrTilde(rr,:)),3)'/(nReps-1);
                % to deal with zeros (see the publication for details):
                aveZeroInds = aveTrainRate==0;
                aveTrainRate(aveZeroInds) = -log(0.5)/(nReps-1);
                logAveTrainRate = log(aveTrainRate);
                sumAveTrainRate = sum(aveTrainRate,1);
                pmfTemp(:)=0;
                for tt = 1:nTimeBins
                    % compute the Poisson log-probability that the data
                    % along the time axis of the current odour presentation
                    % (rr) originates from that of the other presentations
                    % rrTilde. The matrix multiplication sums the
                    % log-probabilities across cells (ie "AND"s the
                    % probabilities).
                    indLogP = extractedData{cc}(tt,:,rr)*logAveTrainRate - ...
                        sum(logFactData{cc}(tt,:,rr)) - sumAveTrainRate;
                    if true
                        indLogP = indLogP - max(indLogP);
                    end
                    % "de-log":
                    nonNormPMF = exp(indLogP);
                    % normalise:
                    pmfTemp(tt,:) = nonNormPMF/sum(nonNormPMF,'all');
                end
                currentPMFs{ww}(cc,:,:,rr) = pmfTemp;
            end
        end
        % compute the average over the five odour presentations for the
        % current permutation
        currentPMFs{ww}(:,:,:,nRepComb) = ...
            mean(currentPMFs{ww}(:,:,:,1:nReps),4);
        if permutationsFlag && bb>0
            rSqDistr(ww,bb,:,:) = ...
                coefficientOfDetermination(...
                currentPMFs{ww},timeCentres,temporalWindow);
        end
        if bb==0
            obsPMFs{ww}=currentPMFs{ww};
        else
            shufPMFs{ww}=shufPMFs{ww}+currentPMFs{ww};
        end
    end
    shufPMFs{ww}=shufPMFs{ww}/optionsStruct.nPermutations;

    % compute the average over the five odour presentations for the 
    % observed and then average of the shuffled data:
    obsPMFs{ww}(:,:,:,nRepComb) = mean(obsPMFs{ww}(:,:,:,1:nReps),4);
    shufPMFs{ww}(:,:,:,nRepComb) = mean(shufPMFs{ww}(:,:,:,1:nReps),4);
    % and computed the coefficients of determination:
    [rSqObs(ww,:,:),RMS_Obs(ww,:,:),RMS_t_Obs{ww}(:,:,:)] = ...
        coefficientOfDetermination(obsPMFs{ww},timeCentres,temporalWindow);
    rSqShufAve(ww,:,:) = ...
        coefficientOfDetermination(shufPMFs{ww},timeCentres,temporalWindow);
end

if any(optionsIndex==[1:3,5:7])
%% Section 2
% To produce the prediction matrix images and null coefficient of
% determination histograms.
%
% Plots are shown in a multi-panel figure. The top row corresponds to the
% matricies in the publication labelled "Observed". The second row are 
% those labelled "Shuffled". The third row are not shown in the 
% publication. These are the average over all of the permutations of the
% shuffling. The fourth row shows the null distributions of the coefficient
% of determination. A set of tables showing R^2 for the five predictor
% matricies of each individual odour presentation, and for the average of
% these predictor matricies (NB this is not the same as the average of the
% R^2 values for each matrix). 

% The plots can be produced from an individual odour presentation, or the
% average prediction matrix. Set repInd to nRepComb for the average
% prediction matrix, or select one of the analysed presentations as stored 
% optionsStruct.repOpts (1 to 5 by default). The average prediction
% matricies are shown in the publication.
repInd = nRepComb;

% Select the temporal binning window option if necessary. Usually leave as
% 1.
windowOpt = 1;

% For figure 2, this will set the order of the matricies (by cell type)
% to that shown in the publication. Please see comments preceeding the
% definition of defaultCellTypes and publicationOrder (around line 121).
if nLineCombs==numel(publicationOrder) && all(arrayfun(...
        @(ii) all(optionsStruct.cellTypes{ii}==defaultCellTypes{ii}), ...
        1:nLineCombs))
    displayOrder = publicationOrder;
else
    displayOrder = 1:nLineCombs;
end
nCTI = numel(displayOrder); % nCTI = number of Cell Type Indicies
nRowsTemp = 1+optionsStruct.shuffleBetweenRepsFlag+2*(permutationsFlag);
hFig = figure('OuterPosition',[50,30,225*nCTI,250*nRowsTemp]);

% preallocate the decoder matricies
[predMatObs, predMatShufExample, predMatShufAve] = ...
    deal(zeros(nCTI,nTimeBins,nTimeBins));

% set up histogram variables
histBinWidth = 0.1;                 % default: 0.1 
rSqHistEdges = [-inf,-2:histBinWidth:1]; %     [-inf,-2:histBinWidth:1]
rSqHistCounts = zeros(numel(rSqHistEdges)-1,nCTI);

% extract the R^2 values for the relevant temporal binning and odour
% presentation options:
rSqObsRep = squeeze(rSqObs(windowOpt,displayOrder,repInd));
if optionsStruct.shuffleBetweenRepsFlag
    rSqShufAveRep = ...
        reshape(rSqShufAve(windowOpt,displayOrder,repInd),nLineCombs,1);
    rSqShufExample = ...
        reshape(rSqDistr(windowOpt,end,displayOrder,repInd),nLineCombs,1);
    rSqDistribution = ...
        reshape(rSqDistr(windowOpt,:,displayOrder,repInd),...
        optionsStruct.nPermutations,nLineCombs);
end

% set up some helper functions to display the matricies
hImShow = @(im) imshow(1-imrotate(squeeze(im),90),'InitialMagnification',700);
hNewAxes = @(ii) subplot(nRowsTemp,nCTI,ii);
hrSqTitle = @(rSq) title(sprintf('R^2=%.2f',rSq));
hDispFig = @(mat,ind,shufFlag,rSq) ...
    disp_predMat(hImShow,hNewAxes,hrSqTitle,rSq,mat,ind,ind+nCTI*shufFlag);

for cc = 1:nCTI
    predMatObs(cc,:,:) = ...
        squeeze(obsPMFs{windowOpt}(displayOrder(cc),:,:,repInd));
    hAx = hDispFig(predMatObs,cc,0,rSqObsRep(cc));
    if nCTI==1
        hAx.Title.String = sprintf('%s\n%s',odourString,hAx.Title.String);
    end
    if cc==1
        hAx.YLabel.String = 'Observed';
    end
    if optionsStruct.shuffleBetweenRepsFlag
        predMatShufExample(cc,:,:) = ...
            squeeze(currentPMFs{windowOpt}(displayOrder(cc),:,:,repInd));
        hAx = hDispFig(predMatShufExample,cc,1,rSqShufExample(cc));
        if cc==1
            hAx.YLabel.String = 'Shuffle example';
        end
        if  permutationsFlag
            predMatShufAve(cc,:,:) = ...
                squeeze(shufPMFs{windowOpt}(displayOrder(cc),:,:,repInd));
            hAx = hDispFig(predMatShufAve,cc,2,rSqShufAveRep(cc));
            if cc==1
                hAx.YLabel.String = 'Shuffle mean (not in figure)';
            end                
            hNewAxes(cc+3*nCTI);
            rSqHistCounts(:,cc) = histcounts(...
                rSqDistribution(:,displayOrder(cc)),rSqHistEdges);
            bar(rSqHistEdges(2:(end))+histBinWidth/2,...
                rSqHistCounts(1:end,cc),1);
            yLims = [0,2500];
            set(gca,'XLim',[rSqHistEdges(2),1]);
            set(gca,'YLim',yLims);

            title(sprintf('R^2=%.2f +/- %.2f', ...
                mean(rSqDistribution(:,displayOrder(cc))), ...
                std(rSqDistribution(:,displayOrder(cc)))));
        end
    end
end

%% Section 3
% To tabulate the coefficient of determination

% Select the temporal binning window option.
windowOpt = 1;

if nLineCombs==numel(publicationOrder)
    displayOrder = publicationOrder;
else
    displayOrder = 1:nLineCombs;
end

hRsqTableHelper = @(str,rSqArr,offset) uitable(figure(...
    'Name',['Coefficient of determination ',str],...
    'OuterPosition',[350,50+offset,1000,250]),...
    'Data',reshape(rSqArr(windowOpt,displayOrder,:),nLineCombs,nRepComb)',...
    'ColumnName',cellfun(@(x) x(1:(end-1)),...
    cellfun(@(x) sprintf('%s,',anatomicalNames{x}),...
    optionsStruct.cellTypes(displayOrder),...
    'UniformOutput',false),'UniformOutput',false),...
    'RowName',[cellstr("repetition "+string(optionsStruct.repOpts)),...
    'R^2 of mean pmf'], ...
    'Position',[20,20,1000,132]);
hRsqTableHelper('observed',rSqObs,500);
if optionsStruct.shuffleBetweenRepsFlag
    hRsqTableHelper('shuffle average',rSqShufAve,0);
    if permutationsFlag
        hRsqTableHelper('shuffle example',reshape(rSqDistr(:,end,:,:),...
            nWindowOpts,nLineCombs,nRepComb),250);
    end
end
end

if optionsIndex == 4
%% Section 4
% To plot the RMS error in prediction time versus time point (figure 3D)

% select the odour presentation (repInd), cell type combination (lineInd)
% and temporal binning width (windowOpt)
repInd = nRepComb;
lineInd = nLineCombs;
windowOpt = 1;
RMS_extracted = RMS_t_Obs{windowOpt}(:,lineInd,repInd);
if optionsStruct.dataSetID==0
    RMS_timeAx = [timeCentres(timeCentres<=2),timeCentres(timeCentres>2)-2];
else
    RMS_timeAx = timeCentres(1:(end-1));
    RMS_extracted(end) = [];
end
hAx = axes(figure());
scatter(RMS_timeAx,RMS_extracted,36,'filled');
hAx.XLabel.String = 'Elapsed time (since odour on/off edge) (s)';
hAx.YLabel.String = 'RMS error in prediction time (s)';
hAx.FontSize = 14;
hAx.YLim = [0,2.5];
annotation('textbox',...
    'String',sprintf('R=%.2f',PCC_helper(RMS_timeAx,RMS_extracted)),...
    'Position',[0.3,0.7,0.1,0.1], ...
    'FitBoxToText','on', ...
    'FontSize', 14);
end

if any(optionsIndex==[8,9])
%% Section 5
% To plot the average RMS error in prediction time versus time window width
% (figure S2C)
repInd = nRepComb;
lineInd = nLineCombs;
RMS_extracted = squeeze(RMS_Obs(:,lineInd,repInd));
hAx = axes(figure());
scatter(optionsStruct.temporalWindowOpts,RMS_extracted,36,'filled');
hAx.XLabel.String = 'Time bin width (s)';
hAx.YLabel.String = 'Mean RMS error in prediction time (s)';
hAx.FontSize = 14;
hAx.YLim = [0,1.5];
end

function [rSqOut, RMS, RMS_t] = ...
    coefficientOfDetermination(pmfsIn,timeCentres,temporalWindow)
% compute the coefficient of determination for prediction matricies for
% each repetition and cell-type combination (ie R^2 over the two inner
% indicies of pmfs). The comparison is to perfect prediction. The common
% terminology is:
%   R^2 = 1 - SS_res/SS_tot
%     SS_res = sum of squared residuals i.e. (pmfs-perfect prediction)^2
%     SS_tot = variance of marginal pmfs along predicted time axis

nLineCombs = size(pmfsIn,1);
nRepComb = size(pmfsIn,4);

% Compute the matrix representing the deviation from perfect prediction
predTimes2D = timeCentres - timeCentres';
% Replicate and permute into the 4-D format of pmfs:
predTSq4D = permute(repmat(predTimes2D.*predTimes2D, ...
    1,1,nLineCombs,nRepComb),[3,1,2,4]);
% compute the deviation from perfect prediction (= SS_res)
temp = sum(pmfsIn.*predTSq4D,3);
SS_res = mean(temp,2);
RMS = squeeze(sqrt(SS_res + temporalWindow*temporalWindow/12));
RMS_t = permute(sqrt(temp + temporalWindow*temporalWindow/12),[2,1,4,3]);
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
end

function hAx = disp_predMat(hImShow,hNewAxes,hTitle,val,mat,ind,axInd)
hAx = hNewAxes(axInd);
hImShow(mat(ind,:,:));
hTitle(val);
end

function optsOut = getPresetOptions(figureID)
switch figureID
    case '2D'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.2;
        optsOut.temporalWindowStepSizeOpts = 0.2;
        optsOut.cellTypes = {1,3,4,5,6,7,[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '2E'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = false;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.2;
        optsOut.temporalWindowStepSizeOpts = 0.2;
        optsOut.cellTypes = {1,3,4,5,6,7,[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '3C'
        optsOut.dataSetID = 1;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,10];
        optsOut.temporalWindowOpts = 0.5;
        optsOut.temporalWindowStepSizeOpts = 0.5;
        optsOut.cellTypes = {1};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '3D'
        optsOut.dataSetID = 1;
        optsOut.significantFlag = false;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,10];
        optsOut.temporalWindowOpts = 0.5;
        optsOut.temporalWindowStepSizeOpts = 0.5;
        optsOut.cellTypes = {1};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '4C top'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 5;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.2;
        optsOut.temporalWindowStepSizeOpts = 0.2;
        optsOut.cellTypes = {[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '4C middle'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 7;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.2;
        optsOut.temporalWindowStepSizeOpts = 0.2;
        optsOut.cellTypes = {[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case '4C bottom'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 1;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.2;
        optsOut.temporalWindowStepSizeOpts = 0.2;
        optsOut.cellTypes = {[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = true;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case 'S2C grey'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = true;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.02:0.01:0.5;
        optsOut.temporalWindowStepSizeOpts = 0.02:0.01:0.5;
        optsOut.cellTypes = {1,3,4,5,6,7,[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = false;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case 'S2C black'
        optsOut.dataSetID = 0;
        optsOut.significantFlag = false;
        optsOut.odourSelection = 1:7;
        optsOut.timeRange = [0,4];
        optsOut.temporalWindowOpts = 0.02:0.01:0.5;
        optsOut.temporalWindowStepSizeOpts = 0.02:0.01:0.5;
        optsOut.cellTypes = {1,3,4,5,6,7,[1,3:7]};
        optsOut.repOpts = 1:5;
        optsOut.shuffleBetweenRepsFlag = false;
        optsOut.rndSeed = 1;
        optsOut.nPermutations = 10000;
    case 'custom'
        optsOut = struct();
    otherwise
        disp('Invalid figure ID.');
        optsOut = [];
end
end

function R = PCC_helper(xx,yy)
mat = corrcoef(xx,yy);
R = mat(1,2);
end