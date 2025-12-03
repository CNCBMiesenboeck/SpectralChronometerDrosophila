%% To reproduce the Bayesian analysis of the behaviour data
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
% This script perfroms the Bayesian modelling of the behaviour data. With
% the default settings, and running with four processes/threads, the code
% should reproduce exactly the results in the publication. With different
% settings, the results should still be approximately the same. The Bayes
% Factors are output to the command window when running the model
% evaluation. The figures from the publication, and more, can also be
% produced.
%
% This script is designed to be run in sections. To run a section,
% highlight it (left click within the section) and then hit Ctrl-Enter
%
% Execute sections 1 and 2 to load the data and set up various parameters.
% The logical variable "testingFlag" in section 2 can be set to true to run
% the modelling at a "very low resolution" in order to test that the code
% runs. 
% Section 3.1 evaluates the 2-component Gaussian models applied to the
% behavioural responses to individual odours (quantify if there is 
% structure to the distributions). This can be run once before running
% sections 3.2 and/or 3.3 with different configurations.
% Section 3.2 and 3.3 produce plots of the output of section 3.1. Please
% see the comments within these sections for how to reproduce the plots.
% Sections 4.1 - 4.3 are analogous to section 3, but for the global models.
%% Section 1
% to load the preprocessed data (i.e. movement initiation times). The
% script "script_behaviour_maxBoutStartTimes.m" demonstrates the conversion
% of the raw data to the movement initiation times and histograms.
clear
clc
dataFileName = 'data\behaviourResponses.mat';
loaded = load(dataFileName);

%% Section 2
% Select the data set and establish the settings for the numerical Bayesian
% estimation of the behavioural models.
% The data sets are the movement initiation times of the flies trained with
% optogentic DAN stimulation, or the optogenetic control as described in the
% publication.
dataSetInd = 1;          % 1=DAN, 2=chrimson control
testingFlag = false;     % set to true to test that the code runs properly
% define a struct to contain the settings for the slice sampling:
modelEvalSettings = struct();
if testingFlag
    modelEvalSettings.slicesamplesPerParam = 100;
    modelEvalSettings.burninPerParam = 10;
    modelEvalSettings.areaSamplesPerParam = 100;
    modelEvalSettings.burninEachFlag = true;
    modelEvalSettings.parallelFlag = false;
    modelEvalSettings.nWorkers = 1;
else
    modelEvalSettings.slicesamplesPerParam = 10000; % publication:  10000
    modelEvalSettings.burninPerParam = 100;         %               100
    modelEvalSettings.areaSamplesPerParam = 2000;   %               2000
    modelEvalSettings.burninEachFlag = true;        %               true
    modelEvalSettings.parallelFlag = true;          %               true
    modelEvalSettings.nWorkers = 6;                 %               6
end
modelEvalSettings.rndSeed = 1;                      %               1

epochNames = {'Early', 'Late', 'Neutral', 'Binned'};
dataSetOpts = {'DAN stim','Chr Ctrl'};
switch dataSetInd
    case 1
        dataStruct = loaded.DANData;
    case 2
        dataStruct = loaded.chrCtrlData;
end
timeRange = dataStruct.movtHistTimeScale([1,end]);
timeDivider = 5.5; % the cutoff between the two Gaussians

% Prepare data structures, parameter search limits, and initialise output
% variables. Do this here so that data set/model options can be run 
% individually without overwritting the others.
%
% for the structure modelling (section 3):
nOdours = 3;
[dataOdour, hMdlSglOdour, structureSamples] = ...
    deal(cell(1,nOdours));
[logNullBaseModelsLiklihood, logBaseModelsLiklihood, logBaseModelsError] = ...
    deal(zeros(1,nOdours));
% for the global modelling (section 4):
nGlobalModels = 4;
[hGlob, hardLimitsGlob, globDescription, paramsMap, globSamples] = ...
    deal(cell(1,nGlobalModels));
[logGlobalModelLiklihood, logGlobalModelError] = deal(zeros(1,nGlobalModels));

% get the data and parameter search ranges for each odour. This is required
% for both sections 3 and 4, but only needs to be done once.
for odInd = 1:nOdours
    [dataOdour{odInd},hardLimits(odInd)] = getDataEpoch(...
        dataStruct,odInd,timeRange,timeDivider);
end

% define helper function to synchronise axes limits across array of hAxes
hAxLimsHelper1 = @(hAxes,lims) [min(arrayfun(@(h) min(get(h,lims)),hAxes)),max(arrayfun(@(h) max(get(h,lims)),hAxes))];
hSetLimsHerlper2 = @(hAxes,range,lims) arrayfun(@(h) set(h,lims,range),hAxes);
hSetAxLims = @(hAxes,lims) hSetLimsHerlper2(hAxes,hAxLimsHelper1(hAxes,lims),lims);
%% Section 3.1
% To test if there is structure to the data by "fitting" Gaussian models to
% the individual data sets.
% This secion outputs the Bayes factors to the command window. These are those
% shown in figure 5B (for dataSetInd=1) and extended 3C (dataSetInd=2).
% This section will take a few minutes to run, depending on your computer
% specs and the modelEvalSettings.
for odInd = 1:nOdours
    % the null model is evenly distributed data points, so the 
    % log-probability of the null model giving rise to the data is simply:
    logNullBaseModelsLiklihood(odInd) = ...
        -numel(dataOdour{odInd})*log(timeRange(2)-timeRange(1));
    % define function handle of the model-data function that takes a set of
    % parameters and computes the probability that the data would be 
    % obtained with those parameters.
    % See also file "multiModeGauss.m"
    hMdlSglOdour{odInd} = @(p) sum(log(BayesianFns.multiModeGauss(dataOdour{odInd}, ...
        p(1:2),p(3:4),p(5:6),timeRange,hardLimits(odInd).array)));
    % perform slice sampling and compute the posterior marginal liklihood
    % (see also function evaluateModel in section 5 below)
    [structureSamples{odInd}, logBaseModelsLiklihood(odInd), logBaseModelsError(odInd)] = ...
        evaluateModel(hMdlSglOdour{odInd}, hardLimits(odInd).array, ...
        modelEvalSettings);
    % compute the Bayes Factor and print to screen:
    logBayesFactor = logBaseModelsLiklihood(odInd)-logNullBaseModelsLiklihood(odInd);
    fprintf('%s data set. %s odour. logBF (compared to null): %.4f +/- %.4e\n', ...
        dataSetOpts{dataSetInd}, epochNames{odInd}, logBayesFactor, ...
        logBaseModelsError(odInd));
end
%% Section 3.2
% To plot the parameter histograms corresponding to the individual models.
% None of the plots produced by this section are shown in the paper.
odInd = 1;  % 1 = early-reinforced odour, 2 = late-reinforced, 3 = non-reinforced
parameterInd = 1; % 1 = Gaussian centres, 2 = widths, 3 = amplitudes
% See comments in function plotHistogram in section 5 below.
plotHistogram(structureSamples{odInd},parameterInd,odInd);

%% Section 3.3
% To reproduce the movement initiation "fits" shown in figure 5:
%   A, dataSetInd=1, odInd=1 (lower-middle panel) or odInd=2 (lower-right)
% Figures produced with other parameter combinations are not shown in the
% paper, but are producible.
odInd = 1; % 1 = early-reinforced odour, 2 = late-reinforced, 3 = non-reinforced
% See comments in function plotFit in section 5 below.
CI = [2.5,97.5]; % confidence interval (shaded region)
plotFit(structureSamples{odInd}, ...
    dataStruct.movtHistInit(odInd,:),dataStruct.movtHistTimeScale,CI,odInd);
%% Section 3.4
% To produce "fits" of the null model (i.e. constant probability) shown in
% figure 5:
%   A, dataSetInd=1, odInd=3 (lower-left panel)
% and figure extended 3:
%   B, dataSetInd=2, odInd=1 (lower-middle), 
%                    odInd=2 (lower-right),
%                    odInd=3 (lower-left)
odInd = 3;
CI = [2.5,97.5]; % confidence interval (shaded region)
plotNullModel(sum(dataStruct.sortedTimes(odInd,:)>dataStruct.times(2)),...
    dataStruct.movtHistInit(odInd,:),dataStruct.movtHistTimeScale,CI,odInd);
%% Section 4.1
% To test if, and which of, the parameters are independent. Corresponding
% to figure 5, panels C-E.
% This section will take several minutes to run, depending on your PC specs
% and the modelEvalSettings.
% Bayes Factors are output to the command window. The two comparisons are
% to the null model (uniform movement initiation) and to the fully linked 
% model (Global model 0). The latter are those presented in figure 5D.

% The following is a combination of the parameter limits required to
% compute the global model parameter search area:
combParamLims = permute(reshape([hardLimits(1:2).array],2,6,2),[3,1,2]);

nGlobalFlies = sum(cellfun(@numel,dataOdour(1:2)));
% The probability of obtaining the data if the responses are randomly 
% distributed. This is not used in the paper. It is only really useful for
% comparison with global model 0 (see below):
logNullGlobModelLiklihood = -nGlobalFlies*log(timeRange(2)-timeRange(1));

% globalModelConfig:
%   0   all parameters are global
%   1   all parameters are independent
%   2   independent amps, global tCentres, global sigmas
%   3   global amps, independent tCentres, independent sigmas
for globalModelConfig = 0:(nGlobalModels-1)
    gg = globalModelConfig+1;
    % obtain the posterior PDF handle and parameter search limits required
    % for the model evaluation. Also obtain string describing the model and
    % the paramsMap required for plotting the output. (See file 
    % buildGlobalModel.m for details):
    [hGlob{gg}, hardLimitsGlob{gg}, globDescription{gg}, paramsMap{gg}] = ...
        BayesianFns.buildGlobalModel(dataOdour(1:2), combParamLims, timeRange, ...
        globalModelConfig);
    % perform slice sampling and compute the posterior marginal liklihood
    % (see also function evaluateModel in section 5 below):
    [globSamples{gg}, logGlobalModelLiklihood(gg), logGlobalModelError(gg)] = ...
        evaluateModel(hGlob{gg}, hardLimitsGlob{gg}, modelEvalSettings);
    % compute Bayes factor(s) and print to screen:
    if gg==1
        logBinnedModelVolume = logGlobalModelLiklihood(gg);
        fprintf(['Global model %d. %s\n\t Log volume: %.4f.\t',...
            'log BF (compared to null): %.4f\n'], ...
            globalModelConfig, globDescription{gg}, logBinnedModelVolume, ...
            logBinnedModelVolume-logNullGlobModelLiklihood);
    else
        logBayesFactor = logGlobalModelLiklihood(gg)-logBinnedModelVolume;
        fprintf(['Global model %d. %s\n\t logBF (compared to null): ',...
            '%.4f\n\t logBF compared to fully linked: %.4f +/- %.4e\n'], ...
            gg-1, globDescription{gg}, ...
            logGlobalModelLiklihood(gg)-logNullGlobModelLiklihood, ...
            logGlobalModelLiklihood(gg)-logBinnedModelVolume, ...
            logGlobalModelError(gg)+logGlobalModelError(1));
    end
end
%% Section 4.2
% To plot the parameter histograms corresponding to the global models shown
% in figure 5:
%   C, dataSetInd=1, ggInd=1, parameterInd=1
%   E, dataSetInd=1, ggInd=2, parameterInd=3
% Plots can be produced with other parameter combinations, but are not shown
% in the paper.
ggInd = 1; % global model selector (see globalModelConfig in section 4.1)
parameterInd = 1; % 1 = Gaussian centres, 2 = widths, 3 = amplitudes
odIndOpts = 1:2; % leave this as is. 1 = early-reinforced odour, 2 = late-reinforced.
% set up the figure
figure('OuterPosition',[50,50,numel(odIndOpts)*560,500]);
hAxTemp = gobjects(numel(odIndOpts),1);
for axInd = 1:numel(odIndOpts)
    odInd = odIndOpts(axInd);
    hAxTemp(axInd) = subplot(1,numel(odIndOpts),axInd);
    % See comments in function plotHistogram in section 5 below.
    plotHistogram(globSamples{ggInd+1},parameterInd,odInd,...
        paramsMap{ggInd+1},globDescription{ggInd+1},hAxTemp(axInd));
end
hSetAxLims(hAxTemp,'YLim');
%% Section 4.3
% To plot the movement initiation "fits" of the global models to the
% individual odour responses (similar to lower panels of figure 5A).
% None of the plots produced by this section are shown in the paper.
% As a guide, have a look at:
%   set A: the null model fits              (ggInd=0)
%   set B: the fully independent model fits (ggInd=1)
%   set C: the intuitive model fits         (ggInd=2)
%   set D: the dummy model fits             (ggInd=3)
% You can see that the data histograms lie much better within the bounds of
% sets B and C compared to set A. Set B may encompass more of the histograms
% than set C, but the latter (the intuitive model) fits tighter around more
% of the histogram. The intuitive model (set C) is more parsimonious, and so
% it outperforms the fully independent model (comparing the marginal
% probabilities as in fig 5D.) The dummy model encompasses the data
% histogram, but is very broad compared to set B or set C, which is why its
% marginal probability is the weakest.
ggInd = 1; % global model selector (see globalModelConfig in section 4.1)
CI = [2.5,97.5]; % confidence interval (shaded region)
odIndOpts = 1:2; % leave this as is. 1 = early-reinforced odour, 2 = late-reinforced.
% set up the figure
figure('OuterPosition',[50,50,1200,500]);
hAxTemp = gobjects(numel(odIndOpts),1);
for odInd = odIndOpts
    hAxTemp(odInd)=subplot(1,2,odInd);
    % Extract the relevant columns of globSamples using the relevant row of 
    % paramsMap. Please see comments within the file buildGlobalModel.m for
    % more detail on paramsMap.
    % See comments in function plotFit in section 5 below.
    plotFit(globSamples{ggInd+1}(:,paramsMap{ggInd+1}(odInd,:)), ...
        dataStruct.movtHistInit(odInd,:),dataStruct.movtHistTimeScale,CI,odInd,...
        globDescription{ggInd+1},hAxTemp(odInd));
end
hSetAxLims(hAxTemp,'YLim');
%% Section 5. Helper functions
function [sortedSamples, logModelLikelihood, logLiklihoodError, ...
    paramSearchArea, sortedLogPs, markovInds] = ...
    evaluateModel(hModel, paramLimits, modelEstimationSettings)
% helper function to prepare for and perform the slice sampling, and obtain
% the likelihood of the model giving rise the data.
% Inputs:
%   hModel        - handle to the posterior PDF, accepting nParams 
%                   parameters and returning the log-probability that the
%                   data would be obtained with those parameters.
%   paramLimits   - the maximum limits of the parameter search area.
% 	modelEstimationSettings 
%                 - the settings to configure the slice sampling. See file
%                   batchSlicesample.m for further details.
%
% Outputs:
%   sortedSamples - the samples of the posterior PDF, sorted in order of
%                   descending probability.
%   logModelLikelihood
%                 - the (log of the) PDF integrated over, and normalised
%                   by, the supported parameter space within paramLimits
%   paramSearchArea
%                 - the hyper area of the supported parameter space within
%                   paramLimits
%   sortedLogPs   - the output of hModel for each of sortedSamples
%   markovInds    - the indicies of the the sortedSamples as they were
%                   generated by the (Markovian) slice sampling procedure.

import BayesianFns.*

% slice sampling can be biased if the starting sample is far from the peak
% of the posterior PDF. We compute the initial point using MLE using the
% function computeInitialGuesses below (see comments within).
initGuesses = computeInitialGuesses(hModel, paramLimits, modelEstimationSettings.rndSeed);

nParams = numel(initGuesses);
nSliceSamples = modelEstimationSettings.slicesamplesPerParam*nParams;
% batchSlicesample simply calls Matlab's slicesample function using
% parallel processing. See comment in file batchSlicesample.m
samplesOut = batchSlicesample(initGuesses, ...
    nSliceSamples, ...
    'logpdf', hModel, ...
    'width', diff(paramLimits), ...
    'burnin', modelEstimationSettings.burninPerParam*nParams, ...
    'seed', modelEstimationSettings.rndSeed, ...
    'parallel', modelEstimationSettings.parallelFlag, ...
    'burnin each', modelEstimationSettings.burninEachFlag, ...
    'nWorkers', modelEstimationSettings.nWorkers);

% estimate the hyper-area of the supported parameter space within
% paramLimits. See comments in file estimatePDFSupportArea.m
paramSearchArea = estimatePDFSupportArea(hModel, paramLimits, ...
    'logpdf flag', true, ...
    'nSamples', modelEstimationSettings.areaSamplesPerParam*nParams, ...
    'rootSeed', modelEstimationSettings.rndSeed, ...
    'parallel', modelEstimationSettings.parallelFlag);

% estimate the (log) hyper-volume of the posterior PDF within paramLimits.
[logModelVolume,~,logLiklihoodError] = estimatePDFVolume(hModel, samplesOut, paramLimits, ...
    modelEstimationSettings.rndSeed, modelEstimationSettings.parallelFlag);
% normalise the PDF hyper-volume by the parameter space hyper-area to
% obtain the model likelihood
logModelLikelihood = logModelVolume - log(paramSearchArea);

% evaluate the PDF at each sample, sort the values, and retain the order of
% the Markovian sampling.
logPs = zeros(nSliceSamples,1);
for ss = 1:nSliceSamples
    logPs(ss) = hModel(samplesOut(ss,:));
end
[sortedLogPs,sortInds] = sort(logPs,'descend');
[~,markovInds] = sort(sortInds);
sortedSamples = samplesOut(sortInds,:);
end

function maxParams = computeInitialGuesses(hFun, hardLimitsArray, rndSeed)
% maximises hFun using the Levenberg-Marquardt algorithm and returns the 
% parameters of hFun at that point.
%
% Inputs:
%   hFun            - handle to a function that accepts nParams inputs and
%                     returns a single value. 
%   hardLimitsArray - the parameter search space. A 2-by-nParams array with
%                     the top and bottom rows containing the lower and
%                     upper limits of each parameter.
%   rndSeed         - the seed value (>= 0) of the random number generator.
%                     Enter -1 to ignore.
% Output:
%   initGuesses     - the parameters at the maximum of hFun.

% A valid initial guess is required for the Levenberg-Marquardt algorithm. 
% Estimate the initial guesses as the centre of the parameter ranges.
initGuesses = mean(hardLimitsArray);
% If the initial guess is not within the supported region of hFun, then
% repeatedly generate a random point within the parameter seach space until
% a valid point is found.
if isinf(hFun(initGuesses))
    if rndSeed >= 0
        rng(rndSeed, 'threefry');
    end
    range = diff(hardLimitsArray);
    nParams = numel(range);
    paramsFail = true;
    while paramsFail
        initGuesses = hardLimitsArray(1,:) + rand(1,nParams).*range;
        paramsFail = isinf(hFun(initGuesses));
    end
end

% settings used for the fitting algorithm. See Matlab help page for
% optimset for details.
lmeOpts = optimset('Algorithm','levenberg-marquardt','MaxFunEvals', 1000,...
    'MaxIter',1000,'TolX',10e-6,'TolFun',1e-6,'display','off');
% run the Levenberg-Marquardt algorithm. See Matlab help page for lsqnonlin
% for further details.
maxParams = lsqnonlin(@(x) -hFun(x),initGuesses,[],[],lmeOpts);
end

function [data, hardLimits] = getDataEpoch(dataStruct, epochInd, timeRange, timeDivider)
% helper function to extract the valid movement initiation times and 
% parameter search limits from dataStruct.
%
% Inputs:
%   dataStruct  - struct with elements:
%       sortedTimes  - containing the initiation time of the most vigorous
%                      movement bout for all flies in a data set in 
%                      response to the three odours.
%       times        - containing the time points of data acquisition
%                      (camera frames).
%       movtHistTimeScale and movtHistInit 
%                    - containing the x and y data for the smoothed
%                      movement initiation histograms.
%   epochInd    - the index to select odour response. 1, 2 and 3 correspond
%                 to early reinforced, late reinforced, and non-reinforced
%                 respectively.
%   timeRange   - a 2-element array containing the lower and upper limits
%                 of the model function support range. 
%   timeDivider - the timepoint midway between the early punishment offset
%                 and the late punishment onset.
%
% Outputs:
%   data        - the valid movement initiation times. 
%   hardLimits  - 2-by-6 array containing the lower and upper limits
%                 of the parameters of the model PDF.

% dataStruct.times starts one time point before the odour onset. Flies with
% no above-threshold movement are assigned a start time corresponding to
% the first index of dataStruct.times. Flies that have their most
% significant movement bout starting before the odour onset will have a
% start time corresponding to the second index of dataStruct.times. 
% We exclude both of these sets of flies (as described in the publication). 
validMoverInds = dataStruct.sortedTimes(epochInd,:)>dataStruct.times(2);
data = dataStruct.sortedTimes(epochInd,validMoverInds)';

% Define the boundaries of the parameter search. It is called hardLimits
% because even within the boundaries of the search area, some parameter
% combinations will be disallowed (if the total area of the double Gaussian
% distribution is greater than one).

% The time limits of the centre of the two Gaussians are:
% 1.5 s -> timeDivider
% timeDivider -> odour offset
hardLimits.timeLimits = [1.5, timeDivider; timeDivider, timeRange(2)]';
hardLimits.timeLimits = [0, timeDivider; timeDivider, timeRange(2)]';
% hardLimits.timeLimits = [0, timeRange(2); 0, timeRange(2)]';
% The limits of the widths of the Gaussians are 0.25 -> 2.0
hardLimits.sigmaLimits = [0.25, 2.0; 0.25, 2.0]';
hardLimits.sigmaLimits = [0.0, 5; 0.0, 5]';

% The upper limit of the amplitude parameter is defined as the maximum of 
% the smoothed histogram (stored in dataStruct.movtHistInit) in the 
% allowed time range of the Gaussian peak, plus one Poisson error based on
% the number of flies with movement initiation within that time range.
%
% First define helper function (handles).
% To extract in-range logical time indicies of the movement initiation 
% histogram:
hPMovTimeInds = @(rng) dataStruct.movtHistTimeScale>rng(1) & ...
    dataStruct.movtHistTimeScale<=rng(2); 
% To extract the maximum height of the histogram in range:
hMaxHistHeight = @(rng) max(dataStruct.movtHistInit(epochInd,hPMovTimeInds(rng)));
% To normalise the histogram and add one Poisson error:
hMaxAmp = @(rng) (hMaxHistHeight(rng)/numel(data))*...
    (1+1/sqrt(sum(data>rng(1)&data<rng(2))));
% Set the amplitude limits:
hardLimits.ampLimits = [...
    0.0, hMaxAmp(hardLimits.timeLimits(:,1)); ...
    0.0, hMaxAmp(hardLimits.timeLimits(:,2))]';

% Concatenate all parameter limits into a single array:
hardLimits.array = [...
    hardLimits.timeLimits,...
    hardLimits.sigmaLimits,...
    hardLimits.ampLimits];
end

function plotFit(samples, movtHist, timeScale, CI, odInd, modelDesc, hAx)
% helper function to plot the "fitted" 2-component Gaussians and confidence
% bounds of the sampled posterior PDF, overlaid onto the smoothed histogram
% of movement initiation times.
%
% Inputs:
%   samples    - an nSamples-by-6 array of (extracted) samples of the 
%                posterior PDF. The samples must have been sorted from most
%                to least likely. nSamples is the number of slice samples
%                (determined outside this function). The columns correspond
%                to the parameters of a 2-component Gaussian (centres,
%                widths, and amplitudes) modelling the response to
%                individual odours. In the case of the global models, the
%                relevant columns of samples must be extracted before 
%                calling this function (please see section 4.3 and the 
%                comments within the file buildGlobalModel.m).
%   movtHist   - the smoothed histogram of the movement initiation times.
%   timeScale  - the time points to plot (the x-axis).
%   CI         - the lower and upper percentiles of the distributions.
%                Default is 2.5 and 97.5.
%
% This function assumes that the top row of samples contains the most
% likely set of parameters, which is plotted in a solid line. The upper and
% lower bounds are computed as follows:
%   randomly sub-sample from rows of samples (default=1000)
%   compute the shape of the 2-component Gaussian for each sub-sample
%   at each point on the time-axis, select the lower and upper bounds using
%     the percentiles in CI.
import BayesianFns.*
if nargin < 4
    CI = [2.5,97.5];
end
if nargin < 5
    odInd = 1;
end
if nargin < 6
    modelDesc = 'Structural model (independent parameters, independent limits)';
end
switch odInd
    case 2
        fitColour = [0.7,0,0];
        odourStr = 'late-reinforced odour';
    case 3
        fitColour = [0.5,0.5,0.5];
        odourStr = 'non-reinforced odour';
    otherwise
        fitColour = [1,0.5,0];
        odourStr = 'early-reinforced odour';
end
nSubSample = inf;
nSamplesIn = size(samples,1);
if nSamplesIn >= nSubSample
    subSamples = samples(randperm(nSamplesIn,nSubSample),:);
else
    subSamples = samples;
end
timeScale = reshape(timeScale,[],1);
nSamples = size(subSamples,1);
nTimePoints = numel(timeScale);
distributions = zeros(nSamples, nTimePoints);
for ss = 1:nSamples
    distributions(ss,:) = multiModeGauss(timeScale, ...
        subSamples(ss,1:2),subSamples(ss,3:4),subSamples(ss,5:6), ...
        timeScale([1,end]),repmat([-inf;inf],1,6));
end
nBest = max(ceil(size(subSamples,1)/100),1);
hBest = @(cols) mean(samples(1:nBest,cols),1);
bestFit = multiModeGauss(timeScale, ...
    ...mean(samples(1:nBest,1:2),1),mean(samples(1:nBest,3:4),1),mean(samples(1:nBest,5:6),1),...
    hBest(1:2),hBest(3:4),hBest(5:6),...
    timeScale([1,end]),repmat([-inf;inf],1,6));
bounds = [prctile(distributions,CI(1),1);prctile(distributions,CI(2),1)];
if nargin<7 || ~isa(hAx,'matlab.graphics.axis.Axes') || ~isvalid(hAx)
    hAx = axes(figure(), 'NextPlot', 'add');
else
    hAx.NextPlot = 'add';
end
patch(hAx,[timeScale;timeScale(end:-1:1)],[bounds(1,:),bounds(2,end:-1:1)], ...
    [0.75,0.75,0.75],'FaceAlpha',0.5,'EdgeColor','none');
plot(hAx,timeScale, ...
    ((nTimePoints-1)/(timeScale(end)-timeScale(1)))*movtHist/sum(movtHist), ...
    'Color',[0,0,0],'LineWidth',3);
plot(hAx,timeScale,bestFit,'Color',fitColour,'LineWidth',3);
title(sprintf('%s\n%s',modelDesc,odourStr));
end

function plotNullModel(nMovers, movtHist, timeScale, CI, fitColourInd)
if nargin < 4
    CI = [2.5,97.5];
end
if nargin < 5
    fitColourInd = 1;
end
switch fitColourInd
    case 2
        odourStr = 'late-reinforced odour';
        fitColour = [0.7,0,0];
    case 3
        odourStr = 'non-reinforced odour';
        fitColour = [0.5,0.5,0.5];
    otherwise
        odourStr = 'early-reinforced odour';
        fitColour = [1,0.5,0];
end
baseline = 1/(timeScale(end)-timeScale(1));
UL = baseline*poissinv(CI(2)/100,nMovers)/nMovers;
LL = baseline*poissinv(CI(1)/100,nMovers)/nMovers;
hAx = axes(figure(), 'NextPlot', 'add');
patch(hAx,timeScale([1,end,end,1]),[LL,LL,UL,UL],...
    [0.75,0.75,0.75],'FaceAlpha',0.5,'EdgeColor','none');
plot(hAx,timeScale, ...
    ((numel(timeScale)-1)/(timeScale(end)-timeScale(1)))*movtHist/sum(movtHist), ...
    'Color',[0,0,0],'LineWidth',3);
plot(hAx,timeScale([1,end]),baseline*[1,1],'Color',fitColour,'LineWidth',3);
title(sprintf('Null model, %s',odourStr));
end

function plotHistogram(samples,parameterInd,odInd,parameterMap,modelDesc,hAx)
% helper function to plot the marginal histograms of individual parameters.
%
% Inputs:
%   samples - the samples of the posterior PDF.
%   parameterInd
%           - to select (1) Gaussian centres, (2) widths, or (3) amplitudes
%   odInd   - for the global models, to select the parameters corresponding
%             to the (1) early- or (2) late-reinforced odour.
%   parameterMap
%           - the mapping of global model parameters to individual odours.
%             please see comments in file buildGlobalModel.m for more
%             details. Do not enter if plotting samples produced from the
%             base models.
%   modelDesc
%           - a string describing the model (title for the plot). Do not
%             enter if plotting samples produced from the base models.
if nargin < 4
    parameterMap = repmat(1:6,3,1);
    modelDesc = 'Structural model (independent parameters, independent limits)';
end
switch odInd
    case 1
        odourStr = 'early-reinforced odour';
        colourOrder = [1,0.5,0;0.5,0.5,0.5];
    case 2
        odourStr = 'late-reinforced odour';
        colourOrder = [0.5,0.5,0.5;0.7,0,0];
    case 3
        odourStr = 'non-reinforced odour';
        colourOrder = [0.5,0.5,0;0.0,0.5,0.5];
end
switch parameterInd
    % 1 for Gaussian centres, 2 for widths, 3 for amplitudes
    case 1
        titles = {'Gaussian centres','Time (s)','Density (/s)'};
        indsTemp = 1:2;
        limsTemp = [0,10];
    case 2
        titles = {'Gaussian widths','Width (s)','Density (/s)'};
        indsTemp = 3:4;
        limsTemp = [0,5];
    case 3
        titles = {'Gaussian peak height','Amplitude (/s)','Density (s)'};
        indsTemp = 5:6;
        limsTemp = [0,0.35];
end
axLims = limsTemp;
nSamples = size(samples,1);
nBins = 125;
edges = linspace(limsTemp(1),limsTemp(2),nBins+1);
dE = edges(2)-edges(1);
edges = [edges(1)-dE,edges,edges(end)+dE];
normConst = nBins/((limsTemp(2)-limsTemp(1))*nSamples);
if nargin<6 || ~isa(hAx,'matlab.graphics.axis.Axes') || ~isvalid(hAx)
    hAx = axes(figure(),'NextPlot','add', 'ColorOrder', ...
        colourOrder);
else
    hAx.NextPlot = 'add';
    hAx.ColorOrder = colourOrder;
end
for ii = indsTemp
    traceTemp = histcounts(samples(:,parameterMap(odInd,ii)),edges)*normConst;
    plot(hAx,edges(2:end)-((edges(2)-edges(1))/2),traceTemp,'LineWidth',2);
end
title(sprintf('%s\n%s, %s',modelDesc,titles{1},odourStr));
xlabel(titles{2});
ylabel(titles{3});
legend({'early component', 'late component'},'Location','best');
hAx.XLim = axLims;
end