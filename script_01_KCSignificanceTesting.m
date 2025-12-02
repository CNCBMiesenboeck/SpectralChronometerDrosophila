%% To perform the significance detection on the KC ephys odour response data. 
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
% This script performs the significance detection on the e-phys recordings
% of Kenyon cells in response to odour pulses. The mathematical details are
% presented in the publication by Kropf, Talbot and Miesenboeck. The
% significance detection is modified from Pouzat et al (2015) as described 
% in the publication.
%
% This script will output a summary of the results to the command window. It
% is configured to run on the data, and using the same settings, as in the 
% publication. To run it, ensure that the file "KC_data.mat", or
% "10s_KC_data.mat" is in the subfolder "data" on the current matlab path,
% then hit run, F5, or Shift-Enter.
%
% In case you want to use this as a basis for classifying your own data, or
% maybe reprocessing our data, please read on.
% 
% Have a look at the readme for the description of the variables in the
% data files. You will have to organise your data equivalent to the data 
% variables "odourOn"/"odourOff" and "baseline". 
%
% If you want to alter any settings, they'll most likely be in the 
% "Classification parameters to modify" section.
%
% Another modification you may want to make is to the criteria for 
% persistent response detection; a response to an odour onset may persist 
% beyond the offset edge, in which case we don't classify it as an off 
% response. I perform this classification within the two sections labelled 
% "Persistent response categorisation"; one is after the for loops and the 
% other is within the function "sigOdours" at the end of this file. You can 
% implement your own criteria within these sections. This will only be of
% concern if recordings post odour offset are being analysed.
%
% Another modification I can think of is to use a different noise 
% stabilisation transformation. I have included the option to switch to the
% alternative methods indicated by Pouzat et al. 
%
% The function "sigOdours" performs the classification of individual
% responses. It performs and outputs the classification of both excitatory 
% and inhibitory responses, but this work only uses the excitatory 
% responses. It should be easy to modify it to access the classification of
% inhibitory responses.
% 
% Once the data is processed, inspect the variables "sigOut", "ppSig", and 
% "currentParamsSig" for more detail than the text output to the command 
% window; these will contain arrays corresponding to individual cell-odour 
% responses. Please see comments below for details of what the variables 
% are.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the data set:
dataSetID = 0;  % enter 0 or 1 for the 4-second, or the 10-second recordings

% clear the workspace, screen, and load the data
clearvars -except dataSetID
clc
if dataSetID==0
    data = load('data\KC_data.mat', 'odourOn', 'odourOff', 'baseline', ...
        'on_significant', 'off_significant', 'anatomicalNames');
    
else
    data = load('data\10s_KC_data.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dataset parameters:
% Please read the readme for a description of the data organisation.
nKCTypes = numel(data.odourOn);     % number of KC subtypes.
KCTypeInds = 1:nKCTypes;            % the KC subtypes to test (all of them).
nCells = cellfun(@(x) size(x,1),data.baseline);   % number of cells in each subtype
nOdours = size(data.odourOn{1},2); % number of stimulii (including air).
nReps = size(data.odourOn{1},3);
ctrlInd = 1;        % the stimulus to use as the control (1=air).
if dataSetID ==0
    % for the 2-second data:
    odourInds = 3:9;    % the stimulii to test (1=air, 2=mineral oil).
    onOffInds = [1,2]; % 1=during odour pulse, 2=after odour pulse
else
    odourInds = [3,5:9]; % there was a problem with odour 4 in these recordings.
    onOffInds = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classification parameters to modify:
if dataSetID==0
    % odour response time relative to the on/off edge for the 4-second data
    % NB 4-seconds recordings = 2-s post odour onset plus 2-s post offset
    % and recall that there was a small amount of jitter in the length of
    % the pulse, hence the time window is truncated to 1.8 seconds (see
    % Methods of the publication)
    responseTimeWindow = [0,1.8];
else
    % or for the 10 second data.
    responseTimeWindow = [0,10];
end
baselineTimeWindow = [-3.5,-0.5];  % the time used for baseline 
alpha = 0.05;           % significance threshold
targetMean = 1;         % the expected number of spikes per bin under spontaneous activity
reps = 1:nReps;         % the odour repetitions to use in the significance detection
stabilisationInd = 1;   % the choice of stabilisation transformation (see below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% some different noise stabilization transformations. In the publication,
% we use option 1:
switch stabilisationInd
    case 1
        hStabilised = @(histIn) 2*sqrt(histIn+0.25);
    case 2
        hStabilised = @(histIn) sqrt(histIn)+sqrt(histIn+1);
    case 3
        hStabilised = @(histIn) 2*sqrt(histIn+0.375);
end
rt2 = sqrt(2);
hStabilisedDifference = @(odour,air) (hStabilised(odour)-hStabilised(air))/rt2;

% sigOut and ppSig are each a cell array of logical arrays. Each cell 
% corresponds to a KC subtype. The logical arrays will be nCell by nOdours 
% in size, indicating if the response crosses the (excitatory) significance
% threshold (sigOut), and if the response is potentially persistent
% (ppSig).
[sigOut,ppSig] = deal(cell(numel(onOffInds),nKCTypes));
% compute some constants for multiple use
responseRange = responseTimeWindow(2)-responseTimeWindow(1);
baselineRange = baselineTimeWindow(2)-baselineTimeWindow(1);
% define some useful helper functions 
hInRange = @(tt,range) tt(tt>range(1) & tt<range(2));
hOdRange = @(tt) hInRange(tt,responseTimeWindow);
hBgRange = @(tt) hInRange(tt,baselineTimeWindow);
hSpont = @(cc) sum(cellfun(@(x) numel(hBgRange(x)),cc))/(nReps*baselineRange);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the data
for onOff = onOffInds % loop through on and off responses
    onFlag = onOff==1;
    if onFlag
        testResponses = data.odourOn;
    else
        testResponses = data.odourOff;
    end
    for tt = KCTypeInds    % loop through KC type
        [sigOut{onOff,tt},ppSig{onOff,tt}] = deal(false(nCells(tt),nOdours));
        for cc = 1:nCells(tt)   % loop over cell number
            airResponse = hOdRange(cell2mat([testResponses{tt}{cc,ctrlInd,reps}]));
            spontAir = hSpont([data.baseline{tt}{cc,ctrlInd,reps}]);
            for oo = odourInds  % loop through odours
                spontOdour = hSpont([data.baseline{tt}{cc,oo,reps}]);
                spont = (spontOdour+spontAir)/2; % average spontaneous rate
                binWidth = targetMean/(nReps*spont); % ideal bin width
                nBins = max(1,floor(responseRange/binWidth)); % number of bins must be an integer >= 1
                bins = (0:nBins)*responseRange/nBins; % bin edges
                cellOdourResponse = hOdRange(cell2mat([testResponses{tt}{cc,oo,reps}]));
                odourPSTH = histcounts(cellOdourResponse,bins);
                airPSTH = histcounts(airResponse,bins);
                stabilisedDiff = hStabilisedDifference(odourPSTH,airPSTH);
                [sigOut{onOff,tt}(cc,oo),ppSig{onOff,tt}(cc,oo)] = ...
                    sigOdours(stabilisedDiff,alpha,onFlag);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Persistent response categorisation
% need to exclude the responses to the odour onset that persist into the 
% odour offset responses.
if dataSetID==0
    persSig = arrayfun(@(cc) ppSig{1,cc}&ppSig{2,cc},KCTypeInds,'UniformOutput',false);
    newOffSig = arrayfun(@(cc) sigOut{2,cc}&~persSig{cc},KCTypeInds,'UniformOutput',false);
    currentParamsSig = [sigOut(1,:);newOffSig];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section outputs a summary of the results.
hReshape = @(vals,arr) reshape(vals,size(arr,1),[]);
hSqueeze = @(arr) hReshape(arrayfun(@(ii) squeeze(arr{ii}(:,:,1)),1:numel(arr),...
    'UniformOutput',false),arr);
hSumArr = @(sig) hReshape(cellfun(@(x) sum(x,'all'),sig),sig);
hComp = @(hOp,AA,BB) hReshape(arrayfun(@(cc) hOp(AA{cc},BB{cc}),1:numel(AA),...
    'UniformOutput',false),AA);
hAAnotBB = @(AA,BB) hComp(@(a,b) a&~b,AA,BB);
if dataSetID==0
    publishedParamsSig = hSqueeze([data.on_significant;data.off_significant]);
    newNotOld = hAAnotBB(currentParamsSig,publishedParamsSig);
    oldNotNew = hAAnotBB(publishedParamsSig,currentParamsSig);
    respBoth = arrayfun(@(cc) sigOut{1,cc}&sigOut{2,cc}&~persSig{cc},KCTypeInds,...
        'UniformOutput',false);

    fprintf('\n\n\n----------------------------------------------------------------------------------------------------\n');
    disp('Numbers of significant responses');
    disp('Columns show different KC types:');
    fprintf('%s\t',data.anatomicalNames{:});fprintf('\n');
    disp('Top rows = ON odour, Bottom rows = OFF odour');
    disp('Using the parameters currently in the script:');
    disp(hSumArr(currentParamsSig));
    disp('Number of persistent responses:');
    disp(hSumArr(persSig));
    disp('Number of responses to both ON and OFF that aren''t persistent:');
    disp(hSumArr(respBoth));
    disp('Using the parameters given in the publication:');
    disp(hSumArr(publishedParamsSig));
    disp('Number of responses present with the current parameters and not with the published:');
    disp(hSumArr(newNotOld));
    disp('Number of responses present with the published parameters and not with the current:');
    disp(hSumArr(oldNotNew));
    fprintf('----------------------------------------------------------------------------------------------------\n');
else
    fprintf('\n\n\n----------------------------------------------------------------------------------------------------\n');
    disp('Numbers of significant cell-odour responses');
    disp('Using the parameters currently in the script:');
    disp(hSumArr(sigOut));
    disp('Using the parameters given in the publication:');
    disp(hSumArr(data.on_significant));
    disp('Number of responses present with the current parameters and not with the published:');
    disp(hSumArr(hAAnotBB(sigOut,data.on_significant)));
    disp('Number of responses present with the published parameters and not with the current:');
    disp(hSumArr(hAAnotBB(data.on_significant,sigOut)));
    fprintf('----------------------------------------------------------------------------------------------------\n');
end

%%
function [excit,PPexcit,inhib,PPinhib] = sigOdours(stabilisedDiff,sigTh,onFlag)
% please read the methods section of the paper for details.
nBins = numel(stabilisedDiff);
erfThresh = (1-sigTh)^(1/nBins);
rtN = sqrt(1:max(2,nBins));
% devMat = the normalised cumulative sums of stabilisedDiff.
% each row corresponds to the different start times of the sum.
% the normalisation is to the significance threshold.
devMat = zeros(nBins);
if onFlag
    binInds = nBins:-1:1;
else
    binInds = 1:nBins;
end
for kk = binInds % kk is the bin in which the cumulative sum starts
    currentBoundary = rtN(2)*erfinv(erfThresh^(1/(rtN(nBins+1-kk))));
    tot = 0;
    for jj = kk:nBins
        tot = tot+stabilisedDiff(jj);
        currMax = tot/(currentBoundary*rtN(jj+1-kk));
        devMat(kk,jj)=currMax;
    end
end
excitMat = devMat>=1;
inhibMat = devMat<=1;
excit = any(excitMat,'all');
inhib = any(inhibMat,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Persistent response categorisation
% To categorise the response as potentially being to the odour onset that
% then persists into the offset.
PPexcit = false;
PPinhib = false;
if excit||inhib
    % classify responses as possibly persistent
    if onFlag
        % for ON responses: if the last bin of any of the cumulative sums
        % is significant, and the deviation is in the correct direction.
        PPexcit = any(devMat(excitMat(:,end),end)>1);
        PPinhib = any(devMat(inhibMat(:,end),end)<1);
    else
        % for OFF responses: if the cumulative sum starting at the first
        % time point reaches the significance threshold and the deviation
        % in the first bin is in the correct direction.
        [firstExcitSigRow, ~] = ind2sub([nBins,nBins],find(excitMat,1));
        [firstInhibSigRow, ~] = ind2sub([nBins,nBins],find(inhibMat,1));

        PPexcit = (~isempty(firstExcitSigRow)&&firstExcitSigRow==1&&devMat(1,1)>0);
        PPinhib = (~isempty(firstInhibSigRow)&&firstInhibSigRow==1&&devMat(1,1)<0);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end