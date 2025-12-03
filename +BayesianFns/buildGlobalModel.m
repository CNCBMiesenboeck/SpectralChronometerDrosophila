function [hFcnOut, hardLimitsOut, description, paramsMap] = ...
    buildGlobalModel(dataOdour, hardLimitsIn, timeRange, config)
% build the global model by linking the parameters of individual models
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
% Inputs:
%   dataOdour 
%       A 2-element cell array, containing the movement initiation times of
%       the flies in response to the early- and late-trained odours.
%   hardLimitsIn 
%       A 2x2x6 array of the two sets of hardLimits for the base models.
%       models. The first index is the odour number. The second index is
%       the lower/upper limit and the third index is the parameter.
%   config 
%       determines the global model:
%           0   all parameters are global
%           1   all parameters are independent
%           2   independent amps, global tCentres, global sigmas
%           3   global amps, independent tCentres, independent sigmas
%
% Outputs:
%   hFcnOut
%       The function handle of the posterior (log-)PDF of the global model.
%       The number of parameters it accepts is dependent on the particular
%       global model configuration. Let's call the array of parameters into
%       hFcnOut "freeParams", and the number of parameters "nFreeParams".
%   hardLimitsOut
%       A 2-by-nFreeParams array containing the lower and upper limits of
%       the parameter search area. nFreeParams is dependent on config.
%   description
%       A string describing the global model.
%   paramsMap
%       A 2-by-6 array mapping freeParams to the separate calls to the base
%       function (two-component Gaussian) for each of the odour responses 
%       in the two elements of dataOdour. The columns correspond to the 
%       parameters of each model (2 each of: time centres, widths, and 
%       amplitudes).
%
% The global model is a combination (sum of log-probability) of the base
% models with some parameters fitted simultaneously. The base models are 
% the multi-component Gaussian models of the distribution of movement
% initiation times in response to the individual odours. The parameters 
% that are fitted simultaneously are refered to as global parameters. The
% global model will be a function handle (hFcnOut) which takes freeParams
% and distributes them to the base models according to paramsMap and then
% sums the output log-probability of the base models. The pseudo code for
% hFcnOut is essentially:
%   hFcnOut = @(P) fcnOut(P,earlyData,lateData,paramsMap)
%   function logP = fcnOut(freeParams,earlyData,lateData,paramsMap)
%       earlyParams = freeParams( paramsMap(top row) )
%       lateParams  = freeParams( paramsMap(bottom row) )
%       prob_earlyData = 2-component Gauss with earlyParams at earlyData time points
%       prob_lateData  = 2-component Gauss with lateParams at lateData time points
%       logP = log( product of prob_earlyData * product of prob_lateData )
%   end

% It is easy to write an algorithm to generate paramsMap. But for this 
% limited set of models for a pre-defined number of base models, defining 
% paramsMap manually is not prohibitive and is easier to understand:
switch config
    case 0
        description = 'Fully linked model (all parameters global)';
        paramsMap = [...
            1,  2,  3,  4,  5,  6; ...
            1,  2,  3,  4,  5,  6];
    case 1
        description = 'Independent model (all parameters independent, uniform limits)';
        paramsMap = [...
            1,  2,  5,  6,  9, 11; ...
            3,  4,  7,  8, 10, 12];
    case 2
        description = 'Intuitive model (global tCentres and sigmas, independent amps)';
        paramsMap = [...
            1,  2,  3,  4,  5,  6; ...
            1,  2,  3,  4,  7,  8];
    case 3
        description = 'Dummy model (independent tCentres, independent sigmas, global amps)';
        paramsMap = [...
            1,  2,  5,  6,  9, 10; ...
            3,  4,  7,  8,  9, 10];
end
nFreeParams = max(paramsMap,[],'all');

% The search area of the global-model parameters is the maximal of the 
% separate structure-models for the two trained odours. This is only
% actually relevant for the amplitude parameters; the upper limit of the 
% search area is the greater of the two independent search areas. 
hardLimitsOut = zeros(2,nFreeParams);
for pp = 1:nFreeParams
    ll = inf;
    uu = -inf;
    for mm = 1:2
        test = paramsMap(mm,:)==pp;
        if any(test)
            ll = min(ll,hardLimitsIn(mm,1,test));
            uu = max(uu,hardLimitsIn(mm,2,test));
        end
    end
    hardLimitsOut(1,pp) = ll;
    hardLimitsOut(2,pp) = uu;
end

hFcnOut = @(pp) sum(arrayfun(@(ii) sum(log(BayesianFns.multiModeGauss(dataOdour{ii},...
    pp(paramsMap(ii,1:2)),pp(paramsMap(ii,3:4)),pp(paramsMap(ii,5:6)),timeRange,...
    hardLimitsOut(:,paramsMap(ii,:))))),1:2));
end