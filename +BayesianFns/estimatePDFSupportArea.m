function [area, isInPDF] = estimatePDFSupportArea(hFun, paramLimits, varargin)
% estimate the hyper-area of a PDF support region within the search
% (hyper-)area defined by paramLimits using Monte Carlo sampling.
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
% This function estimates the hyper-area of a uniform prior distribution
% over a user-defined parameter space. The parameter space may contain
% regions where the model probability density function is zero. This
% function estimates the supported region.
%
% Inputs:
%   hFun
%       - a handle to the PDF/logPDF
%   paramLimits
%       - is a 2-by-nDim array of lower and upper limits of each parameter
%
% Optional Name-Value Inputs:
%   'nSamples'     
%       - number of Monte Carlo samples (default: 1000)
%   'logpdf flag'  
%       - true if hFun is a logPDF (default: true)
%   'rootSeed'     
%       - seed for RNG (default: -1, meaning no seeding)
%   'parallel'     
%       - whether to use parallel processing (default: false)
%   'nWorkers'     
%       - number of parallel workers (default: 4)
%
% Outputs:
%   area
%       - the hyper-area of the PDF support region
%   isInPDF
%       - an array of logical values corresponding to each Monte Carlo
%       sample inidicating if it was in the support area

% parse the optional arguments
args = {'nSamples', 'logpdf flag', 'rootSeed', 'parallel', 'nWorkers'};
defaults = {1000, true, -1, false, 4};
vals = parseVarargin(args, defaults, varargin{:});
[nSamples, logPDFFlag, rootSeed, parallelFlag, nWorkers] = deal(vals{:});

if parallelFlag
    % --- Parallel Execution Block ---
    % Distributes the work between nWorkers, which call this function (the 
    % serial execution block below) and then collates the outputs.
    
    % number of samples to generate per worker, rounding up for one worker
    % to compensate for integer rounding.
    nSamplesWorkers = floor(nSamples/nWorkers)*ones(nWorkers, 1);
    nSamplesWorkers(end) = nSamples - nSamplesWorkers(1)*(nWorkers-1);
    % preallocate futures.
    fObjs(1:nWorkers) = parallel.FevalFuture;
    % loop to dispatch parallel jobs (with parallel flag = default (false))
    for ww = 1:nWorkers
        fObjs(ww) = parfeval(@BayesianFns.estimatePDFSupportArea, 2, ...
            hFun, paramLimits, 'nSamples', nSamplesWorkers(ww), ...
            'logpdf flag', logPDFFlag, 'rootSeed', rootSeed);
        if rootSeed >= 0
            % use different seeds for each worker
            rootSeed = rootSeed+1;
        end
    end
    % wait for workers to complete and collect results
    wait(fObjs);
    isInPDF = false(nSamples,1);
    for ww = 1:nWorkers
        [~, isInPDF((ww-1)*nSamplesWorkers(1)+(1:nSamplesWorkers(ww)))] = ...
            fObjs(ww).fetchOutputs;
    end
    % compute the proportion of samples in the support region
    occupationFraction = mean(isInPDF);
else
    % --- Serial Execution Block ---
    
    % seed the random number generator if necessary
    if rootSeed >= 0
        rng(rootSeed, 'threefry');
    end

    % define function handle to test if a point is within the support
    % region of hFun.
    if logPDFFlag
        hLogicalFun = @(x) ~isinf(hFun(x));
    else
        hLogicalFun = @(x) hFun(x)>0;
    end
    
    nDim = size(paramLimits, 2);
    isInPDF = false(nSamples, 1);
    % Sample uniformly in the bounding hypersquare and test for support
    for ss = 1:nSamples
        isInPDF(ss) = hLogicalFun(...
            (paramLimits(2, :)-paramLimits(1, :)).*rand(1, nDim) + paramLimits(1, :));
    end
    % compute the proportion of samples in the support region
    occupationFraction = mean(isInPDF);
end
% convert occupationFraction to hyper-area.
area = occupationFraction * prod(diff(paramLimits));
end