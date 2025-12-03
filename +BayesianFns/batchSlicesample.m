function samplesOut = batchSlicesample(initial, nsamples, varargin)
% This simply calls the Matlab function slicesample using parallel
% processing unless otherwise specified, with support for seeded random
% number generation. Please first see the Matlab help page for slicesample;
% this function accepts all arguments as slicesample plus the extra 
% optional Name-Value arguments:
%   'nWorkers'
%       - number of parallel workers to use (default: 6)
%   'seed'
%       - base seed for RNGs to ensure reproducibility (default: 1)
%   'parallel'
%       - whether to run in parallel (default: true)
%   'burnin each'
%       - whether each worker performs its own burn-in (default: false)
%
% The returned array is as described in slicesample.

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

% Parse optional arguments
args = {'burnin', 'nWorkers', 'seed', 'parallel', 'burnin each'};
defaults = {0, 6, 1, true, false};
[vals, eargs] = parseVarargin(args, defaults, varargin{:});
[burnin, nWorkers, rootSeed, parallelFlag, burninEachFlag] = deal(vals{:});

if parallelFlag
    if burnin>0
        if burninEachFlag
            % If each worker will do its own burn-in from the same initial
            % point
            initSample = initial;
            eargs = [eargs, 'burnin', burnin];
        else
            % Or if the burn-in is done shared (done once before parallel
            % sampling)
            initSample = seededSliceSample(rootSeed, initial, 1, ...
                'burnin', burnin, eargs{:});
        end
    else
        initSample = initial;
    end

    % partition sampling among workers.
    nPartition = floor(nsamples/nWorkers);
    nExtra = nsamples-nPartition*nWorkers;
    fObjs(1:nWorkers) = parallel.FevalFuture;
    % launch seeded slice sampling on each worker. See the function 
    % seedeSliceSample below.
    for ww = nWorkers:-1:1
        fObjs(ww) = parfeval(@seededSliceSample, 1, ...
            rootSeed+(rootSeed>=0)*ww, initSample, ...
            nPartition+(ww==nWorkers)*nExtra, eargs{:});
    end
    % wait for workers to complete, then collect samples
    wait(fObjs);
    samplesOut = zeros(nPartition*nWorkers, numel(initial));
    for ww = 1:nWorkers
        samplesOut((ww-1)*nPartition+(1:(nPartition+(ww==nWorkers)*nExtra)), :) = ...
            fObjs(ww).fetchOutputs;
    end
else
    % if not using parallel processing, simply call the seeded slice
    % sampling function. See the function seededSliceSample below.
    eargs = [eargs, 'burnin', burnin];
    samplesOut = seededSliceSample(rootSeed, initial, nsamples, eargs{:});
end
end

function samplesOut = seededSliceSample(seed, init, nSamplesOut, varargin)
% wrapper for slicesample that sets the random number generator seed before
% running.
if seed>=0
    rng(seed, 'threefry');
end
samplesOut = slicesample(init, nSamplesOut, varargin{:});
end