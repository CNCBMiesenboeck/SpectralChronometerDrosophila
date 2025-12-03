function out = parArrayFun(hFun, arrayArg, varargin)
% Applies a function (hFun) to rows of an array (arrayArg) in parallel 
% unless otherwise specified. 
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
% hFun can accept optional arguments, but each call to hFun must use the 
% same options, and the Names of the optional arguments to hFun must be 
% different to the optional arguments to this function (i.e. parArrayFun).
% The number of outputs of hFun is limited to one.
%
% Inputs:
%   hFun        
%       - function handle that accepts a vector input as its first argument
%         and optional arguments in Name-Value format, which are passed in
%         Name-Value format to this function (ie parArrayFun). The output
%         of hFun must be a single value.
%   arrayArg    
%       - nIn-by-nArguments array. hFun will be called iteratively with
%         each row of arrayArg.
%
% Optional Name-Value Arguments:
%   'nWorkers'            
%       - number of parallel workers (default: 4)
%   'seed'                
%       - random number generator seed. -1 (default) indicates no seeding. 
%         This will be useful if hFun uses random numbers. The seeding is
%         performed before calling hFun.
%   'parallel processing' 
%       - flag to use parallelism (default: true)
%
% Output:
%   out
%       - column vector of outputs.

% parse optional inputs. Any Name-Value arguments not on this list will be
% stored in eargs and passed to each call of hFun
args = {'nWorkers', 'seed', 'parallel processing'};
defaults = {4, -1, true};
[vals, eargs] = parseVarargin(args, defaults, varargin{:});
[nWorkers, rootSeed, parallelFlag] = deal(vals{:});

nIn = size(arrayArg, 1);
out = zeros(nIn, 1);

if parallelFlag
    % partion the rows of arrayArg between nWorkers, adding rows to the
    % last worker to compensate for integer rounding.
    nPartition = floor(nIn/nWorkers) * ones(nWorkers, 1);
    nPartition(end) = nIn - nPartition(1)*(nWorkers-1);
    partitionInds = arrayfun(@(ww) (ww-1)*nPartition(1)+(1:nPartition(ww)), ...
        1:nWorkers, 'UniformOutput', false);
    % preallocate futures
    fObjs(1:nWorkers) = parallel.FevalFuture;
    % dispatch parallel jobs
    for ww = 1:nWorkers
        fObjs(ww) = parfeval(@hLoop, 1, hFun, rootSeed, ...
            arrayArg(partitionInds{ww}, :), eargs{:});
        if rootSeed >= 0
            rootSeed = rootSeed+1;
        end
    end
    % wait for workers to complete and collect results.
    wait(fObjs);
    for ww = 1:nWorkers
        out(partitionInds{ww}) = fObjs(ww).fetchOutputs;
    end
else
    % specify parallelFlag as false if for example debugging.
    out = hLoop(hFun, rootSeed, arrayArg, eargs{:});
end
end

function outLoop = hLoop(hFun, seed, partArrayArg, varargin)
% helper function to apply hFun to a subset of inputs.

% Seed the random number generator if necessary.
if seed>=0
    rng(seed, 'threefry');
end
% preallocate output and iteratively call hFun.
nOut = size(partArrayArg, 1);
outLoop = zeros(nOut, 1);
for ii = 1:nOut
    outLoop(ii) = hFun(partArrayArg(ii, :), varargin{:});
end
end