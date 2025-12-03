function [logVolOut, volOut, logStdOut, stdOut] = ...
    estimatePDFVolume(hTheta, thetaSamples, parameterLimits, rndSeed, parallelFlag)
% estimate the volume of the log-density function hTheta (NB LOG-density).
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
%   hTheta 
%       - the handle to the log-PDF
%   thetaSamples 
%       - an N-by-P array of samples drawn from hTheta (using e.g. slice
%         sampling sampling). N is the number of samples and P is the 
%         number of parameters taken by hTheta.
%   parameterLimits 
%       - a 2-by-P array with the first row defining the lower limit and
%         the second row defining the upper limit of the search area.
%   rndSeed 
%       - initialises the random number generator to allow for exact
%         reproduction of any results. Enter -1 for no seeding, or an 
%         integer >= 0.
%   parallelFlag 
%       - (boolean) determines whether to use the parallel pool.
%
% Outputs:
%   volOut 
%       - the estimate of the volume of hTheta.
%   logVolOut 
%       - log(volOut)
%   stdOut
%       - the standard deviation for the estimate of the volume.
%   logStdOut
%       - log(stdOut)
%
% the estimation uses equation 19 of the publication:
% volume ~= mean(exp_hTheta(phiSamples)) / mean(exp_hPhi(thetaSamples))
%   where exp_hPhi is a normalised function (i.e. integrates to unity) with 
%   coverage over the same parameter range as hTheta. In this function, we
%   use a normalised Gaussian. phiSamples will be drawn from exp_hPhi. The 
%   notation exp_hTheta and exp_hPhi indicates that these functions are the 
%   linear counterparts to the log-density functions hTheta and hPhi.

% This function has been simplified for the purposes of this publication.
% hTheta must be the logpdf. This function will not generate thetaSamples.
% This function will not automatically determine parameterLimits.

% increasing the number of phiSamples compensates for any possible zones of
% limited coverage between hTheta and hPhi. NB because hPhi is Gaussian,
% generating phiSamples is fast.
nTheta = size(thetaSamples, 1);
nPhi = 5*nTheta;

% calculate the means, standard deviations, and amplitudes of the Gaussians
% defining hPhi. The amplitudes are adjusted using the error function to 
% compensate for chopping off the Gaussian wings with parameterLimits.
mus = mean(thetaSamples);
stds = std(thetaSamples);
rt2 = sqrt(2);
amps = 2./(stds*rt2*sqrt(pi).*(...
    erf((parameterLimits(2, :)-mus)./(stds*rt2)) - ...
    erf((parameterLimits(1, :)-mus)./(stds*rt2))));
logCombAmps = sum(log(amps));

% draw samples from hPhi. See function generateSamples below for details
phiSamples = generateSamples(nPhi, parameterLimits, mus, stds, rndSeed);
% evaluate hTheta at each row of phiSamples. parArrayFun calls hTheta
% multiple times using parallel processing if parallelFlag is true.
thetaFunc_phiSamples = parArrayFun(hTheta, phiSamples, 'parallel processing', parallelFlag);

% evaluate log(hPhi) at each row of thetaSamples
% This is quite quick. Using parallel processing was slower, presumably 
% because of the overheads when transferring data between pools.
vars = 2*stds.*stds;
phiFunc_thetaSamples = zeros(nTheta, 1);
for tt = 1:nTheta
    if all(thetaSamples(tt,:)<=parameterLimits(2,:)) && ...
            all(thetaSamples(tt,:)>=parameterLimits(1,:))
        phiFunc_thetaSamples(tt) = ...
            logCombAmps-sum((thetaSamples(tt, :)-mus).^2./vars);
    else
        phiFunc_thetaSamples(tt) = -inf;
    end 
end

nThetaFunc_PhiSamples_NZ = sum(~isinf(thetaFunc_phiSamples));
nPhiFunc_ThetaSamples_NZ = sum(~isinf(phiFunc_thetaSamples));

% logToLinear "un-logs" probabilities for evaluating means and errors. The
% values are also scaled for compatibility with floating point precision
[expThetaPhi, thetaPhi_Scale] = logToLinear(thetaFunc_phiSamples);
[expPhiTheta, phiTheta_Scale] = logToLinear(phiFunc_thetaSamples);
rescaleFactor = (thetaPhi_Scale-phiTheta_Scale);
meanThetaPhi = mean(expThetaPhi);
varThetaPhi = var(expThetaPhi);
meanPhiTheta = mean(expPhiTheta);
varPhiTheta = var(expPhiTheta);
% log version of equation 19 in the publication.
logVolOut = rescaleFactor + log(meanThetaPhi) - log(meanPhiTheta);
volOut = exp(logVolOut);
% derived using propagation of errors.
logStdOut = sqrt(varPhiTheta/(nPhiFunc_ThetaSamples_NZ*(meanPhiTheta^2)) + ...
    varThetaPhi/(nThetaFunc_PhiSamples_NZ*(meanThetaPhi^2)));
stdOut = volOut*logStdOut;
end

function [scaledArr, logMax] = logToLinear(arr)
% many values of exp(arr) are likely to be too small for double-precision
% floats. This computes a exp(arr)
logMax = max(arr);
scaledArr = exp(arr-logMax);
end

function samples = generateSamples(nPhi, parameterLimits, mus, stds, seed)
% samplesOut is an nPhi-by-nParams array of random numbers drawn from
% Gaussian distributions with means and standard deviations (mus and stds)
% truncated by parameterLimits.
% The size of the parameterLimits must be 2-by-nParams. The sizes of mus
% and stds are 1-by-nParams.
% If seed >=0, the random number generator is seeded for reproducibility
if seed >= 0
    rng(seed, 'threefry');
end
nParams = size(parameterLimits, 2);
% transform parameterLimits to the domain of the error function
erfLim = erf((parameterLimits-repmat(mus, 2, 1))./repmat(stds*sqrt(2), 2, 1));
erfRange = erfLim(2, :) - erfLim(1, :);
% Generate samples using inverse transform sampling in the erf domain:
% 1. Create uniform random values between erfLim(1) and erfLim(2)
% 2. Apply the inverse error function to obtain standard normal samples
% 3. Scale (A) and shift (B) the samples back to the mus and stds
samples = repmat(mus, nPhi, 1) + ...        (step 3B)
    repmat(stds*sqrt(2), nPhi, 1) .* ...    (step 3A)
    erfinv(...                              (step 2)
    rand(nPhi, nParams).*repmat(erfRange, nPhi, 1)+repmat(erfLim(1,:), nPhi, 1)); % (step 1)
end