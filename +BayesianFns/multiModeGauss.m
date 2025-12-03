function pTimePoints = multiModeGauss(dataTimePoints, tCentres, sigmas, amps, timeRange, hardLimits)
% Evaluates a normalized mixture of two truncated Gaussians.
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
%   dataTimePoints
%       - vector of time points at which to evaluate the PDF
%   tCentres       
%       - 1x2 vector of Gaussian means (centres)
%   sigmas        
%       - 1x2 vector of standard deviations
%   amps           
%       - 1x2 vector of amplitudes (heights) of each Gaussian
%   timeRange      
%       - [minTime, maxTime] limits for truncating the distribution
%   hardLimits     
%       - 2x6 matrix defining min and max bounds for each parameter
%
% Output:
%   pTimePoints    
%       - the value of the Gaussian mixture PDF at the given time points.
% % % Copyright (C) 2025 Clifford Talbot

% Check parameter validity and compute the baseline offset of the truncated 
% distribution
offset = checkLimitsAndComputOffset(tCentres,sigmas,amps,timeRange,hardLimits);

% A negative offset indicates that the areas of the two Gaussians sums to 
% more than one, meaning that the parameter combination is not within the 
% support region of the PDF.
if offset>=0
    % data in a single column, parameters in rows.
    
    % compute time difference from each Gaussian centre. Produces an array
    % with size nData-by-2. (nData = number of dataTimePoints).
    deltaT = (dataTimePoints-tCentres);
    % replicate the variances, to the same size array as deltaT:
    varArr = (ones(numel(dataTimePoints),1)*(sigmas.*sigmas));
    % evaluate the Gaussian mixture and add the offset
    pTimePoints = amps * exp(-(deltaT.*deltaT)./(2*varArr))' + offset;
    % truncate values outside the time range.
    pTimePoints(dataTimePoints<timeRange(1)|dataTimePoints>timeRange(2)) = 0;
else
    % invalid parameters. Hence probability = 0.
    pTimePoints = zeros(size(dataTimePoints));
end
end

function offset = checkLimitsAndComputOffset(tCentres, sigmas, amps, timeRange, hardLimits)
% check the limits. NB writing it out long-hand as follows is much faster
% than defining a function handle
TF = tCentres(1)>=hardLimits(1,1) && tCentres(1)<=hardLimits(2,1) && ...
    tCentres(2)>=hardLimits(1,2) && tCentres(2)<=hardLimits(2,2) && ...
    sigmas(1)>=hardLimits(1,3) && sigmas(1)<=hardLimits(2,3) && ...
    sigmas(2)>=hardLimits(1,4) && sigmas(2)<=hardLimits(2,4) && ...
    amps(1)>=hardLimits(1,5) && amps(1)<=hardLimits(2,5) && ...
    amps(2)>=hardLimits(1,6) && amps(2)<=hardLimits(2,6);
if TF
    % compute the areas of each truncated Gaussian using the error function
    sigRt2 = sigmas*sqrt(2);
    areasOut = sqrt(pi/2)*sigmas.*amps.*(...
        erf((timeRange(2)-tCentres)./sigRt2) - erf((timeRange(1)-tCentres)./sigRt2));
    % the offset accounts for the remaining area of the PDF (the PDF must
    % sum to unity). 
    offset = (1-sum(areasOut))/(timeRange(2)-timeRange(1));
else
    offset = -1;
end
end