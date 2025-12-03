function [maxBoutStarts,sortedOrder] = getMaxBoutStarts(dataIn, absThreshold, minBoutLen)
% To compute the most vigorous movement bout of each fly and sort them
% according to the time point of initiation.
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
% The most vigorous movement bout is that in which the individual moved the
% most (i.e. the greatest intentional distance covered). The details are
% described in the publication. A lower threshold is applied to the raw 
% movement traces to remove frame-to-frame camera noise and segment them
% into bouts. A temporal threshold is then applied to remove momentary 
% flicks of the ball. The starting time index of the most vigorous bouts
% are sorted in increasing order. Individuals that exhbited no above 
% threshold movement are assigned a starting time index of -1 and are 
% sorted by their total (below-threshold) movement.
%
% Inputs:
%   dataIn          
%       - array containing the movement of the flies. Rows contain the 
%         movement traces of individual flies.
%   absThreshold    
%       - the camera noise threshold, applied to dataIn.
%   minBoutLen      
%       - the temporal threshold (in units of array elements).
%
% Outputs:
%   maxBoutStarts   
%       - the starting index of the most vigorous movement bout of each
%         fly, sorted as previously described.
%   sortedOrder     
%       - the order of the individuals in dataIn, mapping to maxBoutStarts.

absData = abs(squeeze(dataIn));
nFlies = size(absData,1);
maxBoutStarts = -ones(nFlies,1);
maxBoutVig = zeros(nFlies,1);
for rr = 1:nFlies
    % apply initial threshold and extract movement bouts:
    absThreshInds = absData(rr,:)>=absThreshold;
    % to segment the thresholded movement into bouts:
    [starts,ends,boutLengths]=lIndsToStartsEnds(absThreshInds);
    % remove flicks of the ball that fall below the temporal threshold:
    flickInds = boutLengths<minBoutLen;
    flickThreshInds = absThreshInds;
    for ff = find(flickInds)
        flickThreshInds((starts(ff)+1):ends(ff)) = false;
    end
    flickThreshCurve = absData(rr,:).*flickThreshInds;
    % re-extract the movement bouts:
    [starts,ends]=lIndsToStartsEnds(flickThreshInds);
    % find the movement bout with the most accumulated movement. NB If 
    % there are no movement bouts, the maxBoutStarts(bb) remains at -1.
    for bb = 1:numel(starts)
        boutMag = sum(flickThreshCurve((starts(bb)+1):ends(bb)));
        if boutMag>maxBoutVig(rr)
            maxBoutVig(rr)=boutMag;
            maxBoutStarts(rr)=starts(bb);
        end
    end
end
% sort the bouts
[~,sortedOrder] = sort(maxBoutStarts);
% for flies that exhibit no above threshold movement, sort them by the
% integrated (sub-threshold) movement.
zeroBoutInds = maxBoutStarts==-1;
zeroVig=sum(absData(zeroBoutInds,:),2);
[~,zeroOrder]=sort(zeroVig);
sortedOrder(1:numel(zeroVig))=sortedOrder(zeroOrder);
maxBoutStarts = maxBoutStarts(sortedOrder);
end

function [starts,ends,boutLengths]=lIndsToStartsEnds(logicalInds)
% logicalInds is a 1D logical array. This function segments the logical
% true phases of logicalInds and returns the starting and ending indicies,
% and the length of the bouts. 

changes = logicalInds(2:end)-logicalInds(1:(end-1));
starts = find(changes==1);
ends = find(changes==-1);
if isempty(starts)&&isempty(ends)
    boutLengths = [];
    return;
end
if isempty(starts) || (~isempty(ends) && starts(1)>ends(1))
    starts = [1, starts];
end
if isempty(ends) || ends(end)<starts(end)
    ends = [ends, numel(logicalInds)];
end
boutLengths = ends-starts;
end