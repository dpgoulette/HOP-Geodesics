function [maxpointer, minpointer, hop] = GradientHop(hop,GoodIndex,DT)
%   GradientHop - Runs the HOP algorithm to create a pointer vector.  For
%       each point, the vector contains the neighbor which is in the
%       direction of greatest gradient. This vector makes it easy to
%       calculate hop paths.  The function also catalogs relevent
%       information in the hop data struct. Points that are not
%       "GoodPoints" (i.e. they are on the boundary of the data set) point
%       to NaN.  Maximums and minimums point to Inf. All of the maximums
%       and minimums are marked in the hop data struct.
%
%
%   Inputs:  
%       hop       - the hop data structure
%       GoodIndex - the index of "good" points that are not on the boundary
%                   of the data space.
%       n         - the total number of points in the data set
%   
%   This version of HOP hops to the neighbor of highest gradient. 

maxpointer=NaN(length(DT.X),1);
minpointer=NaN(length(DT.X),1);

% Initialize each point as not being a max or a min.
[hop(:).ismax]=deal(false);
[hop(:).ismin]=deal(false);

%find maximums; A max has density larger than any of it's neigbors. mark in
%hop struct which points are maxs; create pointer list for the HOP paths;
%there will be a NaN at an entry for a bad point. There will be an Inf if
%the point is a max.


for i=1:length(GoodIndex)
    % This next section calculates the rate of change, Delta, with each
    % neighbor of the point i.
    %
    % Delta = (difference in density)/(distance between the points)
    % 
    % The gradient is the greatest positive Delta among the neighbors.
    
    NbsID = hop(GoodIndex(i)).edges(:,2);%indices of i's neighbors
    NbsCoords = DT.X(NbsID,:);%coordinates of i's neigbors
    
    %The next two lines calculate the distance from i to its nbors. 
    %NbsCoords is k x 2  and  DT.X(GoodIndex(i)) is 1 x 2.      
    diffs = bsxfun(@minus,DT.X(GoodIndex(i),:),NbsCoords); %differences
    Dist = sqrt(sum(diffs.^2,2)); %distances
    
    %find the density differences
    DenseDiff = bsxfun(@minus,...
        hop(GoodIndex(i)).edges(:,3),hop(GoodIndex(i)).density);
    
    %Deltas; the rate of change to each neigbor of i
    Deltas = DenseDiff./Dist;
    
    %Now find the positive gradient neighbor and catalog in pointer vector.
    [DeltaMax, DeltaMaxID] = max(Deltas);%Max Delta with Index
    
    % Now catalog whether each point is a max or not and store it in the
    % hop data structure.  If hop is done on an alpha complex then some
    % points may have no edges attached to them so in this case DeltaMax is
    % empty and the point is a max AND a min by definition.  If DeltaMax is
    % negative then all neighbors have negative gradient, so the point is a
    % max.
    if isempty(DeltaMax) || DeltaMax <= 0 
        %then GoodIndex(i) is a max
        hop(GoodIndex(i)).ismax = true;
        maxpointer(GoodIndex(i))=Inf;
    else%then GoodIndex(i) is not a max
        maxpointer(GoodIndex(i)) = hop(GoodIndex(i)).edges(DeltaMaxID,2);
    end
    
    
    %Now do the negative gradient neighbor for HOP to the mins
    [DeltaMin, DeltaMinID] = min(Deltas);
    
    if isempty(DeltaMin) || DeltaMin >= 0
        %then GoodIndex(i) is a min
        hop(GoodIndex(i)).ismin = true;
        minpointer(GoodIndex(i)) = Inf;
    else%GoodIndex(i) is not a min
        minpointer(GoodIndex(i)) = hop(GoodIndex(i)).edges(DeltaMinID,2);
    end
end
end% HighDensityHop function
