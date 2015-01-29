function [maxpointer, minpointer, hop] = HighDensityHop(hop,GoodIndex,DT)
%   hoppointer - Runs the HOP algorithm to create a pointer vector.  This
%       vector contains the highest density neighbor of each point in the
%       data set.  This vector makes it easy to calculate hop paths.  The
%       function also catalogs relevent information in the hop data struct.
%       Points that are not "GoodPoints" (i.e. they are on the boundary of
%       the data set) point to NaN.  Maximums and minimums point to Inf.
%       All of the maximums and minimums are marked in the hop data struct.
%
%   Inputs:  
%       hop       - the hop data structure
%       GoodIndex - the index of "good" points that are not on the boundary
%                   of the data space.
%       n         - the total number of points in the data set
%   
%   This version of HOP hops to the highest density neighbor. This
%   was our original hop approach.

maxpointer=NaN(length(DT.X),1);
minpointer=NaN(length(DT.X),1);

[hop(:).ismax]=deal(false);
[hop(:).ismin]=deal(false);

%find maximums; A max has density larger than any of it's neigbors. mark in
%hop struct which points are maxs; create pointer list for the HOP paths;
%there will be a NaN at an entry for a bad point. There will be an Inf if
%the point is a max.  Note that if we do hop on an alpha complex, then some
%good points have no edges attached to them.  Thus they are maxima AND
%minima by definition.
for z=1:length(GoodIndex)
    if isempty(hop(GoodIndex(z)).edges) || ...
          hop(GoodIndex(z)).edges(end,3)<hop(GoodIndex(z)).density
        %GoodIndex(z) is a max
        hop(GoodIndex(z)).ismax = true;
        maxpointer(GoodIndex(z))=Inf;
    else%GoodIndex(z) is not a max
        maxpointer(GoodIndex(z))=hop(GoodIndex(z)).edges(end,2);
    end
end

%find minimums;mark in hop struct;create pointer list; there will be a NaN
%at an entry for a bad point.  There will be an Inf if the point is a min.
for z=1:length(GoodIndex)
    if isempty(hop(GoodIndex(z)).edges) || ...
          hop(GoodIndex(z)).edges(1,3)>hop(GoodIndex(z)).density
        %GoodIndex(z) is a min
        hop(GoodIndex(z)).ismin = true;
        minpointer(GoodIndex(z))=Inf;
    else%GoodIndex(z) is not a max
        minpointer(GoodIndex(z))=hop(GoodIndex(z)).edges(1,2);
    end
end
end% HighDensityHop function