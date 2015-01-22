function [DT,VV,VC,GoodEdges,GoodIndex] = HOPDataPrepare(Data)
% HOPDataPrepare - prepares the raw 2D or 3D data for geometric and
%     topological analysis.
%
%     input:
%        Data -- k by 2 (or k by 3) array of coordinate points.
%
%     Output:
%        DT -- A Delaunay Triangulation object.  Note: DT.X contains the
%              Data array. DT.Triangulation is the triangulation
%              information.  
%        VV -- n by 2 (or n by 3) array of Voronoi Cell vertices.  n 
%        VC -- k by 1 cell of indices into VV.  The i-th entry in VC
%              corresponds to the i-th data point in DT.X. The i-th entry
%              contains the voronoi cell point indices for the i-th data
%              point.
%
%        GoodIndex -- m by 1 index into DT.X that are the "good" points in
%              the Data.  "Bad" points are on the boundary of the data
%              space and produce pathological voronoi cells.
%        GoodEdges -- m by 2 array of indices into DT.X. Each row in
%              GoodEdges is an edge in the Delaunay Triangulation.
%              GoodEdges is a subset of the Delaunay Triangulation because we
%              have removed edges that contain "bad" data points on the
%              boundary of the data space.
%
%     Note that VV, VC and DT are created by the built-in matlab function
%     DelaunayTri which will be removed in a future release of Matlab.  For
%     more info see:  http://www.mathworks.com/help/matlab/ref/delaunaytri.html

%Error check.  Make sure the data is 2d or 3d
if size(Data,2) < 2 || size(Data,2) > 3
   error('The data must be an Nx2 matrix of points or an Nx3 matrix.')
end

%Create the delaunay triangulation object and the voronoi cell data.
fprintf('\nMaking the Delaunay triangulation object... ')
DT=DelaunayTri(Data);
fprintf('\nDone.\n')
fprintf('\nNow starting the Voronoi diagram.  This may take a while,')
fprintf(' especially if the data is 3D.\n')
[VV, VC]=voronoiDiagram(DT);
fprintf('\nDone making the Delaunay triangulation and the associated,')
fprintf(' Voronoi cell diagram.\n')
fprintf('\n')

% In our work we choose to remove data that is on the boundary of the data
% space.  We do this by removing any data point that has a voronoi cell
% vertex which we deem to be a "bad Voronoi vertex."  Infinity is bad
% because any voronoi cell attached to it has infinite area/volume.  Also
% any voronoi vertex outside of the data space is bad because they are
% often artificially elongated, thus the voronoi cell is less reliable.
% Note, any finite voronoi vertex that is outside of the dataspace
% will not be inside any Delaunay tetra.

% find the indices of the bad Voronoi vertices.  pointLocation returns a nan
% if the voronoi vertex is NOT interior to any triangle/tetra.  We put a NaN
% in the first slot because the point at infinity is listed first in the
% Voronoi diagram.
ID = [NaN; pointLocation(DT,VV(2:length(VV),:))];
% create an index vector of bad voronoi verts
BadVV = find(isnan(ID));

% Now we find the bad data.  "Bad data" is on the boundary of the data space
% so it has a voronoi cell which includes a bad voronoi vertex.  Create an
% index of "BadDataID."
BadDataID = zeros(length(DT.X),1);
a=1;
for b=1:length(VC)
   if ~isempty(intersect(BadVV,VC{b}))
      BadDataID(a)=b;
      a=a+1;
   end
end
BadDataID(BadDataID==0)=[];

% Now we can use this BadDataID vector to remove any Tetras, Triangles and
% Edges that contain any of these bad points.  Then we will only be
% considering cells with GoodPoints for inclusion.

% get ALL Delaunay edges (including bad ones)
DTedges=edges(DT);

% Now remove bad edges.  They are "bad" if they contain a bad data point.
% NOTE!!  The "GoodEdges" matrix that is created here may be altered if
% the user chooses to HOP on an alpha complex instead of the full delaunay.
% GoodIndex contains the indices (into DT.X) of the good data.
[R,~]=find(ismember(DTedges,BadDataID));
R=unique(R);
GoodEdges=DTedges;
GoodEdges(R,:)=[];%delete the bad edges.
GoodIndex=unique(GoodEdges);

end % main function HOPDataPrepare


