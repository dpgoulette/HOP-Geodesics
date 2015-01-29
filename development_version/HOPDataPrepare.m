function [DT,VV,VC,GoodEdges,GoodIndex] = HOPDataPrepare(Data)
%
%  issues - update comments.  added alpha complex here.  Dependent function
%  AlphaCellsSelect.

% HOPDataPrepare - prepares the raw 2D or 3D data for geometric and
%     topological analysis.  In particular, HOPDataPrepare creates a graph
%     using the data points as vertices.  The goal is to apply the HOP
%     algorithm to this graph.  The graph that is created is a subset of
%     the 1-skeleton of the Delaunay triangulation.  We throw away some
%     boundary points (and edges), and the user has the option of
%     calculating an alpha complex as well. (An alpha complex is a subset of
%     the Delaunay triangulation with certain geometric properties; we
%     really only need the 1-skeleton of this alpha complex for HOP.) For
%     more, see the "Details" section of this help section below.
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
%        GoodIndex -- m by 1 index into the DT.X array.  This index
%              represents the "good" points in the Data.  "Bad" points are
%              on the boundary of the data space and produce pathological
%              voronoi cells.
%        GoodEdges -- m by 2 array of indices into DT.X.  Each row in
%              GoodEdges is an edge in the Delaunay Triangulation. These
%              edges define the graph (1-skeleton) that we will use for the
%              HOP algorithm. GoodEdges is a subset of the Delaunay
%              Triangulation because we have removed edges that contain
%              "bad" data points on the boundary of the data space.  Also
%              the user has the option to select an alpha complex which
%              would be a smaller subset of the edges.
%
%
%     Note that VV, VC and DT are created by the built-in matlab function
%     DelaunayTri which will be removed in a future release of Matlab.  For
%     more info see:  http://www.mathworks.com/help/matlab/ref/delaunaytri.html
%     The use of DelaunayTri should be easy to update when that happens.
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Details:
%
%     In the first part of this function, the Delaunay triangulation and
%     the Voronoi diagram are both created.  "Bad data" is removed from the
%     data set. Bad data is data on the boundary of the data space.  This
%     data tends to have unreliable voronoi cell shapes. So any point with
%     a voronoi cell that is at least partially outside the convex hull of
%     the data space is removed from consideration in the model.  The edges
%     (1-cells) in the Delaunay triangulation are collected and the edges
%     with bad points are removed.  This collection of edges makes up the
%     1-skeleton of the Delaunay that we will "HOP" on. (Actually this is a
%     subset of the edges in the Delaunay triangulation because we throw
%     out boundary edges with bad data points.)
%
%     In the second section of this function, the user has the option to
%     further reduce the Delaunay 1-skeleton by calculating an alpha
%     complex.  See AlphaCellsSelect for more on how this function works.

%%%%%%%%%%%%%%%%% MAIN: SECTION 1 %%%%%%%%%%%%%%%%%%%%%%
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

% Finished with the initial prep of the data.  We have created the Delaunay
% triangulation, the Voronoi diagram, returned the edges, and we have
% thrown out boundary probelms.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MAIN: SECTION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we calculate the 1-skeleton of an alpha complex (if the user
% wants to do that).

% User selects whether they want to hop on the complete Delaunay 1-skeleton
% (minus the edges on the boundary of the data space) or an alpha complex
% (subset of Delaunay).
fprintf('\nWould you like to HOP on the full Delaunay 1-skeleton, or the\n')
fprintf('1-skeleton of an alpha complex?\n')
while true
   fprintf('   1) HOP on full Delaunay.\n')
   fprintf('   2) HOP on alpha complex.\n')
   alpha_option = input('Choose one of the above: ');
   if alpha_option == 2 || alpha_option == 1
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
      pause(1)
   end
end
fprintf('\n')

% If the user wants to HOP on the 1-skeleton of the alpha complex (a subset
% of Delaunay), then run the selection scheme.  The GoodEdges array that is
% returned will be a subset of the GoodEdges that is input.
if alpha_option == 1
   % Then we don't need to calculate the alpha complex.
   clear alpha_option
else % We will HOP on an alpha complex 1-skeleton
   clear alpha_option
   GoodEdges = AlphaCellsSelect(DT,GoodEdges,VV,VC,GoodIndex);
   fprintf('\nFinished selecting the alpha complex\n')
end

fprintf('Finished the initial prep of the data.\n')

end % main function HOPDataPrepare


