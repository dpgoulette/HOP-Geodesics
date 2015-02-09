function [DT,VV,VC,GoodEdges,GoodIndex,edge_alpha] =...
   HOPDataPrepare(Data, varargin)
%
%  issues - update comments.  added alpha complex here.  Dependent function
%  AlphaCellsSelect.
%
%  comment varargin and what it is doing.

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
%     IMPORTANT! If you have the variable DT is in your workspace, then it
%     will be pulled into this function and used.  It will not be
%     calculated again.  If DT, VV and VC all exist then all three will be
%     used to speed up the function.
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
%        edge_alpha  -- k by 3 array.  The first two columns hold all of
%              the edges in the complete Delaunay triangulation except for
%              edges that include points on the boundary of the data
%              space (i.e. "bad data," see details below).  The third
%              column holds the value of alpha for the edge in that row.
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
%     that include bad points are removed.  This collection of edges makes
%     up the 1-skeleton of the Delaunay that we will "HOP" on. (Actually
%     this is a subset of the edges in the Delaunay triangulation because
%     we throw out boundary edges with bad data points.)
%
%     In the second section of this function, the user has the option to
%     further reduce the Delaunay 1-skeleton by calculating an alpha
%     complex.  The value of alpha is calculated for each edge in the
%     complex.  See AlphaCellsSelect for more on how this function works.

% Error check.  Make sure the data is 2d or 3d
if size(Data,2) < 2 || size(Data,2) > 3
   error('The data must be an Nx2 matrix of points or an Nx3 matrix.')
end

% check whether this is first run with raw data or a rerun of a previous
% run or not.
first_run = isempty(varargin);

if ~first_run
   fprintf('The variable last_run_data was detected in your workspace.\n')
   fprintf('So HOPGeodesics will use this data to do a rerun of the model.\n')
end

%%%%%%%%%%%%%%%%% MAIN: SECTION 1 %%%%%%%%%%%%%%%%%%%%%%

% Create the delaunay triangulation object and the voronoi cell data if it
% is not already in the workspace.  If DT is found in the workspace then
% this will be used instead of recalculating the DT object.  If VV and VC
% are in the workspace, the Voronoi diagram will not be recalculated.
if first_run
   fprintf('\nMaking the Delaunay triangulation object... ')
   DT=DelaunayTri(Data);
   fprintf('\nDone.\n')
   
   fprintf('\nNow calculating the Voronoi diagram.  This may take a while\n')
   [VV, VC]=voronoiDiagram(DT);
   fprintf('\nDone.\n\n')
   
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
   
else
   % this is a rerun so we don't have to recalculate everything.
   DT = varargin{1}.DT;
   VV = varargin{1}.VV;
   VC = varargin{1}.VC;
   GoodIndex = varargin{1}.GoodIndex;
   if isempty(varargin{1}.edge_alpha)
      GoodEdges = varargin{1}.GoodEdges;
   else
      GoodEdges = varargin{1}.edge_alpha(:, [1,2]);
   end
end

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
% returned will be a subset of the GoodEdges that is input.  edge_alpha is
% a k by 3 array.  It holds all of the original GoodEdges (that were sent
% to the AlphaCellsSelect function) in the first two columns, and the the
% value of alpha in the third column.
if alpha_option == 1
   % Then we don't need to calculate the alpha complex.
   edge_alpha = [];
else % We will HOP on an alpha complex 1-skeleton
   [GoodEdges, edge_alpha] = AlphaCellsSelect(DT,GoodEdges,VV,VC,GoodIndex);
   fprintf('\nFinished selecting the alpha complex\n')
end

fprintf('Finished the initial prep of the data.\n')

end % main function HOPDataPrepare


