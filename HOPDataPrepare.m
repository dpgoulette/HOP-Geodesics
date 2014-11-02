function [DT,VV,VC,BadDataID,DTedges,GoodIndex,...
                            GoodEdges] = HOPDataPrepare(Data)
% comment this and clean the rest up

%Error check.  Make sure the data is 2d or 3d
if size(Data,2) < 2 || size(Data,2) > 3
    error('The data must be an Nx2 matrix of points or an Nx3 matrix.')
end

%Create the delaunay triangulation object and the voronoi cell data.
tic
DT=DelaunayTri(Data);
fprintf('\nDone making the Delaunay triangulation.\n')
fprintf('Now starting the Voronoi diagram.  This may take a while.\n')
[VV, VC]=voronoiDiagram(DT);
fprintf('\nDone making Delaunay and Voronoi which took this long:\n')
toc
fprintf('\n')

%We need to remove any data point that has a voronoi cell that goes to
%infinity.  That is, we need to find any vertex that references a bad
%Voronoi vertex.  Infinity is bad.  Also any voronoi vertex outside of the
%data space.  Any finite voronoi vertex that is outside of the dataspace
%will not be inside any Delaunay tetra. delaunay tetra.

%find the indices of the bad Voronoi vertices.  pointLocation returns a nan
%if the voronoi vertex is NOT interior to any triangle/tetra.  We put a NaN
%in the first slot because the point at infinity is listed first in the
%Voronoi diagram.
ID = [NaN; pointLocation(DT,VV(2:length(VV),:))];
BadVV = find(isnan(ID));%creates index vector of bad vor verts

%Now we find the bad data.  "Bad data" is on the boundary of the data space
%so it has ind the vertices in the data that have a voronoi cell which
%includes a bad voronoi vertex.  Create an index of "BadDataID."
BadDataID = zeros(length(DT.X),1);
a=1;
for b=1:length(VC)
    if ~isempty(intersect(BadVV,VC{b}))
        BadDataID(a)=b;
        a=a+1;
    end
end
BadDataID(BadDataID==0)=[];

%Now we can use this BadDataID vector to remove any Tetras, Triangles and
%Edges that contain any of these bad points.  Then we will only be
%considering cells with GoodPoints for inclusion to the epsilon complex.

%get ALL edges (including bad ones)
DTedges=edges(DT);

%now remove bad edges.  They are "bad" if they contain a bad data point.
[R,~]=find(ismember(DTedges,BadDataID));
R=unique(R);
GoodEdges=DTedges;
GoodEdges(R,:)=[];%delete the bad edges.

%Here is the Good data if you want it. 
Good=Data;
Good(BadDataID,:)=[];

%Index of the good data if you want it.
GoodIndex=unique(GoodEdges);

fprintf('Finished the initial prep of the data.\n')
fprintf('Found the good data and the good edges.\n')

end%function


