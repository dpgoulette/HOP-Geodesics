function [hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT)
% HOPClasses - Creates the maxclass and minclass structs which hold the
%     data for each maxclass and minclass.  A maxclass is the set of all
%     points that hop to the same max (similar for minclass).  The max and
%     min class for each good point in the data set is calculated and a
%     variety of related information is stored in the maxclass, minclass
%     and hop structs (details below).  The key idea is that the maxclasses
%     (and similarly minclasses) form a partition of the data space (so the
%     maxclasses are equivalence classes).
%
%     input:
%        hop - the hop data struct
%        maxindex - the indices of the maxima (index into DT.X)
%        minindex - the indices of the minima (index into DT.X)
%        DT - the delaunay triangulation object
%
%     output:
%        hop - updated with maxclass and minclass information
%        maxclass - the maxclass struct (explained below)
%        minclass - the minclass struct (explained below)
%
%
%     %%%%%%%%%%%%% Explanation of maxclass struct %%%%%%%%%%%%%%
%
%     We explain the maxclass struct here.  The minclass struct is similar.
%
%     "maxclass" is a vector of structs that is as long as the number of
%     maxima in the data set.  So each entry represents a max class.  A max
%     class is the collection of data points (or vertices) that all hop to
%     the same max (so we say they are all a member of the same max class
%     since this is an equivalence relation that partitions the data). Here
%     are the fields in the maxclass structure. So for a max point, m, we
%     want to know:
%        max -    the index of the max point m.
%        points - the points in m's class
%        nbormax - the indices of the maxima that share a boundary with m.
%                  This is an index of the neighboring maxima in the data
%                  matrix DT.X.
%        nbormaxid - the indices of the maxclasses that neighbor m's max class.
%                    This is an index of the neighboring maxima in the
%                    maxclass struct.
%        interiorv - the interior points of m's class
%        boundaryv - the boundary points
%        interiore - the interior edges (I want to improve this)
%        boundarye - boundary edges (I want to improve this)
%        bonds - the edges that attach to the boundary points in m's class to
%                neighboring max classes.  These are well defined.
%        hoptree - all of the edges in each of the hop paths in the max class.
%                  There are no duplicate edges.  This allows plotting the hop
%                  tree.
%        geodesics - the geodesic paths that connect neighboring maxima.
%        geo_tris - the geodesic triangles
%        geo_tetras - the geodesic tetras
%
%     NOTE! The following fields are created in this function but they will
%     be empty when this function terminates: 
%           geodesics
%           geo_tris 
%           geo_tetras  
%     These fields will be filled in later by the following functions:
%           hopmaxconnect
%           SelectGeodesics
%           GeoTris
%           GeoTetras
%           geodesic_classify
%     See HOPGeodesics_main as well as these functions for more
%     information.
%
%     
%     %%%%%%%%%% Explanation of data stored in the hop struct %%%%%%%%%%%
%
%     Many of the fields in the hop struct are updated in this function.
%     For a give point with index p, the following describes which fields
%     are updated and what is assigned to them:
%
%        hop(p).maxclassid          <==  The index of the max class that p
%                                        belongs to.  This is an index into
%                                        the maxclass struct.
%        hop(p).minclassid          <==  The index of the min class that p
%                                        belongs to.  This is an index into
%                                        the minclass struct.
%        hop(p).ismaxboundary       <==  true/false if p is a boundary
%                                        point of its maxclass.
%        hop(p).isminboundary       <==  true/false if p is a boundary
%                                        point of its minclass.
%        hop(p).maxedgesin          <==  Edges attached to p that are
%                                        interior to its maxclass.
%        hop(p).minedgesin          <==  Edges attached to p that are
%                                        interior to its minclass.
%        hop(p).maxbonds            <==  Edges attached to p such that the
%                                        point at the other end of the edge
%                                        is in a different (neighboring)
%                                        maxclass.  These are the bond
%                                        edges.
%        hop(p).minbonds            <==  Edges attached to p such that the
%                                        point at the other end of the edge
%                                        is in a different (neighboring)
%                                        minclass.  These are the bond
%                                        edges.
%        hop(p).maxedgesboundary    <==  Edges attached to p such that the
%                                        point at the other end of the edge
%                                        is in a boundary point of the same
%                                        maxclass.  (Not fully functional
%                                        currently.)
%        hop(p).minedgesboundary    <==  Edges attached to p such that the
%                                        point at the other end of the edge
%                                        is in a boundary point of the same
%                                        maxclass.  (Not fully functional
%                                        currently.)



% Initialize the entries in the maxclass and minclass structs with all
% fields.
maxclass_temp=struct('max',[],'points',[],...
    'nbormax',[],'nbormaxid',[],'interiorv',[],'boundaryv',[],...
    'interiore',[],'boundarye',[],'bonds',[],'hoptree',[],'geodesics',[],...
    'geo_tris',[], 'geo_tetras',[]);
minclass_temp=struct('min',[],'points',[],...
    'nbormin',[],'nborminid',[],'interiorv',[],'boundaryv',[],...
    'interiore',[],'boundarye',[],'bonds',[],'hoptree',[]);

% preallocate the maxclass and minclass structs
maxclass = repmat(maxclass_temp,1,size(maxindex,1));
minclass = repmat(minclass_temp,1,size(minindex,1));

% Get the indices of all the maxima and minima
M=vertcat(hop.maxclass);
m=vertcat(hop.minclass);

% Catalog the max/min point indices in their respective classes.  Also find
% all points that hop to max i and store those points in maxclass(i).points
for i=1:length(maxindex)
    maxclass(i).max=maxindex(i);
    maxclass(i).points=find(M==maxindex(i));
end
for i=1:length(minindex)
    minclass(i).min=minindex(i);
    minclass(i).points=find(m==minindex(i));
end

% Assign to each point in hop the index of it's max/min class (not the
% index of the max itself).  Then hop(i).maxclassid will hold the index into
% the maxclass struct that represents the maxclass for i.  So we can always
% retrieve the struct that holds i's maxclass by doing this:
%     maxclass(hop(i).maxclassid)

for i = 1:length(maxclass)
   hop(maxclass(i).max).maxclassid = i;
end

for i = 1:length(hop)
   % Check that the point is not a bad point (it will have NaN for density)
   % Also check that it is not a max (we have already done the maxima above).
   if ~isnan(hop(i).density) && ~hop(i).ismax
      % Get i's max.  They have the same maxclassid.  Assign it.
      max_of_i = hop(i).maxclass;
      hop(i).maxclassid = hop(max_of_i).maxclassid;
   end
end

% Do the same for mins
for i = 1:length(minclass)
   hop(minclass(i).min).minclassid = i;
end

for i = 1:length(hop)
   if ~isnan(hop(i).density) && ~hop(i).ismin
      % Get i's min.  They have the same minclassid.  Assign it.
      min_of_i = hop(i).minclass;
      hop(i).minclassid = hop(min_of_i).minclassid;
   end
end


%%%%%%%%%%% done assigning max and min classes %%%%%%%%%%%%%%%%


%%%% Outline for the next code section
% Explanation for max classes. minclass is similar.
% 
% loop through hop.  For each good point, k:
%   pull the neighbors of k,
%   check to see if all the neighbors of k are in the same max class as k
%   if yes
%      then k is interior to its maxclass, catalog everything for k
%   else: at least one of its neigbors is in another maxclass
%      then k is a boundary point
%      catalog the bond edges (which are the edges connecting two points in
%                              DIFFERENT max classes)
%      pull the neighbors of k that are not in k's max class
%      find the max class for the neighbors of k; 
%           these max classes are the nbormax for the maxclass of k


% For Max Classes
for k=1:length(hop)
    if ~isempty(hop(k).edges)%then it is a good point
        Edges=hop(k).edges(:,[1 2]);%edges connected to k-th point
        %get the points in k's class
        ClassPoints = maxclass(hop(k).maxclassid).points;
        %find the nbors of k that are in k's class
        b=ismember(Edges(:,2),ClassPoints);            
        if b %all nbors of k are in the k-th point's max class
            hop(k).maxedgesin = Edges;
        else
            % k-th point is a boundary point find the bonds etc. Store the
            % bond edges for the k-th point; I will store the rest in
            % maxedgesin.
            hop(k).ismaxboundary = true;
            hop(k).maxedgesin = Edges(b,:);
            hop(k).maxbonds = Edges(~b,:);
        end
    end
end

%for Min Classes
for k=1:length(hop)
    if ~isempty(hop(k).edges)%then it is a good point
        Edges=hop(k).edges(:,[1 2]);%edges connected to k-th point
        %get the points in k's class
        ClassPoints = minclass(hop(k).minclassid).points;
        %find the nbors of k that are in k's class
        b=ismember(Edges(:,2),ClassPoints);            
        if b 
            %Thus, all nbors of k are in the k-th point's min class.
            %Therefore, k is an interior point            
            hop(k).minedgesin = Edges;
        else 
            hop(k).isminboundary = true;
            hop(k).minedgesin = Edges(b,:);
            hop(k).minbonds = Edges(~b,:);
        end
    end
end

%%%%%%
%  Now find the boundary edges (which we can now do since we have cataloged
%  every point as a boundary point or interior point to each max/min class
%  loop through hop; pull all the maxedgesin (which are both the potential
%  interior edges connected to k and any boundary edges as well. Check
%  whether each neighbor data point is a boundary point or not
%  if an InEdge neighbor is a boundary point then then the edge is a
%           boundary edge
%  else the edge is an interior edge so leave it in maxedgesin
%
% Now we sort the interior edges from the boundary edges for the maxclasses

for k=1:length(hop)
    %we only have to separate edges if k is a good point AND k is a
    %boundary point (if it is an interior point then nothing to do)
    if ~isempty(hop(k).edges) && hop(k).ismaxboundary
        %To be a boundary edge, the neighbor of k needs to be a boundary
        %point.  Also, the k-th point and the neighbor must have a bond
        %neighbor in common.  That is, k and it's boundary neighbor must be
        %adjacent to the same vertex in another max class.
        %So if these are both true, then the neighbor and k are two
        %vertices in a triangle where the third vertex is in another max
        %class.  Thus, the edge between k and it's neighbor is a boundary
        %edge.     
        
        b=false(size(hop(k).maxedgesin,1),1);
        
        %There is one strange case.  If a max is isolated then it's max
        %class has no interior edges so we need to check if b is empty. In
        %any other case, b cannot be empty.
        if ~isempty(b)
            A=hop(k).maxbonds(:,2);%bond neighbors of k
            for i=1:size(b,1)
                if hop(hop(k).maxedgesin(i,2)).ismaxboundary
                    % then the in-class-nbor is a boundary vert
                    % bond neighbors of k's nbor
                    B=hop(hop(k).maxedgesin(i,2)).maxbonds(:,2);
                    if ~isempty(intersect(A,B))
                        %then the in-nbor shares a bond-nbor so it is a
                        %boundary edge.  So mark it in b.
                        b(i)=true;
                    end
                end
            end
        end
        %assign the boundary edges and delete them from the interior
        %set
        hop(k).maxedgesboundary=hop(k).maxedgesin(b,:);
        hop(k).maxedgesin(b,:)=[];
    end
end

% %Now we sort the interior edges from the boundary edges for the minclasses

for k=1:length(hop)
    %we only have to separate edges if k is a good point AND k is a
    %boundary point (if it is an interior point then nothing to do)
    if ~isempty(hop(k).edges) && hop(k).isminboundary
        %To be a boundary edge, the neighbor of k needs to be a boundary
        %point.  Also, the k-th point and the neighbor must have a bond
        %neighbor in common.  That is, k and it's boundary neighbor must
        %adjacent to the same vertex in another min class.
        %So if these are both true, then the neighbor and k are two
        %vertices in a triangle where the third vertex is in another min
        %class.  Thus, the edge between k and it's neighbor is a boundary
        %edge.     
        
        b=false(size(hop(k).minedgesin,1),1);
        
        %There is one strange case.  If a min is isolated then it's min
        %class has no interior edges so we need to check if b is empty. In
        %any other case, b cannot be empty.
        if ~isempty(b)
            A=hop(k).minbonds(:,2);%bond neighbors of k
            for i=1:size(b,1)
                if hop(hop(k).minedgesin(i,2)).isminboundary
                    %then the in-class-nbor is a boundary vert
                    B=hop(hop(k).minedgesin(i,2)).minbonds(:,2);%bond neighbors of k's nbor
                    if ~isempty(intersect(A,B))
                        %then the in-nbor shares a bond-nbor so it is a
                        %boundary edge.  So mark it in b.
                        b(i)=true;
                    end
                end
            end
        end
        %assign the boundary edges and delete them from the interior
        %set
        hop(k).minedgesboundary=hop(k).minedgesin(b,:);
        hop(k).minedgesin(b,:)=[];
    end
end

%  Loop through the maxclass. Catalog the following:
%     maxclass.boundaryv
%     maxclass.interiorv
%     maxclass.boundarye
%     maxclass.interiore
%     maxclass.bonds
%     maxclass.nbormax
%     maxclass.nbormaxid
for k=1:length(maxclass)
    P=maxclass(k).points;
    b=vertcat(hop(P).ismaxboundary);
    maxclass(k).boundaryv = P(b);
    maxclass(k).interiorv = P(~b);
    
    %boundary edges.  There is redundancy because each boundary edge will
    %be cataloged twice, once for each end of the edge.  The set of all
    %boundary edges for the maxclass can't be empty.
    Bedge = vertcat(hop(P).maxedgesboundary);
    if ~isempty(Bedge)
        Bedge = sort(Bedge,2);
        maxclass(k).boundarye = unique(Bedge,'rows');
    end
    %interior edges.  There is redundancy.  There might not be any inter
    Iedge = vertcat(hop(P).maxedgesin);
    if ~isempty(Iedge)
        Iedge = sort(Iedge,2);
        maxclass(k).interiore = unique(Iedge,'rows');
    end
    
    %bonds must be unique.
    maxclass(k).bonds = vertcat(hop(P).maxbonds);
    
    %find the neighboring maximums for each max.  That way we can hop on
    %the maximums. Or find neghboring max class info.
    if ~isempty(maxclass(k).bonds)
        BN = maxclass(k).bonds(:,2);
        NborMaxs = vertcat(hop(BN).maxclass);
        maxclass(k).nbormax = unique(NborMaxs);
        maxclass(k).nbormaxid=vertcat(hop(maxclass(k).nbormax).maxclassid);
    end
end
    
%  Loop through the minclass. Catalog the following:
%     minclass.boundaryv
%     minclass.interiorv
%     minclass.boundarye
%     minclass.interiore
%     minclass.bonds
%     minclass.nbormin
%     minclass.nborminid
for k=1:length(minclass)
    P=minclass(k).points;
    b=vertcat(hop(P).isminboundary);
    minclass(k).boundaryv = P(b);
    minclass(k).interiorv = P(~b);
    
    %boundary edges.  There is redundancy.
    Bedge = vertcat(hop(P).minedgesboundary);
    Bedge = sort(Bedge,2);
    minclass(k).boundarye = unique(Bedge,'rows');
    
    %interior edges.  There is redundancy.
    Iedge = vertcat(hop(P).minedgesin);
    if ~isempty(Iedge)
        Iedge = sort(Iedge,2);
        minclass(k).interiore = unique(Iedge,'rows');
    end
    
    %bonds must be unique.
    minclass(k).bonds = vertcat(hop(P).minbonds);
    
    %find the neighboring minimums
    if ~isempty(minclass(k).bonds)
        BN = minclass(k).bonds(:,2);
        NborMins = vertcat(hop(BN).minclass);
        minclass(k).nbormin = unique(NborMins);
        minclass(k).nborminid=vertcat(hop(minclass(k).nbormin).minclassid);
    end
end

% Now we catalog all of the hop path edges in the maxclass struct. These
% are not ALL of the edges in a max class.  These are only the HOP path
% edges in all of the hop paths for every point in the max class.  The set
% of hop edges in a maxclass creates a tree with the max as the root.  The
% entire data space is then a forest of trees where each connected
% component is a maxclass. Cataloging these edges in one place allows us to
% plot the hop trees for visualisation purposes.
for a=1:length(maxclass)
    EdgesCell = cell(length(maxclass(a).points),1);
    for b=1:length(EdgesCell)
        Path = hop(maxclass(a).points(b)).hopmaxpath;
        if length(Path)<=2
            % there is only one edge in the path
            EdgesCell{b} = Path;
        else % length(Path) > 2
            % then there is more than one edge in the path
            EdgesCell{b} = [Path(1:end-1)', Path(2:end)'];
        end
    end
    E = vertcat(EdgesCell{:});
    E = sort(E,2);
    E = unique(E,'rows');
    maxclass(a).hoptree = E;
end

%Now we do the same for the min classes
for a=1:length(minclass)
    EdgesCell = cell(length(minclass(a).points),1);
    for b=1:length(EdgesCell)
        Path = hop(minclass(a).points(b)).hopminpath;
        if length(Path)<=2
            % there is only one edge in the path
            EdgesCell{b} = Path;
        else % length(Path) > 2
            % then there is more than one edge in the path
            EdgesCell{b} = [Path(1:end-1)', Path(2:end)'];
        end
    end
    E = vertcat(EdgesCell{:});
    E = sort(E,2);
    E = unique(E,'rows');
    minclass(a).hoptree = E;
end

end% main function: HOPClasses
