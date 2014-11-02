function [hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT)

%%%%%%%%%%%%%%%%%% Max and Min class  structure section %%%%%%%%%%%%%%%%%%
%
% Info for maxclass struct; the analogous holds for minclass struct
%
% maxclass is an array of structs that is as long as the number of maximums
% in the data set.  So one entry for each max class, i.e. the "group" of
% vertices that are all a member of the same max class.  Here are the
% fields in the structure.  So for a max point, m, we want to know:
%   max - the index of the max point m.
%   points - the points in m's class
%   nbormax - the indices of the maximums that share a boundary with m.
%           This indexes the points themselves.
%   nbormaxid - the indices of the maxclasses that neighbor m's max class.
%   interiorv - the interior points of m's class
%   boundaryv - the boundary points
%   interiore - the interior edges (I want to improve this)
%   boundarye - boundary edges (I want to improve this)
%   bonds - the edges that attach to the boundary points in m's class to
%       neighboring max classes.  These are well defined.
%   hoptree - all of the edges in each of the hop paths in the max class.
%       There are no duplicate edges.  This allows plotting the hop tree.
maxclass=struct('max',{},'points',{},...
    'nbormax',{},'nbormaxid',{},'interiorv',{},'boundaryv',{},...
    'interiore',{},'boundarye',{},'bonds',{},'hoptree',{});
minclass=struct('min',{},'points',{},...
    'nbormin',{},'nborminid',{},'interiorv',{},'boundaryv',{},...
    'interiore',{},'boundarye',{},'bonds',{},'hoptree',{});

M=vertcat(hop.maxclass);
m=vertcat(hop.minclass);

% deal the max/min indexes to their respective classes
for i=1:length(maxindex)
    maxclass(i).max=maxindex(i);
    maxclass(i).points=find(M==maxindex(i));
end
for i=1:length(minindex)
    minclass(i).min=minindex(i);
    minclass(i).points=find(m==minindex(i));
end

%%%%%%%%%%Last major section %%%%%%
%
%   Now we have all of the information we need to start classifying every
%   point as a boundary point, interior point etc.  Also classify the edges
%   as best we can (more work to be done here regarding boundary edges.
%   Bonds are clear.


% assign to each point, the index of it's max/min class.  That way we can
% find out all of the information about a points class.
temp=0;
[hop(:).maxclassid]=deal(temp);
[hop(:).minclassid]=deal(temp);

temp = zeros(length(DT.X),1);

for i = 1:length(maxindex)
    temp(maxindex(i))=i; 
end

for i = 1:length(hop)
    if ~isnan(hop(i).maxclass)
        hop(i).maxclassid = temp(hop(i).maxclass);
    end
end

temp = zeros(length(DT.X),1);

for i = 1:length(minindex)
    temp(minindex(i))=i; 
end

for i = 1:length(hop)
    if ~isnan(hop(i).minclass)
        hop(i).minclassid = temp(hop(i).minclass);
    end
end

%%%% Outline for max classes:  same for minclasses
%Next, loop through hop.  For each good point, k:
%  pull the neighbors of k,
%  check to see if all the neighbors of k are in the same max class as k
%  if yes, then k is interior, catalog everything for k
% if not, then k is a boundary point
%     catalog the bond edges (can't do boundary edges yet)
%     pull the neighbors of k that are not in k's max class
%     find the max class for the neighbors of k; 
%           these max classes are the nbormax for the maxclass of k
%     

%deal false to hop.isboundaryM and m
[hop.ismaxboundary]=deal(false);
[hop.isminboundary]=deal(false);

%For Max Classes
for k=1:length(hop)
    if ~isempty(hop(k).edges)%then it is a good point
        Edges=hop(k).edges(:,[1 2]);%edges connected to k-th point
        %get the points in k's class
        ClassPoints = maxclass(hop(k).maxclassid).points;
        %find the nbors of k that are in k's class
        b=ismember(Edges(:,2),ClassPoints);            
        if b %all nbors of k are in the k-th point's max class
            hop(k).maxedgesin = Edges;
        else %k-th point is a boundary point find the bonds etc.
            %note, I can't find the boundary edges yet.  I will store the
            %bond edges for the k-th point; I will store the rest in
            %maxedgesin for now (even though some of them are not boundary
            %edges.  I will sort them out right after this block.
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
        else %k-th point is a boundary point find the bonds etc.
            %note, I can't find the boundary edges yet!  I will store the
            %bond edges for the k-th point; I will store the rest in
            %minedgesin for now (even though some of them are not boundary
            %edges.  I will sort them out right after this block.
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
        %neighbor in common.  That is, k and it's boundary neighbor must
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

%  Loop through the maxclass. Catalog everything.
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
    
% Loop through the minclass. Catalog everything.
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

%Now we catalog all of the hop edges in the maxclass structs.
%These are not ALL of the edges in a max class.  These are only the HOP
%edges for the hop paths for the points in the max class.

for a=1:length(maxclass)
%     NumEdges=size(horzcat(hop(maxclass(1).points).hopmaxpath),2)-...
%                                                 size(maxclass(1).points,1);
%     E=zeros(NumEdges,2); %preallocating the space for the edges.
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

%Now we do this for the min classes
for a=1:length(minclass)
%     NumEdges=size(horzcat(hop(minclass(1).points).hopminpath),2)-...
%                                                 size(minclass(1).points,1);
%     E=zeros(NumEdges,2); %preallocating the space for the edges.
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

end%function
