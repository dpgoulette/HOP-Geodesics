function [hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,...
   GoodIndex,DT)
%  HOPStructCreade - creates the main hop data structure and does the
%        geodesic HOP algorithm based on our density function:
%              
%                 1/(volume of the voronoi cell)
%
%        This function does HOP in the direction of greatest positive
%        gradient. Each point "hops" step-by-step until it reaches a max (a
%        point that has greater density than all of its neighbors).  We
%        store all results in the hop struct.  For each point in the data
%        set, with index p, we find the following data and store the
%        results in the following fields:
%
%              hop(p).density    <== the density of p.
%              hop(p).edges      <== the edges connected to p -- k by 2
%                                    array of indices of edge endpoints. 
%              hop(p).maxclass   <== the index of the max p "hops" to.
%              hop(p).minclass   <== the index of the min p "hops" to.
%              hop(p).ismax      <== true/false whether it is a max or not
%              hop(p).ismin      <== true/false whether it is a max or not
%              hop(p).hopmaxpath <== the path p takes to hop to its 
%                                    representative max.
%              hop(p).hopminpath <== the path p takes to hop to its 
%                                    representative min.
%
%     inputs:
%              VV - the voronoi vertices
%              VC - the indices of the voronoi vertices for each cell
%              GoodEdges - the set of edges we will "hop" on
%              GoodIndex - the index (into DT.X) of the "good points"
%              DT - The Delaunay triangulation object
%
%     outputs:
%              hop - the hop data structure updated with the maxclass
%                    assignment for each point and the hop path for each
%                    point.
%              maxindex - the indices (into DT.X) of the hop maxima (for
%                         plotting)
%              minindex - the indices (into DT.X) of the hop minima (for
%                         plotting)
%
%  hop is a vector of structs. Each entry in hop is a struct containing
%  the key information about that point (so the length of hop is as long as
%  the raw data DT.X). The key information that is stored for each point,
%  p, is in the following list.  (Note that we set the key for all of the
%  fields that we will want in hop for our work, but we do not fill in all
%  of the values of these fields in this function.  For example, the
%  maxclass field will be emtpy after this function terminates.  It will
%  not be filled in until we create the maxclass struct later.  See
%  HOPClasses.m and HOPGeadesics_main.m)
%       edges - the edges attached to p
%       density - the density of point p
%       maxclass - the index of the max point in p's max class (i.e. the
%                  max that p "hops" to)
%       minclass - the index of the min point in p's min class (i.e. the
%                  min that p "hops" to via reverse hop)
%       maxclassid - the index of the maxclass in the maxclass struct.  
%       minclassid - the index of the minclass in the minclass struct.
%       ismax - true or false if p is a max or not
%       ismin - true or false
%       ismaxboundary - true if p is a boundary point of it's max class
%       isminboundary - true if p is a boundary point of it's min class
%       maxedgesin - the edges connected to p that are interior edges to
%             p's max class (this has issues.. the code is fine the
%             definition might need some work)
%       minedgesin - same but for Min class
%       maxbonds - maxclass bonds attached to p.  Empty if p is interior
%              point.  If p is a boundary point then one of its neighbors
%              is in a different max class.  The edge that connects these
%              two points is a bond edge.
%       minbonds - Similar to maxbonds.
%       maxedgesboundary - Boundary edges for its max class (empty if it is
%               an interior point)
%       minedgesboundary - Similar to maxedgesboundary.
%       hopmaxpath - The path that p takes to hop to its max.  This is a
%               vector of point indices into DT.X.
%       hopminpath - The path that p takes to hop to its max.
%
%  COMMENTS: We are no longer using "interior edges" and "boundary edges."
%  Consider removing.
hop=struct('edges',{},'density',{},'maxclass',{},'minclass',{},...
   'maxclassid',{},'minclassid',{},...
   'ismax',{},'ismin',{},'ismaxboundary',{},'isminboundary',{},...
   'maxedgesin',{},'minedgesin',{},'maxbonds',{},'minbonds',{},...
   'maxedgesboundary',{},'minedgesboundary',{},'hopmaxpath',{},...
   'hopminpath',{});

% Double the edges for sorting.  We want every edge expressed both ways in
% the matrix for easy searching.
Edges=GoodEdges;
Edges=circshift(Edges,[0 1]);%flip the edges
Edges=[GoodEdges;Edges];
Edges=sortrows(Edges);% So E holds every edge expressed both ways.


%%%%%% Original Scalar function %%%%%%%
% Find the Voronoi cell volumes
DataVols=vor_cell_volumes(VV,VC);

%density values 1/(volume of voronoi).
Densities=1./DataVols;

% Assign to the third column of Edges the value of the density of the
% point in the second column.
Edges(:,3) = Densities(Edges(:,2));

% deal nans to the density slot to prallocate the struct length
[hop(1:size(DT.X,1)).density]=deal(NaN);

fprintf('\nStarting to sort the edges for HOP.\n')

% deal the edges connected to each vertex to the structure. They will be
% sorted later. We only catalog edges that are attached to "good" points,
% hence the need to use the "GoodIndex" vector.  Also, we sort the edges
% based on the density of the neighbor which is in the third column of the
% edges matrix.
for a=1:length(GoodIndex)
   r=Edges(:,1)==GoodIndex(a);
   hop(GoodIndex(a)).edges=sortrows(Edges(r,:),3);
   
   % Progress at the terminal.
   if mod(a,500) == 0 
      fprintf('.')
   end
end
fprintf(' Done.\n\n')

% Catalog the density of each good point in the hop struct.
for a = 1:length(GoodIndex)
   hop(GoodIndex(a)).density = Densities(GoodIndex(a));
end

% % Alternate option for the density function which we have currently
% % abandoned but we may revisit later.
% %
% % % User selects which density function they want to put on the data.
% % fprintf('Which density function would you like to use for the data set?\n')
% % while true
% %    fprintf('   1) 1/(Volume of voronoi).\n')
% %    fprintf('   2) (k + 1)/(Sum of the k+1 voronoi cells).\n')
% %    density_option = input('Choose one of the above: ');
% %    if density_option == 2 || density_option == 1
% %       break
% %    else
% %       fprintf('ERROR! You must enter 1 or 2.\n\n')
% %    end
% % end
% % 
% % % If density_option == 1 then this next block doesn't run.  So 1/Voronoi is
% % % used by default.
% % if density_option == 2
% %    % %%%%  ALTERNATE SCALAR FUNCTION  %%%%%%
% %    %
% %    % Comment out this block to use the original version: 1/vor This scalar
% %    % function depends on the values of 1/vor so it uses the results of the
% %    % code above.  THIS NEXT CODE BLOCK IS INEFFICIENT!!! IF WE DECIDE TO USE
% %    % THIS METHOD MAKE SURE TO REDO THIS SO THAT WE ARE NOT CALCULATING
% %    % DENSITIES TWICE!
% %    %
% %    % Alternate Scalar Function Description:
% %    % Each point has n neighbors.  Instead of 1/(volume of voronoi) we can do
% %    % the sum of the masses of the n neigbors plus the one point, devided by
% %    % the sum of the volumes of the n+1 voronoi cells.
% %    % This version uses (1+n)/(sum of the n+1 voronoi cells)
% %    %
% %    % This version of the density function required us to have the edges sorted
% %    % before this is done since it relies on the neighbor cells.
% %    
% %    for a = 1:length(hop)
% %       if ~isnan(hop(a).density)
% %          Neighbors = hop(a).edges(:,2);
% %          Neighbors_VorVols = DataVols(Neighbors);
% %          sum_vols = DataVols(a) + sum(Neighbors_VorVols);
% %          hop(a).density = (length(Neighbors_VorVols) + 1)/sum_vols;
% %       end
% %    end
% %    
% %    % Since we have recalculated the new densities using the old, we need to
% %    % fix the third column of the edges in the hop datastructure.
% %    % Asssign the density of the vertex in the second column to the third column
% %    for a = 1:length(hop)
% %       if ~isnan(hop(a).density)
% %          temp = hop(a).edges(:,2); % index of neighbors
% %          hop(a).edges(:,3) = vertcat(hop(temp).density);
% %          hop(a).edges = sortrows(hop(a).edges,3);
% %          %       for b = 1:size(hop(a).edges,1)
% %          %          hop(a).edges(b,3) = Densities(hop(GoodIndex(a)).edges(b,2));
% %          %       end
% %       end
% %    end
% % end % alternate density function section
% % Done calculating the density of each point.  Ready for HOP



% Now prepare for the HOP algorithm.  We are currently doing hop in the
% direction of greatest gradient.  (This is, of course, a discrete
% gradient.) The maxpointer and minpointer are n X 1 vectors contaning the
% HOP pointers. So the i-th entry in maxpointer contains the index of the
% point that the i-th point (in DT.X) hops to (it is the neigbor in the
% direction of greatest gradient.  Maximuma and minumuma point to Inf.  Bad
% points point to NaN.
[maxpointer, minpointer, hop] = GradientHop(hop,GoodIndex,DT);

fprintf('Done finding maxs, mins, and HOP path pointers.\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WARNING - this section won't work correctly unless the hop(:).maxclass
%%%% entries are empty to begin with.  So don't rerun this without clearing
%%%% those entries first!!
fprintf('Now finding max and min class assignment for every point.\n\n')

% Now that we know which points are maxima, we perform the HOP algorithm.
% We HOP from each point in the direction of greatest gradient until we hit
% a max.  Thus we find the max class for each point (i.e. we find the max
% that each point hops to).
for k=1:length(maxpointer)
   if isnan(maxpointer(k))
      % bad point so not in good dataset
      hop(k).maxclass=NaN;% Mark it as bad in hop
      % skip to the next.
      continue
   elseif ~isempty(hop(k).maxclass)
      % We already cataloged this point in a previous iteration if this for
      % loop.  Skip to the next.
      continue
   elseif isinf(maxpointer(k))
      % Point k is a max so there is nothing else to do.  It is the max in
      % its max class.
      hop(k).maxclass=k;
   else
      % k is a good point that is not a max and we havent cataloged it.  So
      % we need to run the HOP algorithm on this point until we hit a max
      % or we hit a point that we have already cataloged earlier in the
      % main for loop (thus we can just look up the rest of the path in the
      % hop struct).
      
      path=zeros(1,500); % Preallocate. Way too long.  That is fine.
      
      path(1)=k;%So the path starts on the current point.
      
      % traverse k-th path.  Stop if a path vertex points to inf.  This
      % vertex is the max. Every vertex in the path gets the index of this
      % max for its max class;  Also stop if a path vertex has non empty
      % max class value; we have already computed the rest of the path
      % earlier in the loop. Everything in the path gets the max class of
      % this final vertex.  Also, the path to the max is saved in the hop
      % data struct for each point.
      
      for a=1:length(maxpointer)
         if isinf(maxpointer(path(a))) 
            % then point a is a max.  The hop path is finished.
            path(a+1:end)=[];
            [hop(path).maxclass]=deal(path(a));
            for b=1:(length(path)-1)
               hop(path(b)).hopmaxpath=path(b:end);
            end
            break
            
         elseif ~isempty(hop(maxpointer(path(a))).maxclass)
            % Then we did the rest already.  So the max class for every
            % vertex in the path is the same as maxpointer(path(a))
            
            if hop(maxpointer(path(a))).ismax 
               % Thus path(a) points to a max.
               % Put the max in the last spot.
               path(a+1)=maxpointer(path(a));
               path(a+2:end)=[]; %delete the rest
               
               % Store the index of the maximum for this max class for
               % each point in the path.
               [hop(path(1:end-1)).maxclass]=...
                  deal(hop(maxpointer(path(a))).maxclass);
               
               for b=1:(length(path)-1)
                  hop(path(b)).hopmaxpath=path(b:end);
               end
            else
               path(a+1:end)=[];%delete the rest
               [hop(path(1:end)).maxclass]=...
                  deal(hop(maxpointer(path(a))).maxclass);
               
               % Now concatinate the current path with the path already
               % cataloged for the point that path(a) points to.  Then
               % store the paths in the uncataloged slots.
               path = [path, hop(maxpointer(path(end))).hopmaxpath];
               for b=1:a
                  hop(path(b)).hopmaxpath=path(b:end);
               end
            end
            break
         else %the path continues;
            path(a+1)=maxpointer(path(a));
         end
      end
      
   end
end

%find the min class for every point

for k=1:length(minpointer)
   if isnan(minpointer(k))%bad point so not in good dataset
      hop(k).minclass=NaN;%Mark it
      continue
   elseif ~isempty(hop(k).minclass)
      %We already cataloged this point in a previous iteration
      continue
   elseif isinf(minpointer(k))
      %point k is a min
      hop(k).minclass=k;
   else
      %k is not a min and we havent cataloged it
      path=zeros(1,1000);%way too long.  That is fine.
      
      path(1)=k;%So the path starts on the current point.
      
      %traverse k-th path.  Stop if a path vertex points to inf.  This
      %vertex is the min. Every vertex in the path gets the index of this
      %min for its min class;  Also stop if a path vertex has non empty
      %min class value; we have already computed the rest of the path
      %earlier in the loop. Everything in the path gets the min class of
      %this final vertex.  Also, the path to the min is saved in the hop
      %data struct for each point.
      
      for a=1:length(minpointer)
         if isinf(minpointer(path(a))) %then point a is a min
            path(a+1:end)=[];
            [hop(path).minclass]=deal(path(a));
            for b=1:(length(path)-1)
               hop(path(b)).hopminpath=path(b:end);
            end
            break
            
         elseif ~isempty(hop(minpointer(path(a))).minclass)
            %then we did the rest already.  So the min class for every
            %vertex in the path is the same as minpointer(path(a))
            
            if hop(minpointer(path(a))).ismin %thus path(a) points to a min
               path(a+1)=minpointer(path(a));%put the min in the last spot
               path(a+2:end)=[];%delete the rest
               
               %Store the index for the minimum for this min class for
               %each point in the path.
               [hop(path(1:end-1)).minclass]=...
                  deal(hop(minpointer(path(a))).minclass);
               for b=1:(length(path)-1)
                  hop(path(b)).hopminpath=path(b:end);
               end
            else
               path(a+1:end)=[];%delete the rest
               [hop(path(1:end)).minclass]=...
                  deal(hop(minpointer(path(a))).minclass);
               %Now concatinate the current path with the path already
               %cataloged for the point that path(a) points to.  Then
               %store the paths in the uncataloged slots.
               path = [path, hop(minpointer(path(end))).hopminpath];
               for b=1:a
                  hop(path(b)).hopminpath=path(b:end);
               end
            end
            break
         else %the path continues;
            path(a+1)=minpointer(path(a));
         end
      end
      
   end
end

fprintf('Done finding the max and min class assignments.\n\n')

% make an index array of the maxima and minima in the data set.  We will
% need this for plotting purposes.
maxindex=find(isinf(maxpointer));
minindex=find(isinf(minpointer));

end% main function
