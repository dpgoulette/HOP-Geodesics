function [hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,...
   GoodIndex,DT)
% Fix the comments for this function!!!
%
%
%  Issues 
%     1)  This code calculates 1/vor and then adjusts if the user wants the
%     other density function.  Inefficient.


%CREATE THE HOP DATA STRUCTURE
%  hop is an array of structs.  The length is as long as the raw data. The
%  information I want to know regarding any point, p, in the data set is
%  the following (in the order the are in the fields):
%       edges - the edges attached to p
%       density - the density of point p   i.e. 1/vorvolume
%       maxclass - the index of the max point in p's max class
%       minclass - the index of the min point in p's min class
%       maxclassid - the index of the maxclass in the maxclass struct
%       minclassid - the index of the minclass in the minclass struct
%       ismax - true or false if p is a max or not
%       ismin - true or false
%       ismaxboundary - true if p is a boundary point of it's max class
%       isminboundary - true if p is a boundary point of it's min class
%       maxedgesin - the edges connected to p that are interior edges to
%             p's max class (this has issues.. the code is fine the
%             definition might need some work)
%       minedgesin - same but for Min class
%       maxbonds - maxclass bonds attached to p.  Empty if p is interior
%              point.
%       minbonds - bonds for min class attached to p
%       maxedgesboundary - Boundary edges for its max class (empty if it is
%               an interior point)
%       minedgesboundary - Boundary edges for its min class (empty if it is
%               an interior point)
%
%  COMMENTS: We are no longer using "interior edges" and "boundary edges."
hop=struct('edges',{},'density',{},'maxclass',{},'minclass',{},...
   'maxclassid',{},'minclassid',{},...
   'ismax',{},'ismin',{},'ismaxboundary',{},'isminboundary',{},...
   'maxedgesin',{},'minedgesin',{},'maxbonds',{},'minbonds',{},...
   'maxedgesboundary',{},'minedgesboundary',{},'hopmaxpath',{},...
   'hopminpath',{});

% double the edges for sorting.  We want every edge expressed both ways in
% the matrix for easy searching.
Edges=GoodEdges;
Edges=circshift(Edges,[0 1]);%flip the edges
Edges=[GoodEdges;Edges];
Edges=sortrows(Edges);%So E holds every edge expressed both ways.


%%%%%% Original Scalar function %%%%%%%
% Find the Voronoi cell volumes
DataVols=vcellvolumes(VV,VC);

%density values 1/(volume of voronoi).
Densities=1./DataVols;

% Assign to the third column of Edges the value of the density of the
% point in the second column.
Edges(:,3) = Densities(Edges(:,2));

% deal nans to the density slot to prallocate the struct length
[hop(1:size(DT.X,1)).density]=deal(NaN);

tic
fprintf('Starting to sort the edges for HOP.  This may take a while.\n')

% deal the edges connected to each vertex to the structure. They will be
% sorted later. We only catalog edges that are attached to "good" points,
% hence the need to use the "GoodIndex" vector.  Also, we sort the edges
% based on the density of the neighbor which is in the third column of the
% edges matrix.
for a=1:length(GoodIndex)
   r=Edges(:,1)==GoodIndex(a);
   hop(GoodIndex(a)).edges=sortrows(Edges(r,:),3);
   %fprintf('%d out of %d.\n',a,length(GoodIndex)-1)
end
fprintf('Finished sorting the hop edges into the hop data struct.\n')
toc
fprintf('\n')

% Catalog the density of each good point
for a = 1:length(GoodIndex)
   hop(GoodIndex(a)).density = Densities(GoodIndex(a));
end

fprintf('The densities are now cataloged in the hop data structure.\n')
toc

fprintf('\n')

% User selects which density function they want to put on the data.
fprintf('Which density function would you like to use for the data set?\n')
while true
   fprintf('   1) 1/(Volume of voronoi).\n')
   fprintf('   2) (k + 1)/(Sum of the k+1 voronoi cells).\n')
   density_option = input('Choose one of the above: ');
   if density_option == 2 || density_option == 1
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end

% If density_option == 1 then this next block doesn't run.  So 1/Voronoi is
% used by default.
if density_option == 2
   % %%%%  ALTERNATE SCALAR FUNCTION  %%%%%%
   %
   % Comment out this block to use the original version: 1/vor This scalar
   % function depends on the values of 1/vor so it uses the results of the
   % code above.  THIS NEXT CODE BLOCK IS INEFFICIENT!!! IF WE DECIDE TO USE
   % THIS METHOD MAKE SURE TO REDO THIS SO THAT WE ARE NOT CALCULATING
   % DENSITIES TWICE!
   %
   % Alternate Scalar Function Description:
   % Each point has n neighbors.  Instead of 1/(volume of voronoi) we can do
   % the sum of the masses of the n neigbors plus the one point, devided by
   % the sum of the volumes of the n+1 voronoi cells.
   % This version uses (1+n)/(sum of the n+1 voronoi cells)
   %
   % This version of the density function required us to have the edges sorted
   % before this is done since it relies on the neighbor cells.
   
   for a = 1:length(hop)
      if ~isnan(hop(a).density)
         Neighbors = hop(a).edges(:,2);
         Neighbors_VorVols = DataVols(Neighbors);
         sum_vols = DataVols(a) + sum(Neighbors_VorVols);
         hop(a).density = (length(Neighbors_VorVols) + 1)/sum_vols;
      end
   end
   
   % Since we have recalculated the new densities using the old, we need to
   % fix the third column of the edges in the hop datastructure.
   % Asssign the density of the vertex in the second column to the third column
   for a = 1:length(hop)
      if ~isnan(hop(a).density)
         temp = hop(a).edges(:,2); % index of neighbors
         hop(a).edges(:,3) = vertcat(hop(temp).density);
         hop(a).edges = sortrows(hop(a).edges,3);
         %       for b = 1:size(hop(a).edges,1)
         %          hop(a).edges(b,3) = Densities(hop(GoodIndex(a)).edges(b,2));
         %       end
      end
   end
end % alternate density function section
% Done calculating the density of each point.  Ready for HOP

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               HOP SECTION
% %   THIS SECTION HAS BEEN REMOVED TO A SEPARATE FUNCTION!!!
% %   By separating it, we can alter HOP more easily.
% %
% %   This version of HOP hops to the highest density neighbor. This
% %   was our original hop approach.
% %
% % maxpointer=NaN(length(DT.X),1);
% % minpointer=NaN(length(DT.X),1);
% %
% % %find maximums; A max has density larger than any of it's neigbors. mark in
% % %hop struct which points are maxs; create pointer list for the HOP paths;
% % %there will be a NaN at an entry for a bad point. There will be an Inf if
% % %the point is a max.
% % for z=1:length(GoodIndex)
% %     if hop(GoodIndex(z)).edges(end,3)<hop(GoodIndex(z)).density
% %         %GoodIndex(z) is a max
% %         hop(GoodIndex(z)).ismax = true;
% %         maxpointer(GoodIndex(z))=Inf;
% %     else%GoodIndex(z) is not a max
% %         maxpointer(GoodIndex(z))=hop(GoodIndex(z)).edges(end,2);
% %     end
% % end
% %
% % %find minimums;mark in hop struct;create pointer list; there will be a NaN
% % %at an entry for a bad point.  There will be an Inf if the point is a min.
% % for z=1:length(GoodIndex)
% %     if hop(GoodIndex(z)).edges(1,3)>hop(GoodIndex(z)).density
% %         %GoodIndex(z) is a min
% %         hop(GoodIndex(z)).ismin = true;
% %         minpointer(GoodIndex(z))=Inf;
% %     else%GoodIndex(z) is not a max
% %         minpointer(GoodIndex(z))=hop(GoodIndex(z)).edges(1,2);
% %     end
% % end
% % %%%%%%%%%%%%%%%%%%%%%%%  END OF HOP SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%
% % clear z
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%This function replaces the removed block above.  The maxpointer and
%minpointer are n X 1 vectors contaning the direction HOP hops from each
%point.  Maximums and minumums point to Inf.  Bad points point to NaN.
[maxpointer, minpointer, hop] = GradientHop(hop,GoodIndex,DT);


fprintf('Done finding maxs, mins, and HOP path pointers.\n')
toc
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WARNING - this section won't work correctly unless the hop(:).maxclass
%%%% entries are empty to begin with.  So don't rerun this without clearing
%%%% those entries first!!
fprintf('Now finding max and min class assignment for every point.\n')

%find the max class for each point
tic
for k=1:length(maxpointer)
   if isnan(maxpointer(k))%bad point so not in good dataset
      hop(k).maxclass=NaN;%Mark it
      continue
   elseif ~isempty(hop(k).maxclass)
      %We already cataloged this point in a previous iteration
      continue
   elseif isinf(maxpointer(k))
      %point k is a max
      hop(k).maxclass=k;
   else
      %k is not a max and we havent cataloged it
      path=zeros(1,100);%way too long.  That is fine.
      
      path(1)=k;%So the path starts on the current point.
      
      %traverse k-th path.  Stop if a path vertex points to inf.  This
      %vertex is the max. Every vertex in the path gets the index of this
      %max for its max class;  Also stop if a path vertex has non empty
      %max class value; we have already computed the rest of the path
      %earlier in the loop. Everything in the path gets the max class of
      %this final vertex.  Also, the path to the max is saved in the hop
      %data struct for each point.
      
      for a=1:length(maxpointer)
         if isinf(maxpointer(path(a))) %then point a is a max
            path(a+1:end)=[];
            [hop(path).maxclass]=deal(path(a));
            for b=1:(length(path)-1)
               hop(path(b)).hopmaxpath=path(b:end);
            end
            break
            
         elseif ~isempty(hop(maxpointer(path(a))).maxclass)
            %then we did the rest already.  So the max class for every
            %vertex in the path is the same as maxpointer(path(a))
            
            if hop(maxpointer(path(a))).ismax %thus path(a) points to a max
               path(a+1)=maxpointer(path(a));%put the max in the last spot
               path(a+2:end)=[];%delete the rest
               
               %Store the index for the maximum for this max class for
               %each point in the path.
               [hop(path(1:end-1)).maxclass]=...
                  deal(hop(maxpointer(path(a))).maxclass);
               for b=1:(length(path)-1)
                  hop(path(b)).hopmaxpath=path(b:end);
               end
            else
               path(a+1:end)=[];%delete the rest
               [hop(path(1:end)).maxclass]=...
                  deal(hop(maxpointer(path(a))).maxclass);
               %Now concatinate the current path with the path already
               %cataloged for the point that path(a) points to.  Then
               %store the paths in the uncataloged slots.
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
toc
fprintf('\n')

tic
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

fprintf('Done finding the max and min class assignments.\n')
toc
fprintf('\n')

%make an index of the maximums and mins in the data set
maxindex=find(isinf(maxpointer));
minindex=find(isinf(minpointer));

end%function
