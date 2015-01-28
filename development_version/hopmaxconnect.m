function maxclass = hopmaxconnect(points,maxclass,hop)
%  FIX - if two neighbor max paths have the same max-length edge, then we
%     need to check the second longest, and then the third etc, until the tie
%     is broken.  Currently the code only checks the longest edge.
%
%  FIX? Conider refactoring code to avoid using maxconnect
%
%  FIX - update the comment/help section
%
%  FIX? - (I think this is actually good as is.  Probably don't need to fix
%     it.  But there may be a way to speed it up if needed.  Might be slow
%     on the full SDSS.) Speed up cataloging the reverse geodesic
%     direction.  Avoid the "tempedges" matrix and the "ismember" call.
%
%
%  COMPLETED THIS FIX (10/13) ==>  THIS FUNCTION CAN BE IMPROVED FOR SPEED.
%     Instead of calculating the hop path, just pull it from
%     hop(a).hopmaxpath.
%
% hopmaxconnectV2 - finds the geodesic path which connects two neighboring
%       maxima.  Two maxima are neigbors if they share a bond edge, which
%       means that they have points on their boundary that are adjacent in
%       the Delaunay graph (1-skeleton).  The geodesic is the path which
%       goes through a common bond (between the two neigboring max classes)
%       and minimizes the length of the longest edge. The geodesics between
%       each pair of neigboring maxima are stored in the maxconnect data
%       struct.

%this counts the number of geodesics we will need to store so we can
%preallocate the maxconnect cell.
count=0;
for a=1:length(maxclass)
   count = count+size(maxclass(a).nbormax,1);
end
maxconnect = cell(count,2);
[maxconnect{:,1}]=deal([0,0,0]);
tempedges = zeros(count,2);
k=1;

for a=1:length(maxclass)
   %the next if statement guards against the case where a max class is
   %isolated (rare artifact on the boundary of the data).
   if ~isempty(maxclass(a).bonds)
      % Put a placeholder in maxclass(a).geodesics where the geodesics will
      % be stored.  There are four columns. 
      %  Column 1: the indices of the points in the geodesic.
      %  Column 2: the length of the geodesic
      %  Column 3: the rank of the geodesic
      %  Column 4: whether the geodesic was chosen for inclusion
      %         (calculated later in 'SelectGeodesics.m' see HOPGeodesics_main
      maxclass(a).geodesics = cell(size(maxclass(a).nbormaxid,1),4);
      B=maxclass(a).bonds;%B holds the bonds of maxclass a
      for b=1:size(maxclass(a).nbormaxid,1)
         if maxclass(a).nbormaxid(b) < a
            %thus we have already cataloged this geodesic in a previous
            %iteration.  So we simply need to flip it around and save it in
            %the current max class iteration.
            r = ismember(tempedges,[maxclass(maxclass(a).nbormaxid(b)).max,...
               maxclass(a).max],'rows');
            % Store the results in maxconnect and maxclass(a).geodesics
            maxconnect{k,1}(3)=maxconnect{r,1}(3);
            maxconnect{k,2}=fliplr(maxconnect{r,2});
            maxconnect{k,1}([1 2])= fliplr(maxconnect{r,1}([1 2]));
            maxclass(a).geodesics{b,1} = maxconnect{k,2};
            maxclass(a).geodesics{b,2} = maxconnect{k,1};
            tempedges(k,:) = maxconnect{k,1}([1 2]);
            k=k+1;
         else
            %Bn will hold the bonds of maxclass b, where maxclass b is
            %a neighbor of maxclass a.
            Bn=maxclass(maxclass(a).nbormaxid(b)).bonds;
            Bn=fliplr(Bn);
            r=ismember(B,Bn,'rows');%find a and b's common bond edges
            common_bonds=B(r,:);%temp holds the common bonds.
            Paths = cell(size(common_bonds,1),1);
            Plengths = zeros(size(common_bonds,1),1);
            
            % TEST VERSION FOR TIE BREAKING
            % This cell will store the lengths of all the edges in all the
            % max paths.  We need them in case there is a tie for the
            % shortest maximum edge length.
            Plengths_all = cell(size(common_bonds,1),1);
            
            for c=1:size(common_bonds,1)
              
               P1 = hop(common_bonds(c,1)).hopmaxpath;
               P2 = hop(common_bonds(c,2)).hopmaxpath;
               if isempty(P1)
                  % then temp(c,1) is a boundary point that is also a
                  % max! So it's hop path has no edges.  Thus we must
                  % tack it on to the geodesic by making P1 contain
                  % this point only.  i.e. P1 is a single vertex "path."
                  P1 = common_bonds(c,1);
               end
               if isempty(P2)
                  %same as for P1 above.
                  P2 = common_bonds(c,2);
               end
               
               Paths{c}=[fliplr(P1), P2];
               % now find the euclidean length of the longest edge in
               % the path from max to max.  Also return the lengths of all
               % the edges 
               [Plengths(c), Plengths_all{c}] =Plength(Paths{c},points);
               
            end
            % Now remove the bond paths we just calculated from the B
            % matrix for the next iteration of this for loop.
            B(r,:)=[];
            [min_length, id] = min(Plengths);   
            
            % Even if we have to go to a tie breaker for the geodesic, we
            % already have the length of it, so store it now.
            maxconnect{k,1}(3) = min_length;
            
            slots_with_min = find(Plengths == min_length);
            
            if length(slots_with_min) == 1
               % Then there is only one geodesic that attains the minimum
               % max-edge length, so store the geodesic in the maxconnect and
               % maxclass structures
               maxconnect{k,2}=Paths{id};
               maxconnect{k,1}([1 2])= [maxconnect{k,2}(1),maxconnect{k,2}(end)];
               maxclass(a).geodesics{b,1} = Paths{id};
               maxclass(a).geodesics{b,2} = maxconnect{k,1};
            else
               % There are 2 or more paths that tie for the minimum max
               % edge length. We need to check the second longest edge,
               % then the third longest etc. until we break the tie.
               
               % This index means to compare the 2nd longest edge.  Then
               % index if there is still a tie, etc.
               edge_index = 2; 
               while length(slots_with_min) > 1                  
                  Temp = zeros(size(slots_with_min));
                  for P = 1:length(slots_with_min)
                     Temp(P) = Plengths_all{slots_with_min(P)}(edge_index);
                  end
                  [tie_breaker, id] = min(Temp);
                  Temp2 = find(Temp == tie_breaker);
                  
                  edge_index = edge_index + 1;
                  r = Temp ~= tie_breaker;
                  slots_with_min(r) = [];
               end
%                maxconnect{k,1}(3) = min_length;
%                [maxconnect{k,1}(3), id]=min(Plengths);%minimum max edge length path.
               
               % Store the geodesic in maxconnect and maxclass
               maxconnect{k,2}=Paths{slots_with_min};
               maxconnect{k,1}([1 2])= [maxconnect{k,2}(1),maxconnect{k,2}(end)];
               maxclass(a).geodesics{b,1} = Paths{slots_with_min};
               maxclass(a).geodesics{b,2} = maxconnect{k,1};
            end
            
            %the next line is to catalog which pairs of neighboring maxima
            %we have already calculated.  We need this to catalog it in the
            %other direction efficiently in the first if statement after
            %the b index for loop.
            tempedges(k,:) = [maxconnect{k,2}(1),maxconnect{k,2}(end)];
            k=k+1;
         end
      end
   end
   if mod(a,100)==0
      fprintf('%d out of %d completed.\n',a,length(maxclass))
   end
end

%%%%%%%%% Sort the maxconnect cell contents %%%%%%%%%%
% Sort the geodesics connected to each max; so when we are considering the
% kth max in maxclass, the geodesics connected to k are in order from
% shortest to longest.  Also sort maxclass(:).geodesics at the same time.
GeoIndex = 1;
for i = 1:length(maxclass)
   if ~isempty(maxclass(i).nbormax)
      NumMaxNeighbors = length(maxclass(i).nbormaxid);
      common_bonds = vertcat(maxconnect{GeoIndex:GeoIndex+NumMaxNeighbors-1,1});
      [~,sortID] = sort(common_bonds(:,3));
      maxclass(i).geodesics = maxclass(i).geodesics(sortID,:);
      sortID = sortID+(GeoIndex-1);
      maxconnect(GeoIndex:GeoIndex+NumMaxNeighbors-1,:) = ...
         maxconnect(sortID,:);
      GeoIndex = GeoIndex + NumMaxNeighbors;
   end
end

% Now add a third column to maxconnect that holds the rank of the
% geodesic in that row.  The rank is from shortest to longest.  So the
% shortest geodesic will get a 1, the next longest a 2, etc.  This will be
% used for the second step function in the GeodesicPlot code.
SortTemp = vertcat(maxconnect{:,1});
SortTemp(:,1) = 1:length(SortTemp);
SortTemp=sortrows(SortTemp,3);
SortTemp(1:2:end-1,2)=1:size(SortTemp,1)/2;
SortTemp(2:2:end,2)=1:size(SortTemp,1)/2;
SortTemp=sortrows(SortTemp,1);
maxconnect(:,3)=num2cell(SortTemp(:,2));

GeoIndex = 1;
for i = 1:length(maxclass)
   if ~isempty(maxclass(i).nbormax)
      NumMaxNeighbors = size(maxclass(i).nbormaxid,1);
      maxclass(i).geodesics(:,3) = maxconnect(GeoIndex:...
                                       GeoIndex+NumMaxNeighbors-1,3);
      GeoIndex = GeoIndex + NumMaxNeighbors;
   end
end

end% maxconnect function

function [L, T] = Plength(P,X)
% find the length of the longest edge in the path return that length as L.
% Also return the lengths of all of the edges sorted return as T.
D = X(P(2:end),:)-X(P(1:end-1),:);
T = sort(sqrt(sum(D.^2,2)),'descend');
L = max(T);
end %Plength function