function maxconnect = hopmaxconnect(points,maxclass,hop)
%  FIX - if two neighbor max paths have the same max-length edge, then we
%     need to check the second longest, and then the third etc, until the tie
%     is broken.  Currently the code only checks the longest edge.
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
      B=maxclass(a).bonds;%B holds the bonds of maxclass a
      for b=1:size(maxclass(a).nbormaxid,1)
         if maxclass(a).nbormaxid(b) < a
            %thus we have already cataloged this geodesic in a previous
            %iteration.  So we simply need to flip it around and save it in
            %the current max class iteration.
            r = ismember(tempedges,[maxclass(maxclass(a).nbormaxid(b)).max,...
               maxclass(a).max],'rows');
            maxconnect{k,1}(3)=maxconnect{r,1}(3);
            maxconnect{k,2}=fliplr(maxconnect{r,2});
            maxconnect{k,1}([1 2])= fliplr(maxconnect{r,1}([1 2]));
            tempedges(k,:) = maxconnect{k,1}([1 2]);
            k=k+1;
         else
            %Bn will hold the bonds of maxclass b, where maxclass b is
            %a neighbor of maxclass a.
            Bn=maxclass(maxclass(a).nbormaxid(b)).bonds;
            Bn=fliplr(Bn);
            r=ismember(B,Bn,'rows');%find a and b's common bond edges
            temp=B(r,:);%temp holds the common bonds.
            Paths = cell(size(temp,1),1);
            Plengths=zeros(size(temp,1),1);
            
            for c=1:size(temp,1)
               
               %%%%%%%%%%%%%%%%%%
               % This next block was removed along with the "Maxpath"
               % function below.  Not needed after I added the
               % 'hopmaxpath' field in the hop data structure
               %%%%%%%%%%%%%%%%%%
               %                     P1a = Maxpath(temp(c,1),maxpointer,maxclass(a).max);
               %                     P2a = Maxpath(temp(c,2),maxpointer,...
               %                         maxclass(maxclass(a).nbormaxid(b)).max);
               
               P1 = hop(temp(c,1)).hopmaxpath;
               P2 = hop(temp(c,2)).hopmaxpath;
               if isempty(P1)
                  %then temp(c,1) is a boundary point that is also a
                  %max! So it's hop path has no edges.  Thus we must
                  %tack it on to the geodesic by making P1 contain
                  %this point only.  i.e. P1 is a single vertex "path."
                  P1 = temp(c,1);
               end
               if isempty(P2)
                  %same as for P1 above.
                  P2 = temp(c,2);
               end
               
               Paths{c}=[fliplr(P1), P2];
               %now find the euclidean length of the longest edge in
               %the path from max to max.
               Plengths(c)=Plength(Paths{c},points);
            end
            %Now remove the bond paths we just calculated from the B
            %matrix for the next iteration of this for loop.
            B(r,:)=[];
            [maxconnect{k,1}(3), id]=min(Plengths);%minimum max edge length path.
            
            maxconnect{k,2}=Paths{id};
            maxconnect{k,1}([1 2])= [maxconnect{k,2}(1),maxconnect{k,2}(end)];
            
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
% shortest to longest.
GeoIndex = 1;
for i = 1:length(maxclass)
   if ~isempty(maxclass(i).nbormax)
      NumMaxNeighbors = length(maxclass(i).nbormaxid);
      temp = vertcat(maxconnect{GeoIndex:GeoIndex+NumMaxNeighbors-1,1});
      [~,sortID] = sort(temp(:,3));
      sortID = sortID+(GeoIndex-1);
      maxconnect(GeoIndex:GeoIndex+NumMaxNeighbors-1,:) = ...
         maxconnect(sortID,:);
      GeoIndex = GeoIndex + NumMaxNeighbors;
   end
end
% Now add a third column to maxconnectSort that holds the rank of the
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

end% maxconnect function

function L = Plength(P,X)
% find the length of the longest edge in the path
D = X(P(2:end),:)-X(P(1:end-1),:);
T = sqrt(sum(D.^2,2));
L = max(T);
end %Plength function

%%%% Old section.  Not needed any more.  We save the hop paths when doing
%%%% hop so we can just look this up in the data base.  No need to
%%%% calculate it here.
% function Path = Maxpath(p,maxpointer,Max)
% % find the maxpath for the point p.
%
% if p==Max
%     %point k is a max
%     Path = p;
% else
%     Path=zeros(1,100);%way too long.  That is fine.
%     Path(1)=p;
%     for a=2:length(maxpointer)
%         if maxpointer(Path(a-1))==Max %done
%             Path(a)=maxpointer(Path(a-1));
%             Path(a+1:end)=[];
%             break
%         else %the path continues;
%             Path(a)=maxpointer(Path(a-1));
%         end
%     end
% end
% end
