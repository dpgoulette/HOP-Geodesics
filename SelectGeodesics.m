function GoodMaxGeodesics = SelectGeodesics(maxconnect,maxclass)
%  SelectGeodesics - selects which geodesics will be included.  The
%     function returns a vector with dimension: length(maxclass) x 1.  The
%     ith entry in this vector is the number of geodesics connected to the
%     ith max in maxclass that we will keep (the rest we will throw away).
%     So if the entry is 3, then we will keep the three shortest geodesics
%     incident to the ith max. The function analyzes the differences
%     between the ranks of the geodesics (from shortest to longest) that
%     are connected to each max.  To see a plot of the rank order of each
%     geodesic, run GeodesicPlot.m after running this function and see the
%     bottom right plot.  (Note that the geodesics in max class are in
%     blocks corresponding to each max.  And within each block they are
%     sorted in order of length from shortest to longest. Note that the
%     third column of maxclass holds the rank order of the geodesics; they
%     are ranked by length. Remember that the length of a geodesic is equal
%     to the length of the longest edge in the geodesic path.)  The rank of
%     the shortest geodesic is 1, the next shortest one is 2, etc.
%
%     This function searches for the longest bar among bars 1 through n
%     (where n is the number of geodesics connected to the current max).
%     This function does NOT consider the 0 bar.  So at least one geodesic
%     will be preliminarily selected for each max.  But we will then check
%     to see if the geodesic was preliminarily selected by both maxima at
%     its endpoints.  If so, then we will include the geodesic.  If only
%     one end selects the geodesic then we will not.
%

GoodMaxGeodesics = zeros(length(maxclass),1);
GeoIndex = 1;
for Max = 1:length(maxclass)
   if ~isempty(maxclass(Max).nbormaxid)
      NumMaxNeighbors = length(maxclass(Max).nbormaxid);
      %%%%  I'm saving the following code for now because I revert to this.
      %
      %    temp = zeros(NumMaxNeighbors+1,1);
      %    for Nbor = 1:NumMaxNeighbors+1
      %       if Nbor == 1 %the first bar
      %          temp(Nbor) = maxconnect{GeoIndex,3};
      %          GeoIndex = GeoIndex +1;
      %       elseif Nbor < NumMaxNeighbors + 1 %not the first bar AND not the last
      %          temp(Nbor) = maxconnect{GeoIndex,3} - maxconnect{GeoIndex-1,3};
      %          GeoIndex = GeoIndex +1;
      %       else %the last bar.  So Nbor == NumMaxNeighbors +1
      %          temp(Nbor) = size(maxconnect,1)/2 - maxconnect{GeoIndex-1,3};
      %       end
      %    end
      %    [~,slot] = max(temp);
      %    GoodMaxGeodesics(Max) = slot - 1;
      %%%%%%%%%%%%%%%%%%
      
      % This block selects the longest of bars 1 through n (not 0).
      temp = zeros(NumMaxNeighbors,1);
      for Nbor = 1:NumMaxNeighbors
         if Nbor < NumMaxNeighbors % not the last bar
            temp(Nbor) = maxconnect{GeoIndex+1, 3} - maxconnect{GeoIndex, 3};
            GeoIndex = GeoIndex + 1;
         else %the last bar.
            % So Nbor == NumMaxNeighbors, thus the length of the last bar is
            % equal to the number of geodesics minus the rank of the current
            % geodesic.
            temp(Nbor) = size(maxconnect,1)/2 - maxconnect{GeoIndex,3};
            GeoIndex = GeoIndex + 1;
         end
      end
      [~,slot] = max(temp);
      GoodMaxGeodesics(Max) = slot;
   end
end

end%function
