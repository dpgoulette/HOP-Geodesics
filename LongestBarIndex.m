function PersistMaxNbors = LongestBarIndex(maxconnectSort,maxclass)
%  LongestBarIndex - COMMENT THIS!!

PersistMaxNbors = zeros(length(maxclass),1);
GeoIndex = 1;
for Max = 1:length(maxclass)
   NumMaxNeighbors = length(maxclass(Max).nbormaxid);
   temp = zeros(NumMaxNeighbors+1,1);
   for Nbor = 1:NumMaxNeighbors+1
      if Nbor == 1 %the first bar
         temp(Nbor) = maxconnectSort{GeoIndex,3};
         GeoIndex = GeoIndex +1;
      elseif Nbor < NumMaxNeighbors + 1 %not the first bar AND not the last
         temp(Nbor) = maxconnectSort{GeoIndex,3} - maxconnectSort{GeoIndex-1,3};
         GeoIndex = GeoIndex +1;
      else %the last bar.  So Nbor == NumMaxNeighbors +1
         temp(Nbor) = size(maxconnectSort,1)/2 - maxconnectSort{GeoIndex-1,3};
%          GeoIndex = GeoIndex +1;
      end
   end
   [~,slot] = max(temp);
   PersistMaxNbors(Max) = slot - 1;
end

end%function
