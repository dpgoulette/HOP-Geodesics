function GoodMaxGeodesics = SelectGeodesics(maxconnect,maxclass)
%  LongestBarIndex - COMMENT THIS!!

GoodMaxGeodesics = zeros(length(maxclass),1);
GeoIndex = 1;
for Max = 1:length(maxclass)
   NumMaxNeighbors = length(maxclass(Max).nbormaxid);
   temp = zeros(NumMaxNeighbors+1,1);
   for Nbor = 1:NumMaxNeighbors+1
      if Nbor == 1 %the first bar
         temp(Nbor) = maxconnect{GeoIndex,3};
         GeoIndex = GeoIndex +1;
      elseif Nbor < NumMaxNeighbors + 1 %not the first bar AND not the last
         temp(Nbor) = maxconnect{GeoIndex,3} - maxconnect{GeoIndex-1,3};
         GeoIndex = GeoIndex +1;
      else %the last bar.  So Nbor == NumMaxNeighbors +1
         temp(Nbor) = size(maxconnect,1)/2 - maxconnect{GeoIndex-1,3};
%          GeoIndex = GeoIndex +1;
      end
   end
   [~,slot] = max(temp);
   GoodMaxGeodesics(Max) = slot - 1;
end

end%function
