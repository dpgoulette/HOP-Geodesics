function [GoodMaxGeodesics, maxconnect] = SelectGeodesics(maxconnect,maxclass)
% Fix the comments!
%    Need:
%        1) Pick the longest bar from 0 to n
%        2) Pick the longest from 1 to n (not 0).
%    And make the function flexible to accomodate other selection schemes.
%    Make the output consistent.  This can be done with a 0, 1 or 2 in the
%    fourth column of maxconnect.
%

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

% Get the plot options from the user
fprintf('Which geodesic selection scheme do you want:\n')
while true
   fprintf('    1) Longest bar from 0 to n.\n');
   fprintf('    2) Longest bar from 1 to n.\n')
   select_option = input('Which of the above options do you prefer: ');
   if select_option==1 || select_option==2
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end


GoodMaxGeodesics = zeros(length(maxclass),1);
GeoIndex = 1;

switch select_option % two cases for select_option
   case select_option == 1
      
      for Max = 1:length(maxclass)
         % There are rare cases where a maxclass has no geodesics connected
         % to it.  This only happens on the boundary of the data set.  This
         % next if statment avoids that case.
         if ~isempty(maxclass(Max).nbormaxid)
            NumMaxNeighbors = length(maxclass(Max).nbormaxid);
            % this block selecs the longest bar from 0 to n (so there can
            % be no geodesics selected by the max.
            temp = zeros(NumMaxNeighbors+1,1);
            
            % Note that the for loop starts at 0 for the zero bar
            for Bar = 0:NumMaxNeighbors
               if Bar == 0 %the first bar
                  %bar length is the rank of the current geodesic
                  temp(Bar+1) = maxconnect{GeoIndex,3};
%                   GeoIndex = GeoIndex +1;
               elseif Bar < NumMaxNeighbors %not the first bar AND not the last
                  %bar length is the rank of the next geodesic minus the
                  %rank of the current geodesic.
                  temp(Bar+1) = maxconnect{GeoIndex+1,3} - maxconnect{GeoIndex,3};
                  GeoIndex = GeoIndex + 1;
               else %the last bar.  So Nbor == NumMaxNeighbors +1
                  %bar length is the rank of the longest geodesic minus the
                  %rank of the current geodesic.
                  temp(Bar+1) = size(maxconnect,1)/2 - maxconnect{GeoIndex,3};
                  GeoIndex = GeoIndex + 1;
               end
            end
            [~,slot] = max(temp);
            NumGeosSelected = slot - 1;
            % Currently not using this.  Consider removing.
            GoodMaxGeodesics(Max) = NumGeosSelected;
            
            % Now throw a 1 in the fourth column for each geodesic that was
            % selected in the above process.
            for i = GeoIndex-NumMaxNeighbors:GeoIndex-1
               if i < GeoIndex - NumMaxNeighbors + NumGeosSelected
                  maxconnect{i,4} = 1;
               else
                  maxconnect{i,4} = 0;
               end
            end
         end % if ~empty
      end %loop throught maxclass
      
   otherwise % select_option == 2
      
      for Max = 1:length(maxclass)
         if ~isempty(maxclass(Max).nbormaxid)
            NumMaxNeighbors = length(maxclass(Max).nbormaxid);
            % This block selects the longest of bars 1 through n (not 0).
            % So this max must select at least the shortest geodesic
            % attached to it.
            temp = zeros(NumMaxNeighbors,1);
            for Bar = 1:NumMaxNeighbors
               if Bar < NumMaxNeighbors % not the last bar
                  temp(Bar) = maxconnect{GeoIndex+1, 3} - maxconnect{GeoIndex, 3};
                  GeoIndex = GeoIndex + 1;
               else %the last bar.
                  % So Nbor == NumMaxNeighbors, thus the length of the last bar is
                  % equal to the number of geodesics minus the rank of the current
                  % geodesic.
                  temp(Bar) = size(maxconnect,1)/2 - maxconnect{GeoIndex,3};
                  GeoIndex = GeoIndex + 1;
               end
            end
            [~,NumGeosSelected] = max(temp);
            
            % Not using this currently.  Consider removing.
            GoodMaxGeodesics(Max) = NumGeosSelected;
            
            % Now throw a 1 in the fourth column for each geodesic that was
            % selected in the above process.
            for i = GeoIndex-NumMaxNeighbors:GeoIndex-1
               if i < GeoIndex-NumMaxNeighbors+NumGeosSelected
                  maxconnect{i,4} = 1;
               else
                  maxconnect{i,4} = 0;
               end
            end
         end% if ~empty max neighbors
      end% loop through maxclass
      
end %switch select_option

% Issue!  How do we handle the one sided geodesics?  Currently, in a case
% where max A and B are neighobrs and A selects their geodesic but B does
% not, there will be a 1 in the fourth column for max A and a 0 in the
% fourth column for max B.  If both A and B select the geodesic, then there
% will be a 2 on both ends.

% Now find which geodesics were selected from both ends.  If a geodesic has
% been selected from both ends, store a 2 in maxconnect{i,4} and also do
% the same for the max at the other end of the geodesic.
temp = vertcat(maxconnect{:,1});
temp(:,3) = [];
for i = 1:size(maxconnect,1)
   if (temp(i, 1) < temp(i, 2)) && maxconnect{i,4} == 1
      b = ismember(temp, circshift(temp(i,:), [0, 1]),'rows');
      if maxconnect{b,4} == 1 % then both ends selected this geodesic
         maxconnect{i,4} = 2;
         maxconnect{b,4} = 2;
      end
   end
end

end%function

