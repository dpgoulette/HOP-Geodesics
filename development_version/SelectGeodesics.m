function maxclass = SelectGeodesics(maxclass,hop)
% Fix the comments!
%    Need:
%        1) Pick the longest bar from 0 to n
%        2) Pick the longest from 1 to n (not 0).
%    And make the function flexible to accomodate other selection schemes.
%    Make the output consistent.  This can be done with a 0, 1 or 2 in the
%    fourth column of maxconnect.


% Issue!  How do we handle the one sided geodesics?  Currently, in a case
% where max A and B are neighobrs and A selects their geodesic but B does
% not, there will be a 1 in the fourth column for max A and a 0 in the
% fourth column for max B.  If both A and B select the geodesic, then there
% will be a 2 on both ends.
%

%  SelectGeodesics - selects which geodesics will be included in the model.
%     The user is prompted to choose a selection scheme they wish to use.
%     The maxclass struct is updated with the selection information.
%
%        input:
%              maxclass - the maxclass struct
%              hop - the hop data struct
%
%        output:
%              maxclass
%
%     Currently this function has three selection scheme options:
%
%        1) The cutoff for geodesic inclusion is based on the longest bar in
%           the step function (the longest bar from 0 to n).
%
%        2) The cutoff for geodesic inclusion is based on the longest bar
%           in the step function NOT including the zero bar. (the longest
%           bar from 1 to n).
%
%        3) Keep all geodesics.
%
%     The function adds entries to the fourth column of
%     maxclass(i).geodesics, where i is an index of a particular max. In
%     particular, let i represent a max and let j represent the j-th
%     geodesic connected to max i, then when the function terminates:
%
%        maxclass(i).geodesics{j, 4} == 0 
%            if this geodesic is not included by i.
%
%        maxclass(i).geodesics{j, 4} == 1 
%            if this geodesic IS included by i but NOT included by the max at
%            the other end of the geodesic.
%
%        maxclass(i).geodesics{j, 4} == 2
%            if this geodesic is included by i AND also by the max at the
%            other end.
%


% Get the geodesic selection scheme from the user.  This will mostly effect
% plotting, but it could be important for data analysis after running the
% model simulation.
fprintf('\n\nWhich geodesic selection scheme do you want:\n')
while true
   fprintf('    1) Longest bar from 0 to n.\n');
   fprintf('    2) Longest bar from 1 to n.\n');
   fprintf('    3) Keep all geodesics.\n');
   select_option = input('Which of the above options do you prefer: ');
   if select_option==1 || select_option==2 || select_option==3
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end

% get the total number of geodesics
temp = vertcat(maxclass.geodesics);
total_num_geodesics = size(temp, 1) / 2;
      
% The geodesic selection functions are separated into dependent functions
% (which are in this file below the main parent function).  Except the last
% case where we include EVERY geodesic.
switch select_option 
   case 1
      maxclass = zero_to_n(maxclass, total_num_geodesics);
   case 2
      maxclass = one_to_n(maxclass, total_num_geodesics);
   
   case 3
      % we are keeping every geodesic. So put a 2 in the fourth column for
      % each geodesic.
      for a = 1:length(maxclass)
         [maxclass(a).geodesics{:,4}] = deal(2);
      end
      
end %switch select_option

% This final code block finds which geodesics were selected from both ends.
% If a geodesic has been selected from both ends, store a 2 in
% maxconnect{i,4} and also do the same for the max at the other end of the
% geodesic.
if select_option == 1 || select_option == 2
   for i = 1:length(maxclass)
      if isempty(maxclass(i).nbormax)
         % Then there are no geodesics connected to this max.  Skip to the
         % next.
         continue
      end
      
      for j = 1:size(maxclass(i).geodesics, 1)
         % Get the current neighbor max index.  We need the index of the
         % point as well as the index of the max.
         
         % Get the point-index of the neighbor max (index into DT.X)
         neighbor_point = maxclass(i).geodesics{j, 2}(2);
         % Get the maxclass index of the neighbor max (index into maxclass)
         curr_n = hop(neighbor_point).maxclassid;
         
         % The next block checks to see if the maxima at both ends of the
         % geodesic included the current geodesic j. Check that the current
         % max has index smaller than neigbor_point (that way we only do
         % the checking from one end of the geodesic, when we get to the
         % max at the other end we will know that we already did it).  Also
         % check that the 4th column has a 1, i.e. the geodesic was
         % included by the i-th max.
         if maxclass(i).geodesics{j,2}(1) < neighbor_point &&...
               maxclass(i).geodesics{j,4} == 1
            % So we haven't done this geodesic yet and this geodesic was
            % selected.
            
            % Find the row, r, that curr_n holds this same geodesic in.
            temp = vertcat( maxclass(curr_n).geodesics{:, 2} );
            r = temp(:,2) == maxclass(i).max; % logical index
            
            if maxclass(curr_n).geodesics{r, 4} == 1
               % Thus, both ends of the geodesic selected this geodesic.
               maxclass(i).geodesics{j, 4} = 2;
               maxclass(curr_n).geodesics{r, 4} = 2;
            end
         end
      end
   end
end

end % main function: SelectGeodescis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% DEPENDENT FUNCTIONS %%%%%%%%%%%%%%
%               | |       | |
%               | |       | |
%              \   /     \   /
%               \ /       \ /

function maxclass = zero_to_n(maxclass, total_num_geodesics)
%  zero_to_n -- loops through each max in maxclass and determines which
%        geodesics to keep and which to throw away.  We rank all of the
%        geodesics in order of their length (based on our metric).  The
%        rank of the shortest geodesic in ALL of maxclass has rank 1 and
%        the longest geodesic has rank r = "the total number of geodesics."
%        The selection scheme for the i-th max is based on the rank order
%        of the lengths of the geodesics connected to max i.  We use a
%        persistence-like scheme to find out the cut-off point for which
%        geodesics are connected to the i-th max and which are removed. For
%        example, if all of the geodesics connected to i are very short (in
%        comparison to all of the geodesics in maxclass), then this max is
%        close to it's neighbor max. So it is in a denser region, and hence
%        we keep all of it's geodesics.  But suppose, for a different
%        example, max i has 6 geodesics connected to it.  Further suppose
%        that the two shortest geodesics are very short, but the longest 4
%        are very long, then there will be a large gap in the rank between
%        the 2nd geodesic and the 3rd (when taken in order of length).  Our
%        algorithm finds the largest gap in the ranks and keeps all of the
%        geodesics before that gap.  So in this example we would keep the 2
%        shortest and not include the 4 longest.  NOTE! Since there are n
%        geodesics connected to max i, there are n+1 gaps to consider.
%        This is because we compare the rank of the shortest geodesic to 1
%        and we compare the rank of the longest geodesic to r = "the total
%        number of geodesics.
%
%        zero_to_n updates the fourth column in maxclass(i).geodescis for
%        each geodesic connected to i.

for Max = 1:length(maxclass)
   % There are cases where a maxclass has no geodesics connected to
   % it.  This only happens on the boundary of the data set or if you
   % hop on an alpha complex. Skip to the next max if this is the
   % case.
   if isempty(maxclass(Max).nbormaxid)
      continue
   end
   NumMaxNeighbors = length(maxclass(Max).nbormaxid);
   
   % The next section selects the longest bar from 0 to n (so there might
   % be no geodesics selected by the max).
   temp = zeros(NumMaxNeighbors + 1, 1);
   % Note that the for loop starts at 0 for the zero bar
   for Bar = 0 : NumMaxNeighbors
      if Bar == 0 %the first bar
         %bar length is the rank of the current geodesic
         temp(Bar + 1) = maxclass(Max).geodesics{Bar + 1, 3};
         %                   GeoIndex = GeoIndex +1;
      elseif Bar < NumMaxNeighbors %not the first bar AND not the last
         %bar length is the rank of the next geodesic minus the
         %rank of the current geodesic.
         temp(Bar+1) = maxclass(Max).geodesics{Bar + 1, 3} - ...
            maxclass(Max).geodesics{Bar, 3};
      else %the last bar.  So Nbor == NumMaxNeighbors +1
         %bar length is the rank of the longest geodesic minus the
         %rank of the current geodesic.
         temp(Bar+1) = total_num_geodesics - ...
            maxclass(Max).geodesics{Bar, 3};
      end
   end
   [~,slot] = max(temp);
   NumGeosSelected = slot - 1;
   
   % Now put a 1 in the fourth column for each geodesic that was
   % selected in the above process.  Put 0 for the rest that weren't
   % selected.
   iter = 1;
   while iter <= NumGeosSelected
      maxclass(Max).geodesics{iter, 4} = 1;
      iter = iter + 1;
   end
   while iter <= NumMaxNeighbors
      maxclass(Max).geodesics{iter, 4} = 0;
      iter = iter + 1;
   end
   
end %loop through maxclass

end % function: zero_to_n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maxclass = one_to_n(maxclass, total_num_geodesics)
%  one_to_n -- loops through each max in maxclass and determines which
%        geodesics to keep and which to throw away.  We rank all of the
%        geodesics in order of their length (based on our metric).  The
%        rank of the shortest geodesic in ALL of maxclass has rank 1 and
%        the longest geodesic has rank r = "the total number of geodesics."
%        The selection scheme for the i-th max is based on the rank order
%        of the lengths of the geodesics connected to max i.  We use a
%        persistence-like scheme to find out the cut-off point for which
%        geodesics are connected to the i-th max and which are removed. For
%        example, if all of the geodesics connected to i are very short (in
%        comparison to all of the geodesics in maxclass), then this max is
%        close to it's neighbor max. So it is in a denser region, and hence
%        we keep all of it's geodesics.  But suppose, for a different
%        example, max i has 6 geodesics connected to it.  Further suppose
%        that the two shortest geodesics are very short, but the longest 4
%        are very long, then there will be a large gap in the rank between
%        the 2nd geodesic and the 3rd (when taken in order of length).  Our
%        algorithm finds the largest gap in the ranks and keeps all of the
%        geodesics before that gap.  So in this example we would keep the 2
%        shortest and not include the 4 longest.  NOTE! Since there are n
%        geodesics connected to max i, there are n+1 gaps to consider. This
%        is because we compare the rank of the shortest geodesic to 1 and
%        we compare the rank of the longest geodesic to r = "the total
%        number of geodesics.  HOWEVER, in this function we do not consider
%        gap zero (the difference in rank between the shortest geodesic
%        connected to i and the shortest geodesic in all of maxclass).
%        Thus, this function will include AT LEAST ONE GEODESIC FOR EVERY
%        MAX IN MAXCLASS.  After this function terminates, the parent
%        function, SelectGeodesics, will only include geodesics that were
%        selected from BOTH ends (meaning the maxima at both ends of the
%        geodesic selected the geodesic connecting them).
%
%        one_to_n updates the fourth column in maxclass(i).geodescis for
%        each geodesic connected to i.


for Max = 1:length(maxclass)
   % There are cases where a maxclass has no geodesics connected to
   % it.  This only happens on the boundary of the data set or if you
   % hop on an alpha complex. Skip to the next max if this is the
   % case.
   if isempty(maxclass(Max).nbormaxid)
      continue
   end
   NumMaxNeighbors = length(maxclass(Max).nbormaxid);
   % This block selects the longest of bars 1 through n (not 0).
   % So this max must select at least the shortest geodesic
   % attached to it.
   temp = zeros(NumMaxNeighbors, 1);
   % Note that the for loop starts at 1 since we are not considering
   % the zero bar. (Note: this means that every max will select at
   % least one geodesic.  This is o.k. because, in the end, we will
   % only include geodesics that are selected from BOTH ends.)
   for Bar = 1 : NumMaxNeighbors
      if Bar < NumMaxNeighbors %not the first bar AND not the last
         %bar length is the rank of the next geodesic minus the
         %rank of the current geodesic.
         temp(Bar) = maxclass(Max).geodesics{Bar + 1, 3} - ...
            maxclass(Max).geodesics{Bar, 3};
      else %the last bar.  So Nbor == NumMaxNeighbors +1
         %bar length is the rank of the longest geodesic minus the
         %rank of the current geodesic.
         temp(Bar) = total_num_geodesics - ...
            maxclass(Max).geodesics{Bar, 3};
      end
   end
   [~,slot] = max(temp);
   NumGeosSelected = slot;
   
   % Now put a 1 in the fourth column for each geodesic that was
   % selected in the above process.  Put 0 for the rest that weren't
   % selected.
   iter = 1;
   while iter <= NumGeosSelected
      maxclass(Max).geodesics{iter, 4} = 1;
      iter = iter + 1;
   end
   while iter <= NumMaxNeighbors
      maxclass(Max).geodesics{iter, 4} = 0;
      iter = iter + 1;
   end
   
end %loop through maxclass

end % function: one_to_n

