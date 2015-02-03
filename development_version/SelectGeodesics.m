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
% (except the last case where we include EVERY geodesic).
switch select_option 
   case 1
      maxclass = geo_select_0_to_n(maxclass, total_num_geodesics);
   case 2
      maxclass = geo_select_1_to_n(maxclass, total_num_geodesics);
   
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
