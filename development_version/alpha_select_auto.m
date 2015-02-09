function edge_alpha = alpha_select_auto(edge_alpha, DT, GoodIndex)
% This function is experimental and not fully tested yet.  Results
% not guaranteed.

fprintf('\nTEST VERSION PLOTTING!\n')
fprintf('This will cycle through each point 1-by-1 and plot the edges ')
fprintf('in order of alpha from smallest to largest, and it will also ')
fprintf('plot the persistence step functions if you want it to.\n')
fprintf('\n   Do you want to plot the step functions for each point?\n')
while true
   fprintf('      1) Yes\n')
   fprintf('      2) No\n')
   plot_select = input('Choose one of the above: ');
   if plot_select == 2 || plot_select == 1
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
fprintf('\n')

if plot_select == 1
   % Do the parameter free selection scheme.
   AlphaEdgesWithStepFunc(DT,GoodIndex,edge_alpha)
end
% We need each edge expressed both ways.
temp = [fliplr(edge_alpha(:, [1 2])), edge_alpha(:,3)];
simplices_1D = sortrows([edge_alpha; temp],1);

% Now we sort the edges in simplices_1D so that the edges attached
% to each good data point is in order from smallest alpha to largest
% alpha. We store the number of neighbors each point has in
% num_neighbors. NOTE! the ith entry in num_neighbors corresponds to
% the ith entry in GoodIndex, not DT.X.
num_neighbors = NaN(size(DT.X,1),1);
edge_index = 1;
temp = simplices_1D(:,1);
for i = 1:length(GoodIndex)-1
   num_neighbors(GoodIndex(i)) =...
      find(temp(edge_index:end) ~= temp(edge_index), 1)-1;
   [~,sortID] = sort(simplices_1D(edge_index : ...
      edge_index + num_neighbors(GoodIndex(i)) - 1, 3));
   sortID = sortID + edge_index - 1;
   simplices_1D(edge_index :...
      edge_index + num_neighbors(GoodIndex(i)) - 1,:) = ...
      simplices_1D(sortID,:);
   edge_index = edge_index + num_neighbors(GoodIndex(i));
end
% The last entry in GoodIndex has to be outside the loop.
num_neighbors(GoodIndex(end)) = length(temp) + 1 - edge_index;
[~,sortID] = sort(simplices_1D(edge_index :...
   edge_index + num_neighbors(GoodIndex(end)) - 1, 3));
sortID = sortID + edge_index - 1;
simplices_1D(edge_index : end,:) = simplices_1D(sortID,:);
%%% done sorting alpha edges %%%%

% Add a fourth column to simplices_1D to hold the rank of alpha for each
% edge.  The rank is in order from smallest alpha to largest.  Fill
% the fourth column with the rank of that row.
temp = [simplices_1D(:,3), (1:size(simplices_1D, 1))',...
   zeros(size(simplices_1D(:,3)))];
temp = sortrows(temp,1);
temp(1:2:end-1,3) = (1: size(simplices_1D, 1)/2)';
temp(2:2:end,3) = (1: size(simplices_1D, 1)/2)';
temp = sortrows(temp,2);
simplices_1D(:,4) = temp(:,3);
%%% done calculating the ranks and storing them %%%

% Add a fifth column to simplices_1D to hold the selection scheme
% index.
simplices_1D(:,5) = zeros(size(simplices_1D,1),1);


% Create a cell with an entry for each good point and store the rows
% of simplices_1D that correspond to each point
simplex_attachments_1D = cell(size(DT.X(:,1)));
edge_index = 1;
for i = 1:length(GoodIndex)
   simplex_attachments_1D{GoodIndex(i)} = simplices_1D(edge_index : ...
      edge_index + num_neighbors(GoodIndex(i)) - 1, :);
   edge_index = edge_index + num_neighbors(GoodIndex(i));
end

% The selection scheme.  Longest bar from 1 to n.  Store which edges
% were chosen in the 5th column.
for p = 1:length(simplex_attachments_1D)
   % The empty entries are bad points so skip them
   if ~isempty(simplex_attachments_1D{p})
      bar_lengths = [simplex_attachments_1D{p}(2:end,4);...
         size(edge_alpha, 1)] - ...
         simplex_attachments_1D{p}(:,4);
      [~,num_edges_selected] = max(bar_lengths);
      simplex_attachments_1D{p}(1:num_edges_selected,5) = 1;
   end
end

% Now have an option that plots all of the ONCE chosen edges as well
% as the edges that are twice chosen.

for p = 1: length(simplex_attachments_1D)
   for e = 1:size(simplex_attachments_1D{p},1)
      if simplex_attachments_1D{p}(e,1)<simplex_attachments_1D{p}(e,2)...
            && simplex_attachments_1D{p}(e,5) == 1
         nborID = simplex_attachments_1D{p}(e,2);
         nbor_nbors = simplex_attachments_1D{nborID}(:,2);
         r = find(nbor_nbors == simplex_attachments_1D{p}(1,1));
         if simplex_attachments_1D{nborID}(r,5) == 1
            simplex_attachments_1D{p}(e,5) = 2;
            simplex_attachments_1D{nborID}(r,5) = 2;
         end
      elseif simplex_attachments_1D{p}(e,5) == 0
         break
      else
         % Then we have already cataloged this one earlier so skip.
         continue
      end
   end
end

T = vertcat(simplex_attachments_1D{:});
T(T(:,5)==0,:) = [];
TT = T(:,[1,2]);
TT = sort(TT,2);

edge_alpha = unique(TT,'rows');

end % alpha_select_auto