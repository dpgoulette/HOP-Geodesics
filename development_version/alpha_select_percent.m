function [alpha_edges, keep_percent] = alpha_select_percent(edge_alpha)

% In this selection scheme we choose a global alpha threshold as a
% percentage.  So, for example, .90 would keep the shortest 90% of
% the alpha edges.

% Get the percentage from the user.
while true
   fprintf('\nThis selection option allows you to keep the shortest')
   fprintf(' alpha values up to a certain percentile.')
   fprintf('\nWhat percent of the smallest edge alpha values do you want ')
   fprintf('to keep? \n')
   keep_percent = input('Enter a decimal from 0 to 1: ');
   if keep_percent < 0 || keep_percent > 1
      fprintf('\nERROR: The number must be in the interval [0,1].\n')
   else
      fprintf('\n\nSo the smallest %f percent',keep_percent * 100)
      fprintf(' of alpha values will be chosen.\n\n')
      pause(1)
      break
   end
end

% Find how many of the edges pass the threshold.
edge_alpha = sortrows(edge_alpha,3);
num_keep = floor(size(edge_alpha,1) * keep_percent);

% Now delete the edges that do not pass the threshold. (Note: it
% is fine that edge_alpha was row-sorted on the alpha values in
% the 3rd column. Later functions that depend on these edges will
% resort this matrix anyway.)
if num_keep <= 0
   fprintf('\nAll edges will be removed!!\n')
   alpha_edges =[];
elseif num_keep >= size(edge_alpha,1)
   fprintf('\nNOTE!\nNo edges will be deleted. This is the full Delaunay. ')
   fprintf('triangulation.\n')
   pause(1)
   alpha_edges = edge_alpha(:, [1,2]);
else
   % Delete the rows that don't pass the threshold cutoff.
   edge_alpha(num_keep + 1:end, :) = [];
   
   % Return the edges that passed the threshold.
   edge_alpha(:,3) = [];
   alpha_edges = edge_alpha(:, [1,2]);
end

end % alpha_select_percent