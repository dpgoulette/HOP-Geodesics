function [alpha_edges, standard_devs] = alpha_select_stdev(edge_alpha)
% alpha_select_stdev -- Removes the edges that have an alpha value larger
%     than the user defined threshold.  This function prompts the user to
%     enter a threshold, p, as a decimal percentage (number in
%     interval[0,1]).  Then the below that percentile are kept.  Those
%     above are deleted.  If the user inputs .8, then the edges with the
%     smallest 80% of alpha values are kept and the largest 20% are thrown
%     away.  
%
%     input:
%        edge_alhpa - k by 3 array. First two columns hold the enpoints of
%                     each edge. The third column holds the alpha value for
%                     the edge in that row.
%
%        alpha_edges - 

fprintf('\nHow many standard deviations above the mean do you want')
fprintf(' to keep? Entering a negative number will be below the mean.')
fprintf(' Every edge with alpha less than or equal to the number of')
fprintf(' standard deviations will be kept.\n')
standard_devs = input('   Enter any real number: ');

edge_alpha = sortrows(edge_alpha,3);
SD = std(edge_alpha(:,3));
cutoff = mean(edge_alpha(:,3)) + standard_devs * SD;
cutoff_start = find(edge_alpha(:,3) > cutoff, 1);

if cutoff_start == 1
   fprintf('\nAll edges will be removed!! So evry point is isolated.\n')
   alpha_edges =[];
   pause(2)
   
elseif isempty(cutoff_start)
   fprintf('\nNOTE!!!\nNo edges will be deleted. This is equivalent to')
   fprintf(' the full Delaunay.\n')
   alpha_edges = edge_alpha(:, [1,2]);
   pause(2)
   
else
   % Delete the rows that don't pass the threshold cutoff.
   edge_alpha(cutoff_start:end,:) = [];
   
   % Return the edges that passed the threshold.
   alpha_edges = edge_alpha(:, [1,2]);
   
end

end % std_dev_threshold