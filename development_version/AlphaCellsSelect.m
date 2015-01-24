function alpha_complex = AlphaCellsSelect(DT,GoodEdges,VV,VC,GoodIndex)
% AlphaCellsSelect -- Selects a subset of the Delaunay 1-skeleton (i.e. the
%     Delaunay graph) based on user input.  The resulting edges in the
%     alpha complex will be returned.  The user can select one of the
%     selection schemes (described below).  The resulting complex is then
%     plotted before the function terminates so the user can decide whether
%     to redo the complex with different options.
%
%     inputs:
%              DT - The delaunay triangulation object
%              GoodEdges - the "good edges" in the delaunay triangulation.
%                          These are edges that do not have data points on
%                          the boundary of the data space. (see comments in
%                          HOPDataPrepare.m for more)
%              VV - Voronoi diagram vertices (the coordinate points)
%              VC - Voronoi cell vertices (indices into VV for each cell)
%              GoodIndex - The "good points" in DT.X (points that arent on
%                          the boundary of the data space.
%     output:
%              alpha_complex - k x 2 matrix with the edges that were
%                          selected.  These edges have an alpha value that
%                          passed the threshold.
%     
%     Note! If you would like the function to return ALL of the edges along
%     with the associated value of alpha for each edge, then use the
%     AlphaOneCells2d or AlphaOneCells3d functions directly.
%
%     These are the selection schemes available in this function (the code
%     for these selection schemes can be found below the main function in
%     this file):
%           1) percentage_threshold
%           2) st_dev_threshold
%           3) step_function_selection
%
%     Details: 
%     Scheme 1 and 2 are the more standard approach of using a
%     global alpha value and including the 1-cells (edges) that pass that
%     threshold.  The third is an attempt at a parameter free, dynamic,
%     local alpha selection process to determine which edges are included.
%
%     1) percentage_threshold
%           Prompts the user for a percentage (decimal in [0,1]).  The
%           edges with alpha below that percentile are kept. Those above
%           are thrown away.  Example: if the user inputs 0.8, then the
%           smallest 80% of alpha values are kept. 
%     2) st_dev_threshold
%           Prompts the user for a percentage (decimal in [0,1]).  The
%           edges with alpha below that percentile are kept. Those above
%           are thrown away.  Example: if the user inputs 0.8, then the
%           smallest 80% of alpha values are kept. 
%     3) step_function_selection
%           This is an experimental parameter free selection scheme that is
%           under development.

dimension = size(DT.X,2);

% Calculate the alhpa value for each 1-cell in GoodEdges.  edge_alhpa adds
% a third column to GoodEdges.  The third column holds the alpha value for
% the edge in that row (the edge endpoints is is column 1 and 2).
fprintf('Calculating alpha for the 1-cells.\n\n')
if dimension == 2
   edge_alpha = AlphaOneCells2d(DT,GoodEdges,VV,VC);
else
   edge_alpha = AlphaOneCells3d(DT,GoodEdges,VV,VC);
end

alpha_complex_option = 1;
while alpha_complex_option == 1;
   
   %       %%%%  THIS BLOCK IS DISABLED CURRENTLY %%%%%%
   %       % If we want to use the 2-cells at some time enable this block and
   %       % make sure that the function returns the cells2 matrix
   %
   %       % We already have the 1-cells (GoodEdges).  We need the 2-cells
   %       DTtris=DT.Triangulation; % 2-cells
   %       % Now remove bad triangles.  They are "bad" if they contain
   %
   %       [R,~] = find(ismember(DTtris, BadDataID));
   %       R = unique(R);
   %       cells2 = DTtris;
   %       cells2(R,:) = [];
   %
   %       % Find the circumcenter radius for the triangles.  This is the value of
   %       % alpha for all of teh 2 cells.
   %       [~,RCC] = circumcenters(DT);
   %       % Now delete the radii from bad triangles
   %       RCC(R)=[];
   %       % Save the 2 cells along with their alpha value in the third column.
   %       cells2=[cells2, RCC];
   %
   %%%%% WE WOULD ALSO NEED SOME SORT OF A SELECTION SCHEME!!! RIGHT??!!
   
   
   % User selects which alpha complex selection scheme they want to use.
   fprintf('Which alpha 1-cell selection scheme do you want? (NOTE: ')
   fprintf('option 3 is an experimental scheme under development.)\n\n')
   while true
      fprintf('   1) Choose a global alpha threshold for all edges')
      fprintf(' (as a percentage).\n')
      fprintf('   2) Standard deviations above the mean edge alpha.\n')
      fprintf('   3) Our parameter free local alpha select.\n\n')
      alpha_select_option = input('Choose one of the above: ');
      if alpha_select_option == 2 || alpha_select_option == 1 ||...
            alpha_select_option == 3
         break
      else
         fprintf('ERROR! You must enter 1 or 2.\n\n')
      end
   end
   fprintf('\n')
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%
   
   switch alpha_select_option
      case 1
         [alpha_complex, keep_percent] = percentage_threshold(edge_alpha);
         
      case 2
         [alpha_complex, standard_devs] = st_dev_threshold(edge_alpha);
         
      otherwise % 3
         alpha_complex = step_function_selection(edge_alpha, DT, GoodIndex);
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN %%%%%%%%%%%%%%%%%%%%%%%%%
   
   % The rest of this function just plots the result of the above alpha complex
   % selection scheme (based on the user input).  The user can then decide
   % whether to keep the chosen complex or redo it with different options.
   
   % plot the results for the user.
   alhpa_complex_plot(DT, GoodIndex, alpha_complex)
   
   % Create and print the title of the alpha complex plot.  This will
   % remind the user what their selection choice was.  Title depends on
   % choice of selection scheme.
   if alpha_select_option == 1
      S = sprintf('Percentile of alpha edges kept: %0.1f',keep_percent);
   elseif alpha_select_option == 2
      S = sprintf('Kept all edges below %0.2f std. devs. from the mean.',...
         standard_devs);
   else
      S = sprintf('Result of the current parameter free selection scheme.');
   end
   title(S)
   
   % This section allows the user to decide if they want to redo the
   % alhpa complex selection or not.
   fprintf('\n\n\n\nBefore we proceed with applying HOP to the current alpha ')
   fprintf('complex, you now have the option to redo the ')
   fprintf('alpha complex selection if you wish.\n\n')
   fprintf('Do you want to keep the current alpha complex?\n')
   while true
      fprintf('   1) Yes. Apply HOP to this alpha complex.\n')
      fprintf('   2) No. I want to redo the alpha complex with different ')
      fprintf('settings.\n\n')
      Choice = input('Choose one of the above: ');
      if Choice == 2 || Choice == 1
         break
      else
         fprintf('ERROR! You must enter 1 or 2.\n\n')
      end
   end
   if Choice == 1
      alpha_complex_option = 0;
   end
end% while alhpa_complex_option == 1

end % main function -- AlphaCellsSelect



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      DEPENDENT FUNCTIONS     %%%%%%%%%%%%%%%


function [edge_alpha_trimmed, keep_percent] = percentage_threshold(edge_alpha)

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
      fprintf('\n\nSo the smallest %f percent of alpha values will be chosen.\n\n',...
         keep_percent * 100)
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
   edge_alpha_trimmed =[];
elseif num_keep >= size(edge_alpha,1)
   fprintf('\nNOTE!\nNo edges will be deleted. This is the full Delaunay. ')
   fprintf('triangulation.\n')
   pause(1)
   edge_alpha_trimmed = edge_alpha;
else
   % Delete the rows that don't pass the threshold cutoff.
   edge_alpha(num_keep + 1:end, :) = [];
   
   % Return the edges that passed the threshold.
   edge_alpha(:,3) = [];
   edge_alpha_trimmed = edge_alpha(:, [1,2]);
end

end % percentage_threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edge_alpha_trimmed, standard_devs] = st_dev_threshold(edge_alpha)
% st_dev_threshold -- Removes the edges that have an alpha value larger
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
%        edge_alpha_trimmed - 

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
   edge_alpha_trimmed =[];
   pause(2)
   
elseif isempty(cutoff_start)
   fprintf('\nNOTE!!!\nNo edges will be deleted. This is equivalent to')
   fprintf(' the full Delaunay.\n')
   edge_alpha_trimmed = edge_alpha;
   pause(2)
   
else
   % Delete the rows that don't pass the threshold cutoff.
   edge_alpha(cutoff_start:end,:) = [];
   
   % Return the edges that passed the threshold.
   edge_alpha_trimmed = edge_alpha(:, [1,2]);
   
end

end % std_dev_threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edge_alpha = step_function_selection(edge_alpha, DT, GoodIndex)
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

end % step_function_selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alhpa_complex_plot(DT, GoodIndex, alpha_complex)
%

dimension = size(DT.X,2);

if isempty(alpha_complex)
   pause(.5)
   fprintf('\n\n THERE ARE NO EDGES!')
   pause(.5)
   fprintf('   THERE ARE NO EDGES!')
   pause(.5)
   fprintf('   THERE ARE NO EDGES!\n\n')
   pause(.5)
   fprintf('Your complex has no edges (1-cells).\n\n')
   pause(2)
   figure('name','The 1-skeleton of the alpha complex you chose.')
   if dimension == 2
      P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
   else
      P1 = plot3(DT.X(GoodIndex,1),DT.X(GoodIndex,2),...
         DT.X(GoodIndex,3),'k.'); % raw data
   end
   axis equal
   axis tight
else
   % plot the 1 skeleton of the alpha complex so the user can see the result
   figure('name','The 1-skeleton of the alpha complex you chose.')
   if dimension == 2
      P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
   else
      P1 = plot3(DT.X(GoodIndex,1),DT.X(GoodIndex,2),...
         DT.X(GoodIndex,3),'k.'); % raw data
   end
   axis equal
   hold on
   axis tight
   if dimension == 2
      X=[DT.X(alpha_complex(:,1),1)';DT.X(alpha_complex(:,2),1)'];
      Y=[DT.X(alpha_complex(:,1),2)';DT.X(alpha_complex(:,2),2)'];
      alpha_plot = plot(X,Y,'b');
   else
      X=[DT.X(alpha_complex(:,1),1)';DT.X(alpha_complex(:,2),1)'];
      Y=[DT.X(alpha_complex(:,1),2)';DT.X(alpha_complex(:,2),2)'];
      Z=[DT.X(alpha_complex(:,1),3)';DT.X(alpha_complex(:,2),3)'];
      alpha_plot = plot3(X,Y,Z,'b');
   end
end

end % alhpa_complex_plot