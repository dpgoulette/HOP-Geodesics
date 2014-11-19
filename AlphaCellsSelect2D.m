function [GoodEdges,one_cells] = AlphaCellsSelect2D(DT,GoodEdges,VV,VC,GoodIndex)
%

% Already have good delaunay edges (1-cells).  Calculate the value of
% alpha for each good edge

alpha_complex_option = 1;
while alpha_complex_option == 1;
   AllGoodEdges = GoodEdges;
   fprintf('Calculating alpha for the 1-cells.\n\n')
   edge_alpha = EpsilonOneCells2d(DT,AllGoodEdges,VV,VC);
   one_cells = edge_alpha;
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
   
   % User selects which alpha complex selection scheme they want.
   fprintf('Which alpha 1-cell selection scheme do you want?\n')
   while true
      fprintf('   1) Choose a global alpha for all edges (as a percentage).\n')
      fprintf('   2) Standard deviations above the mean edge alpha.\n')
      fprintf('   3) Our parameter free local alpha select.\n')
      alpha_select_option = input('Choose one of the above: ');
      if alpha_select_option == 2 || alpha_select_option == 1 ||...
            alpha_select_option == 3
         break
      else
         fprintf('ERROR! You must enter 1 or 2.\n\n')
      end
   end
   fprintf('\n')
   switch alpha_select_option
      case 1
         % Here we choose a global alpha as a percentage.  So .90 would keep
         % the shortest 90% of the alpha edges.
         while true
            fprintf('\n This selection option allows you to keep the shortest')
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
         edge_alpha = sortrows(edge_alpha,3);
         num_keep = floor(size(edge_alpha,1) * keep_percent);
         if num_keep <= 0
            fprintf('\nAll edges will be removed!!\n')
            AllGoodEdges =[];
         elseif num_keep >= size(edge_alpha,1)
            fprintf('\nNo edges will be deleted. This is the full Delaunay.\n')
         else
            edge_alpha(num_keep + 1:end,:) = [];
            AllGoodEdges = edge_alpha(:,[1,2]);
         end
         
      case 2
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
            AllGoodEdges =[];
         elseif isempty(cutoff_start)
            fprintf('\nNo edges will be deleted. This is equivalent to')
            fprintf(' the full Delaunay.\n')
         else
            edge_alpha(cutoff_start:end,:) = [];
            AllGoodEdges = edge_alpha(:,[1,2]);
         end
      otherwise
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
                  size(AllGoodEdges, 1)] - ...
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
         E = unique(TT,'rows');
         
         AllGoodEdges = E;
         
         %       T(T(:,5)==1,:) = [];
         %       TT = T(:,[1,2]);
         %       TT = sort(TT,2);
         %       E = unique(TT,'rows');
   end
   
   if isempty(AllGoodEdges)
      fprintf('\n\n THERE ARE NO EDGES!\n\n')
      pause(.5)
      fprintf('\n\n THERE ARE NO EDGES!\n\n')
      pause(.5)
      fprintf('\n\n THERE ARE NO EDGES!\n\n')
      pause(2)
      figure('name','The 1-skeleton of the alpha complex you chose.')
      P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
   else 
      % plot the 1 skeleton of the alpha complex so the user can see the result
      figure('name','The 1-skeleton of the alpha complex you chose.')
      P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
      axis equal
      hold on
      axis tight
      X=[DT.X(AllGoodEdges(:,1),1)';DT.X(AllGoodEdges(:,2),1)'];
      Y=[DT.X(AllGoodEdges(:,1),2)';DT.X(AllGoodEdges(:,2),2)'];
      alpha_plot = plot(X,Y,'b');
   end
   if alpha_select_option == 1
      S = sprintf('Percentile of alpha edges kept: %0.1f',keep_percent);
   elseif alpha_select_option == 2
      S = sprintf('Kept all edges below %0.2f std. devs. from the mean.',...
         standard_devs);
   else
      S = sprintf('Result of the current selection scheme.');
   end
   title(S)
   
   fprintf('\nDo you want to redo the alpha complex selection?\n')
   while true
      fprintf('   1) Yes.\n')
      fprintf('   2) No.\n')
      Choice = input('Choose one of the above: ');
      if Choice == 2 || Choice == 1
         break
      else
         fprintf('ERROR! You must enter 1 or 2.\n\n')
      end
   end
   if Choice == 2
      break
   end
end% while alhpa_complex_option
end % main function