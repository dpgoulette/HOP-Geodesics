function [alpha_complex, edges_with_alpha_all] =...
   AlphaCellsSelect(DT,GoodEdges, VV,VC,GoodIndex)
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
%              alpha_complex - k by 2 matrix with the edges that were
%                          selected.  These edges have an alpha value that
%                          passed the threshold.
%     
%              edges_with_alpha_all - k by 3 matrix.  The first two columns
%                          hold all of the edges from GoodEdges.  The third
%                          column holds the alpha value for each edge.
%
%
%     These are the selection schemes currently available in this function: 
%           1) alpha_select_percent
%           2) alpha_select_stdev
%           3) alpha_select_auto
%
%     Details: 
%     Scheme 1 and 2 are the more standard approach of using a
%     global alpha value and including the 1-cells (edges) that pass that
%     threshold.  The third is an attempt at a parameter free, dynamic,
%     local alpha selection process to determine which edges are included.
%
%     1) alpha_select_percent
%           Prompts the user for a percentage (decimal in [0,1]).  The
%           edges with alpha below that percentile are kept. Those above
%           are thrown away.  Example: if the user inputs 0.8, then the
%           smallest 80% of alpha values are kept. 
%     2) alpha_select_stdev
%           Prompts the user for a percentage (decimal in [0,1]).  The
%           edges with alpha below that percentile are kept. Those above
%           are thrown away.  Example: if the user inputs 0.8, then the
%           smallest 80% of alpha values are kept. 
%     3) alpha_select_auto
%           This is an experimental parameter free selection scheme that is
%           under development.

dimension = size(DT.X,2);

% Calculate the alhpa value for each 1-cell in GoodEdges.  edge_alhpa adds
% a third column to GoodEdges.  The third column holds the alpha value for
% the edge in that row (the edge endpoints is is column 1 and 2).
fprintf('Calculating alpha for the 1-cells.\n\n')
if dimension == 2
   edges_with_alpha_all = AlphaOneCells2d(DT,GoodEdges,VV,VC);
else
   edges_with_alpha_all = AlphaOneCells3d(DT,GoodEdges,VV,VC);
end

done = false;
while ~done
   
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
         [alpha_complex, keep_percent] =...
            alpha_select_percent(edges_with_alpha_all);
         
      case 2
         [alpha_complex, standard_devs] =...
            alpha_select_stdev(edges_with_alpha_all);
         
      otherwise % 3
         alpha_complex =...
            alpha_select_auto(edges_with_alpha_all,DT,GoodIndex);
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN %%%%%%%%%%%%%%%%%%%%%%%%%
   
   % The rest of this function just plots the result of the above alpha complex
   % selection scheme (based on the user input).  The user can then decide
   % whether to keep the chosen complex or redo it with different options.
   
   % plot the results for the user. (dependent function found in this file
   % below main.)
   alhpa_complex_plot(DT, GoodIndex, alpha_complex)
   
   % Create and print the title of the alpha complex plot.  This will
   % remind the user what their selection choice was.  Title depends on
   % choice of selection scheme.
   if alpha_select_option == 1
      S = sprintf('Percentile of alpha edges kept: %0.3f',keep_percent);
   elseif alpha_select_option == 2
      S = sprintf('Kept all edges below %0.3f std. devs. from the mean.',...
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
      done = true;
   end
end% while alhpa_complex_option == 1
end % main function -- AlphaCellsSelect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       ######  #     #  ####        #     #     #     #####  #     #
%       #       ##    #  #   #       ##   ##    # #      #    ##    #
%       #       # #   #  #    #      # # # #   #   #     #    # #   #
%       ####    #  #  #  #    #      #  #  #  #######    #    #  #  #
%       #       #   # #  #    #      #     #  #     #    #    #   # #
%       #       #    ##  #   #       #     #  #     #    #    #    ##
%       ######  #     #  ####        #     #  #     #  #####  #     #
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      DEPENDENT FUNCTION     %%%%%%%%%%%%%%%%

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