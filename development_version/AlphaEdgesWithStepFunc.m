function AlphaEdgesWithStepFunc(DT, GoodIndex, cells1)
%  COMMENT THIS!!
%
% - Express each edge both ways for plotting and selection.
% - plot in 4 panes using template below (just like geodesics)
% - plot raw data and edges of the current point in top left
% - plot zoomed in in bottom left
% - bar code 1: Don't do anything now
% - bar code 2: rank order of alpha

% need the max length geodesic in the dataset
MaxAlpha = max(cells1(:,3));

% We need each edge in cells1 expressed both ways.
temp = [fliplr(cells1(:, [1 2])), cells1(:,3)];
double_cells1 = sortrows([cells1; temp],1);

% Now we sort the edges is double_cells1 so that the edges attached to each
% good data point is in order from smallest alpha to largest alpha.
% We store the number of neighbors each point has in num_neighbors.
% NOTE! the ith entry in num_neighbors corresponds to the ith entry in
% GoodIndex, not DT.X.
num_neighbors = zeros(size(GoodIndex));
edge_index = 1;
temp = double_cells1(:,1);
for i = 1:length(GoodIndex)-1
   num_neighbors(i) = find(temp(edge_index:end) ~= temp(edge_index), 1) - 1;
   [~,sortID] = sort(double_cells1(edge_index : ...
      edge_index + num_neighbors(i) - 1, 3));
   sortID = sortID + edge_index - 1;
   double_cells1(edge_index : edge_index + num_neighbors(i) - 1,:) = ...
      double_cells1(sortID,:);
   edge_index = edge_index + num_neighbors(i);
end
% The last entry in GoodIndex has to be outside the loop.
num_neighbors(end) = length(temp) + 1 - edge_index;
[~,sortID] = sort(double_cells1(edge_index :...
   edge_index + num_neighbors(end) - 1, 3));
sortID = sortID + edge_index - 1;
double_cells1(edge_index : end,:) = double_cells1(sortID,:);

% Add a fourth column to double_cells1 to hold the rank of alpha for each
% edge.  The rank is in order from smallest alpha to largest.
temp = [double_cells1(:,3), (1:size(double_cells1, 1))',...
   zeros(size(double_cells1(:,3)))];
temp = sortrows(temp,1);
temp(1:2:end-1,3) = (1: size(double_cells1, 1)/2)';
temp(2:2:end,3) = (1: size(double_cells1, 1)/2)';
temp = sortrows(temp,2);
double_cells1(:,4) = temp(:,3);


%%%%%%%%%%%%% This section is disabled for now but saved in case we use it.
% To simplify the code later we set plot_option to 1 so that we are
% plotting every edge connected to each point.
plot_option = 1;
% Get the plot options from the user
% fprintf('This function loops through each hop max and plots the geodesics ')
% fprintf('connected to it.  It can plot all geodesics attached to each max.')
% fprintf('Or it can plot only the geodesics that were selected.\n\n')
% fprintf('Which do you want to see:\n')
% while true
%    fprintf('    1) Every geodesic.\n');
%    fprintf('    2) Selected geodesics.\n')
%    plot_option = input('Which of the above options do you prefer: ');
%    if plot_option==1 || plot_option==2
%       break
%    else
%       fprintf('\nERROR: You must enter 1 or 2.\n')
%    end
% end
fprintf('\n')
fprintf('The edges can be deleted after each iteration or they ')
fprintf('can aggregate.\n\n')
fprintf('Which do you prefer:\n')
while true
   fprintf('    1) Delete edges after each iteration.\n');
   fprintf('    2) Aggregate the edges.\n')
   delete_option = input('Which of the above options do you prefer: ');
   if delete_option == 1 || delete_option == 2
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%% MAIN PLOTTING SECTION %%%%%%%%%%%%%
figure(1)
% force the subplot to be in a certain position.  The position vector is
% [left bottom width height]
DataPlot1 = subplot('Position',[0.01 0.505 0.52 0.49]);
P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
axis equal
hold on
axis tight
set(DataPlot1, 'XTick', []);
set(DataPlot1, 'YTick', []);

DataPlot2 = subplot('Position',[0.01 0.005 0.52 0.49]);
P3 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.'); % raw data
axis equal
hold on
% axis tight
set(DataPlot2, 'XTick', []);
set(DataPlot2, 'YTick', []);

BarCode1 = subplot('Position',[0.575, 0.537, 0.4, 0.4]);
BarCode2 = subplot('Position',[0.575 0.04 0.4 0.4]);

edge_index = 1;
for i=1:length(GoodIndex)
      
   %Plot the current point in top left
   subplot(DataPlot1)
   Point_plot_1 = plot(DT.X(GoodIndex(i),1),DT.X(GoodIndex(i),2),'r*');
   axis tight
   
   
   % Zoom the plot for both plots
   subplot(DataPlot2)
   Point_plot_2 = plot(DT.X(GoodIndex(i),1),DT.X(GoodIndex(i),2),'r*');
   ZoomIndex = [GoodIndex(i);...
      double_cells1(edge_index:edge_index - 1 + num_neighbors(i), 2)];
   ZoomCoords = [DT.X(ZoomIndex,1),DT.X(ZoomIndex,2)];
   PointCoord = DT.X(GoodIndex(i),:);
   MaxXdiff = max(abs(ZoomCoords(:,1)-PointCoord(1)));
   MaxYdiff = max(abs(ZoomCoords(:,2)-PointCoord(2)));
   Maxdiff = max(MaxXdiff,MaxYdiff);
   scale = 1.25;
   axis([PointCoord(1)-Maxdiff*scale, PointCoord(1)+Maxdiff*scale,...
      PointCoord(2)-Maxdiff*scale, PointCoord(2)+Maxdiff*scale])
   
   pause
   
   subplot(DataPlot1)
   scale = 7;
   axis([PointCoord(1)-Maxdiff*scale, PointCoord(1)+Maxdiff*scale,...
      PointCoord(2)-Maxdiff*scale, PointCoord(2)+Maxdiff*scale])
   
   %make arrays for storing plot handles; for batch deleting later
   EdgesPlotHandles1 = cell(num_neighbors(i),1);
   EdgesPlotHandles2 = cell(num_neighbors(i),1);
   Bar1 = cell(num_neighbors(i)+1,1);
   Bar2 = cell(num_neighbors(i)+1,1);
   
%    %plot the zero bar first before the loop.  This is for the
%    %persistence by geodesic LENGTH.
%    subplot(BarCode1)
%    BarX = [0; double_cells1(edge_index,3)];
%    BarY = [0; 0];
%    Bar1{1} = plot(BarX,BarY,'r-');
%    hold on
%    title('Persistance by value of alpha.')
%    axis([0 MaxAlpha -0.3 num_neighbors+0.3])
%    set(BarCode1,'XTickLabel',[])
%    set(BarCode1,'Ytick',0:num_neighbors)
   
   %plot the zero bar first before the loop.  This is for the
   %persistence by edge alpha RANK.
   subplot(BarCode2)
   BarX = [0; double_cells1(edge_index,4)];
   BarY = [0; 0];
   Bar2{1} = plot(BarX,BarY,'r-');
   hold on
   title('Persistance by rank order of edge alpha.')
   axis([0, size(double_cells1,1)/2, -0.3, num_neighbors(i)+0.3])
   %       set(BarCode2,'XTickLabel',[])
   set(BarCode2,'Ytick',0:num_neighbors)
   
   pause
   
   for j=edge_index:(edge_index+num_neighbors(i)-1)
      % In this section:
      %     iterate through the eges of the ith data point
      %     plot the edges connected to the ith point in order from
      %     smallest to largest alpha.  Two plot options
      %     here based on the PlotOption input.
      %     plot the persistence diagram as well
      
      % plot the edges connected to the current max.
      if plot_option == 2 %then only plot selected geodesics
         %            if j-edge_index+1 <= GoodMaxGeodesics(i)
%          if maxconnect{j,1}(1) < maxconnect{j,1}(2) && ...
%                maxconnect{j,4} == 2
%             % Plot only the geodesics that are included from both ends.
%             % Plot in both the zoomed out and zoomed in plots.
            Edge = double_cells1(j,[1,2]);
            X=[DT.X(Edge(1,1),1)';DT.X(Edge(1,2),1)'];
            Y=[DT.X(Edge(1,1),2)';DT.X(Edge(1,2),2)'];
            subplot(DataPlot1)
            EdgesPlotHandles1{j-edge_index+1} = plot(X,Y,'r-');
            subplot(DataPlot2)
            EdgesPlotHandles2{j-edge_index+1} = plot(X,Y,'r-');
%          end
      else % plot all edges
         %Plot this edge connected to the current point.  Plot in both
         %the zoomed out and zoomed in plots.
         Edge = double_cells1(j,[1,2]);
         X=[DT.X(Edge(1,1),1)';DT.X(Edge(1,2),1)'];
         Y=[DT.X(Edge(1,1),2)';DT.X(Edge(1,2),2)'];
         subplot(DataPlot1)
         EdgesPlotHandles1{j-edge_index+1} = plot(X,Y,'r-');
         subplot(DataPlot2)
         EdgesPlotHandles2{j-edge_index+1} = plot(X,Y,'r-');
      end
      
%       %%%%% Bar code 1. Persistance by the length of geodesic.  %%%%
%       if num_neighbors == 1
%          %This case forces the entire persistance plot to be plotted all
%          %at once.
%          BarX = [maxconnect{j,1}(3);...
%             MaxAlpha];
%          BarY = [1;...
%             1];
%          
%       elseif j==edge_index% j== edge_index && NumMaxNeighbors > 1
%          
%          % then we are on the first iteration of this for loop (index j)
%          % and the number of neighbors is more than one.  Thus not on
%          % the last neighbor so the persistance bar will start at point
%          % (0,0) and end at the first change point.
%          
%          BarX = [maxconnect{j,1}(3);
%             maxconnect{j+1,1}(3)];
%          BarY = [j-edge_index+1;
%             j-edge_index+1];
%       elseif j < edge_index+num_neighbors-1
%          %j is greater than edge_index but not on the last geodesic yet
%          BarX = [maxconnect{j,1}(3);
%             maxconnect{j+1,1}(3)];
%          BarY = [j-edge_index+1;
%             j-edge_index+1];
%       else %we are on the last geodesic for this max
%          %plot the last persistance bar.
%          BarX = [maxconnect{j,1}(3);...
%             MaxAlpha];
%          
%          BarY = [num_neighbors;...
%             num_neighbors];
%       end
%       
%       subplot(BarCode1)
%       Bar1{j-edge_index+2} = plot(BarX,BarY,'r-');
%       %          hold on
%       %          title('Persistance by length of geodesic.')
%       axis([0 MaxAlpha 0-0.3 num_neighbors+0.3])
%       %          set(BarCode1,'XTickLabel',[])
%       %          set(BarCode1,'Ytick',0:NumMaxNeighbors)
%       %%%%%%%%%%%%%  end of Bar code 1 section %%%%%%%%%%%%%%%%%
      
      %%%%%%%%%%%%%  Bar code 2 section %%%%%%%%%%%%%%%%%%%%%%%%
      % most of this code is similar to Bar code 1 section.  The only
      % difference is with the x coordinate of the changepoints in the
      % step graph.
      if num_neighbors == 1
         %This case forces the entire persistance plot to be plotted all
         %at once.
         BarX = [double_cells1(j,4);...
                 size(double_cells1,1)/2];
         BarY = [1;...
                 1];
         
%       elseif j==edge_index% j== edge_index && num_neighbors > 1
%          
%          % then we are on the first iteration of this for loop (index j)
%          % and the number of neighbors is more than one.  Thus not on
%          % the last neighbor so the persistance bar will start at point
%          % (0,0) and end at the first change point.
%          
%          BarX = [double_cells1(j,4);
%             maxconnect{j+1,3}];
%          BarY = [j-edge_index+1;
%             j-edge_index+1];
      elseif j < edge_index+num_neighbors(i)-1
         %we are not on the last edge yet
         BarX = [double_cells1(j,4);
            double_cells1(j+1,4)];
         BarY = [j-edge_index+1;
            j-edge_index+1];
      else %we are on the last geodesic for this max
         %plot the last persistance bar.
         BarX = [double_cells1(j,4);...
            size(double_cells1,1)/2];
         
         BarY = [num_neighbors(i);...
            num_neighbors(i)];
      end
      subplot(BarCode2)
      Bar2{j-edge_index+2} = plot(BarX,BarY,'r-');
      %          hold on
      %          title('Persistance by rank order of geodesic.')
      axis([0, size(double_cells1,1)/2, -0.3 num_neighbors(i)+0.3])
      %          set(BarCode2,'XTickLabel',[])
      %          set(BarCode2,'Ytick',0:NumMaxNeighbors)
      
      %          pause(0.5)
      pause
      
   end
   
   %       pause
   
   %Update the edge_index
   edge_index = edge_index+num_neighbors(i);
   if delete_option == 1
      %delete all of the plots
      temp = vertcat(EdgesPlotHandles1{:});
      delete(temp);
      temp = vertcat(EdgesPlotHandles2{:});
      delete(temp);
   end
   temp = vertcat(Bar1{:});
   delete(temp);
   temp = vertcat(Bar2{:});
   delete(temp);
   delete(Point_plot_1);
   delete(Point_plot_2);
   
end
end% function: GeodesicPersistencePlot

