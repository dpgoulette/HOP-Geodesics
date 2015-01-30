function GeosWithStepFunc2D(DT, GoodIndex, maxclass, maxindex)
%  COMMENT THIS!!
%
% Loops through the maxima and plots the geodesics connected to that
% particular max along with the persistence diagrams.


% Get the plot options from the user
fprintf('This function loops through each hop max and plots the geodesics ')
fprintf('connected to it.  It can plot all geodesics attached to each max.')
fprintf('Or it can plot only the geodesics that were selected.\n\n')
fprintf('Which do you want to see:\n')
while true
   fprintf('    1) Every geodesic.\n');
   fprintf('    2) Selected geodesics.\n')
   plot_option = input('Which of the above options do you prefer: ');
   if plot_option==1 || plot_option==2 
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end
fprintf('\n')
fprintf('The geodesics can be deleted after each iteration or they ')
fprintf('can aggregate.\n\n')
fprintf('Which do you prefer:\n')
while true
   fprintf('    1) Delete geodesics after each iteration.\n');
   fprintf('    2) Aggregate the geodesics.\n')
   delete_option = input('Which of the above options do you prefer: ');
   if delete_option == 1 || delete_option == 2 
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end
fprintf('\n')
while true
   fprintf('Would you like to see the HOP trees?\n')
   fprintf('  1) Yes\n')
   fprintf('  2) No\n')
   tree_option = input('(Enter 1 or 2): ');
   if tree_option == 1 || tree_option == 2
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end

figure
% force the subplot to be in a certain position.  The position vector is 
% [left bottom width height]
DataPlot1 = subplot('Position',[0.01 0.505 0.52 0.49]);
P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.');
axis equal
hold on
axis tight
set(DataPlot1, 'XTick', []);
set(DataPlot1, 'YTick', []);
if tree_option == 1
   E = vertcat(maxclass.hoptree);
   X=[DT.X(E(:,1),1)';DT.X(E(:,2),1)'];
   Y=[DT.X(E(:,1),2)';DT.X(E(:,2),2)'];
   hoptree_plot = plot(X,Y,'b:');
end
P2 = plot(DT.X(maxindex,1),DT.X(maxindex,2),'g.');


DataPlot2 = subplot('Position',[0.01 0.005 0.52 0.49]);
P3 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.');
axis equal
hold on
% axis tight
set(DataPlot2, 'XTick', []);
set(DataPlot2, 'YTick', []);
if tree_option == 1
   E = vertcat(maxclass.hoptree);
   X=[DT.X(E(:,1),1)';DT.X(E(:,2),1)'];
   Y=[DT.X(E(:,1),2)';DT.X(E(:,2),2)'];
   hoptree_plot = plot(X,Y,'b:');
end
P4 = plot(DT.X(maxindex,1),DT.X(maxindex,2),'g.');

BarCode1 = subplot('Position',[0.575, 0.537, 0.4, 0.4]);
BarCode2 = subplot('Position',[0.575 0.04 0.4 0.4]);

%need the max length geodesic in the dataset
temp = vertcat(maxclass.geodesics);
num_geodesics = size(temp, 1) / 2;
temp = vertcat(temp{:,2});
MaxGeoLength = max(temp(:,3));

for i=1:length(maxclass)
   % Now check to see that the max class has max neighbors. (There are odd
   % cases on the boundary of the data space where a max class is isolated.
   % These cases are pathological and can be ignored.)
   if isempty(maxclass(i).nbormaxid)
      continue
   else
      %Now get the number of max neighobrs for the ith max.
      NumMaxNeighbors = length(maxclass(i).nbormaxid);
      
      %Plot the current max in top left
      subplot(DataPlot1)
      Pmax1 = plot(DT.X(maxclass(i).max,1),DT.X(maxclass(i).max,2),'r*');
      
      % Zoom the closeup geodesics plot for the current max (bottom left plot)
      subplot(DataPlot2)
      Pmax2 = plot(DT.X(maxclass(i).max,1),DT.X(maxclass(i).max,2),'r*');
      ZoomIndex = [maxclass(i).max; maxclass(i).nbormax];
      ZoomCoords = [DT.X(ZoomIndex,1),DT.X(ZoomIndex,2)];
      MaxCoord = DT.X(maxclass(i).max,:);
      MaxXdiff = max(abs(ZoomCoords(:,1)-MaxCoord(1)));
      MaxYdiff = max(abs(ZoomCoords(:,2)-MaxCoord(2)));
      Maxdiff = max(MaxXdiff,MaxYdiff);
      scale = 1.25;
      axis([MaxCoord(1)-Maxdiff*scale, MaxCoord(1)+Maxdiff*scale,...
         MaxCoord(2)-Maxdiff*scale, MaxCoord(2)+Maxdiff*scale])

      %make arrays for storing plot handles; for batch deleting later
      GeoPlotHandles1 = cell(NumMaxNeighbors,1);
      GeoPlotHandles2 = cell(NumMaxNeighbors,1);
      Bar1 = cell(NumMaxNeighbors+1,1);
      Bar2 = cell(NumMaxNeighbors+1,1);
      
      %plot the zero bar first before the loop.  This is for the
      %persistence by geodesic LENGTH.
      subplot(BarCode1)
% %       BarX = [0; maxconnect{GeoIndex,1}(3)];
      BarX = [0; maxclass(i).geodesics{1,2}(3)];
      BarY = [0; 0];
      Bar1{1} = plot(BarX,BarY,'r-');
      hold on
      title('Persistance by length of geodesic.')
      axis([0 MaxGeoLength 0-0.3 NumMaxNeighbors+0.3])
      set(BarCode1,'XTickLabel',[])
      set(BarCode1,'Ytick',0:NumMaxNeighbors)
      
      %plot the zero bar first before the loop.  This is for the
      %persistence by gedesic length RANK.
      subplot(BarCode2)
      BarX = [0; maxclass(i).geodesics{1,3}];
      BarY = [0; 0];
      Bar2{1} = plot(BarX,BarY,'r-');
      hold on
      title('Persistance by rank order of geodesic.')
      axis([0, num_geodesics, -0.3 NumMaxNeighbors+0.3])
%       set(BarCode2,'XTickLabel',[])
      set(BarCode2,'Ytick',0:NumMaxNeighbors)
      
      pause
      
% %       for j=GeoIndex:(GeoIndex+NumMaxNeighbors-1)
      for j=1:NumMaxNeighbors
         % In this section:
         %     iterate through the geodesics of the ith max
         %     plot the geodesics connected to the ith max in order from
         %     shortest to longest (in the geodesic metric).  Two plot options
         %     here based on the PlotOption input.
         %     plot the persistence diagram as well
         
         % plot the geodesics connected to the current max.
         if plot_option == 2 %then only plot selected geodesics
            if maxclass(i).geodesics{j,4} == 2
               % Plot only the geodesics that are included from both ends.
               % Plot in both the zoomed out and zoomed in plots.
               GeoPath = maxclass(i).geodesics{j,1};
               GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
               X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
               Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
               subplot(DataPlot1)
               GeoPlotHandles1{j} = plot(X,Y,'r-');
               subplot(DataPlot2)
               GeoPlotHandles2{j} = plot(X,Y,'r-');
            end
         else % plot all geodesics
            %Plot this geodesic connected to the current max.  Plot in both
            %the zoomed out and zoomed in plots.
            GeoPath = maxclass(i).geodesics{j,1};
            GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
            X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
            Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
            subplot(DataPlot1)
            GeoPlotHandles1{j} = plot(X,Y,'r-');
            subplot(DataPlot2)
            GeoPlotHandles2{j} = plot(X,Y,'r-');
         end
         
         %%%%% Bar code 1. Persistance by the length of geodesic.  %%%%
         if j < NumMaxNeighbors
            
            BarX = [maxclass(i).geodesics{j,2}(3);
                    maxclass(i).geodesics{j+1,2}(3)];
            BarY = [j;
                    j];
         else %we are on the last geodesic for this max
            %plot the last persistance bar.
            BarX = [maxclass(i).geodesics{j,2}(3);...
                    MaxGeoLength]; 
               
            BarY = [NumMaxNeighbors;...
                    NumMaxNeighbors];
         end
         
         subplot(BarCode1)
         Bar1{j+1} = plot(BarX,BarY,'r-');
%          hold on
%          title('Persistance by length of geodesic.')
         axis([0 MaxGeoLength 0-0.3 NumMaxNeighbors+0.3])
%          set(BarCode1,'XTickLabel',[])
%          set(BarCode1,'Ytick',0:NumMaxNeighbors)
         %%%%%%%%%%%%%  end of Bar code 1 section %%%%%%%%%%%%%%%%%
         
         %%%%%%%%%%%%%  Bar code 2 section %%%%%%%%%%%%%%%%%%%%%%%% 
         % most of this code is similar to Bar code 1 section.  The only
         % difference is with the x coordinate of the changepoints in the
         % step graph.
         
% %          if NumMaxNeighbors == 1
% %             %This special case forces the entire persistance plot to be plotted all
% %             %at once.
% %             BarX = [maxconnect{j,3};...
% %                     size(maxconnect,1)/2];
% %             BarY = [1;...
% %                     1];
% %             
% %          elseif j==1% j== GeoIndex && NumMaxNeighbors > 1
% %             
% %             % then we are on the first iteration of this for loop (index j)
% %             % and the number of neighbors is more than one.  Thus not on
% %             % the last neighbor so the persistance bar will start at point
% %             % (0,0) and end at the first change point.
% %             
% %             
         if j < NumMaxNeighbors
            % Then there is more than one neigbor (hence more than one
            % geodesic).
            BarX = [maxclass(i).geodesics{j,3};
                    maxclass(i).geodesics{j+1,3}];
            BarY = [j;
                    j];
         else % j == NumMaxNeighbors. We are on the last geodesic for this max
            % plot the last persistance bar.
            BarX = [maxclass(i).geodesics{j,3};...
                    num_geodesics]; 
               
            BarY = [NumMaxNeighbors;...
                    NumMaxNeighbors];
         end
         
         
         
         
         
         
         
         subplot(BarCode2)
         Bar2{j+1} = plot(BarX,BarY,'r-');
%          hold on
%          title('Persistance by rank order of geodesic.')
         axis([0, num_geodesics, 0-0.3 NumMaxNeighbors+0.3])
%          set(BarCode2,'XTickLabel',[])
%          set(BarCode2,'Ytick',0:NumMaxNeighbors)
         
%          pause(0.5)
         pause
         
      end
      
%       pause
      
      %Update the GeoIndex
      if delete_option == 1
         %delete all of the plots
         temp = vertcat(GeoPlotHandles1{:});
         delete(temp);
         temp = vertcat(GeoPlotHandles2{:});
         delete(temp);
      end
      temp = vertcat(Bar1{:});
      delete(temp);
      temp = vertcat(Bar2{:});
      delete(temp);
      delete(Pmax1);
      delete(Pmax2);
   end
end
end% function: GeodesicPersistencePlot

