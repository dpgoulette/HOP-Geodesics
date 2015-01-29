% function GeosOnly(DT, GoodIndex, maxconnect, maxclass, maxindex,...
%    GoodMaxGeodesics, pause_option)
%  COMMENT THIS!!
%
%  GoodMaxGeodesics is unused and may be unneeded.  It comes from
%  SelectGeodesics.  Consider renaming it to: num_selected_geos (or
%  something like that)
%
% Loops through the maxima and plots the geodesics connected to that
% particular max along with the persistence diagrams.
pause_option = 2;

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

% while true
%    fprintf('Pause between maxima?\n')
%    fprintf('  1) Yes\n')
%    fprintf('  2) No\n')
%    pause_option = input('(Enter 1 or 2): ');
%    if pause_option == 1 || pause_option == 2
%       break
%    else
%       fprintf('\nERROR: You must enter 1 or 2.\n')
%    end
% end
% fprintf('\n')

if pause_option == 1
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
else
   delete_option = 2;
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
points_plot = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.');
axis equal
hold on
axis tight
set(gca, 'XTick', []);
set(gca, 'YTick', []);

if tree_option == 1
   E = vertcat(maxclass.hoptree);
   X=[DT.X(E(:,1),1)';DT.X(E(:,2),1)'];
   Y=[DT.X(E(:,1),2)';DT.X(E(:,2),2)'];
   hoptree_plot = plot(X,Y,'b:');
end

maxima_plot = plot(DT.X(maxindex,1),DT.X(maxindex,2),'g.');


GeoIndex = 1;
for i=1:length(maxclass)
   %Now check to see that the max class has max neighbors. (There are odd
   %cases on the boundary of the data space where a max class is isolated.
   %These cases are pathological and can be ignored.)
   if ~isempty(maxclass(i).nbormaxid)
      % Now get the number of max neighobrs for the ith max.
      NumMaxNeighbors = length(maxclass(i).nbormaxid);
      Pmax1 = plot(DT.X(maxclass(i).max,1),DT.X(maxclass(i).max,2),'r*');
      
      % make arrays for storing plot handles; for batch deleting later
      GeoPlotHandles1 = cell(NumMaxNeighbors,1);
      
      if pause_option ==1
         pause
      end
      
      for j=GeoIndex:(GeoIndex+NumMaxNeighbors-1)
         % In this section:
         %     iterate through the geodesics of the ith max plot the
         %     geodesics connected to the ith max based on the user options
         %     input.
         
         % plot the geodesics connected to the current max.
         if plot_option == 2 %then only plot selected geodesics
            %  if j-GeoIndex+1 <= GoodMaxGeodesics(i)
            if maxconnect{j,1}(1) < maxconnect{j,1}(2) && ...
                  maxconnect{j,4} == 2
               % Plot only the geodesics that are included from both ends.
               % Plot in both the zoomed out and zoomed in plots.
               GeoPath = maxconnect{j,2};
               GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
               X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
               Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
               GeoPlotHandles1{j-GeoIndex+1} = plot(X,Y,'r-');
            end
         else % plot all geodesics
            %Plot this geodesic connected to the current max.  Plot in both
            %the zoomed out and zoomed in plots.
            GeoPath = maxconnect{j,2};
            GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
            X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
            Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
            GeoPlotHandles1{j-GeoIndex+1} = plot(X,Y,'r-');
         end
         
         if pause_option == 1
            pause
         end
      end
      
      %Update the GeoIndex
      GeoIndex = GeoIndex+NumMaxNeighbors;
      if delete_option == 1
         %delete all of the plots
         temp = vertcat(GeoPlotHandles1{:});
         delete(temp);
      end
      delete(Pmax1);
   end
end
% end% function

