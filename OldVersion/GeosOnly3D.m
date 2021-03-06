function GeosOnly3D(DT, GoodIndex, maxclass, maxindex, pause_option, Geodesic_Tris)
%  COMMENT THIS!!
%
%  GoodMaxGeodesics is unused and may be unneeded.  It comes from
%  SelectGeodesics.  Consider renaming it to: num_selected_geos (or
%  something like that)
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

%%%%%%%%%%%%%%%%%%%%%%%%%               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%               %%%%%%%%%%%%%%%%%%%%%

figure
points_plot3 = plot3(DT.X(GoodIndex,1), DT.X(GoodIndex,2),...
   DT.X(GoodIndex,3),'k.');
axis equal
hold on
axis tight
set(gca, 'XTick', []);
set(gca, 'YTick', []);

if tree_option == 1
   E = vertcat(maxclass.hoptree);
   X=[DT.X(E(:,1),1)';DT.X(E(:,2),1)'];
   Y=[DT.X(E(:,1),2)';DT.X(E(:,2),2)'];
   Z=[DT.X(E(:,1),3)';DT.X(E(:,2),3)'];
   hoptree_plot = plot3(X,Y,Z,'b:');
end

maxima_plot = plot3(DT.X(maxindex,1),DT.X(maxindex,2),DT.X(maxindex,3),'g.');



for i=1:length(maxclass)
   %Now check to see that the max class has max neighbors. (There are odd
   %cases on the boundary of the data space where a max class is isolated.
   %These cases are pathological and can be ignored.)
   if ~isempty(maxclass(i).nbormaxid)
      % Now get the number of max neighobrs for the ith max.
      if pause_option == 1
         % Emphasize the current max while paused
         Pmax1 = plot(DT.X(maxclass(i).max,1),DT.X(maxclass(i).max,2),'r*');
         pause
      end
      
      % make arrays for storing plot handles; for batch deleting later
      if delete_option == 1
         GeoPlotHandles1 = cell(size(maxclass(i).nbormaxid,1),1);
      end
            
      for j=1:size(maxclass(i).geodesics,1)
         % In this section:
         %     iterate through the geodesics of the ith max plot the
         %     geodesics connected to the ith max based on the user options
         %     input.
         
         % plot the geodesics connected to the current max.
         if plot_option == 2 %then only plot selected geodesics
            %  if j-GeoIndex+1 <= GoodMaxGeodesics(i)
            if maxclass(i).geodesics{j,1}(1) < maxclass(i).geodesics{j,1}(2) && ...
                  maxclass(i).geodesics{j,4} == 2
               
               GeoPath = maxclass(i).geodesics{j,1};
               GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
               X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
               Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
               Z=[DT.X(GeoEdges(:,1),3)';DT.X(GeoEdges(:,2),3)'];
               
               % TEST ONLY IF STATEMENT TO ONLY PLOT 1D
               color_style = color_select(maxclass,i,j);
               GeoPlotHandles1{j} = plot3(X,Y,Z,color_style);
            end
         else % plot all geodesics
            %Plot this geodesic connected to the current max.
            GeoPath = maxclass(i).geodesics{j,1};
            GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
            X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
            Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
            Z=[DT.X(GeoEdges(:,1),3)';DT.X(GeoEdges(:,2),3)'];
            color_style = color_select(maxclass,i,j);
            GeoPlotHandles1{j} = plot3(X,Y,Z,color_style);
            
         end
         
         if pause_option == 1
            pause
         end
      end
      
      if delete_option == 1
         %delete all of the plots
         temp = vertcat(GeoPlotHandles1{:});
         delete(temp);
      end
      if pause_option == 1
         delete(Pmax1);
      end
   end
end


% for a = 1:size(Geodesic_Tris,1)
%    for b = 1:size(Geodesic_Tris{a,3},1)
%       if a == Geodesic_Tris{a,1}(b,1)
%          %Then we haven't plotted it yet.
%          if plot_option == 2
%             % Here we only plot triangles created by SELECTED geodesics.
%             if Geodesic_Tris{a,3}{b,2} == 1
%                geo_tri_temp = Geodesic_Tris{a,3}{b,1};
%                temp = [geo_tri_temp', [(geo_tri_temp(2:end))';geo_tri_temp(1)]];
%                X=[DT.X(temp(:,1),1)';DT.X(temp(:,2),1)'];
%                Y=[DT.X(temp(:,1),2)';DT.X(temp(:,2),2)'];
%                Z=[DT.X(temp(:,1),3)';DT.X(temp(:,2),3)'];
%                plot3(X,Y,Z,'r-');
%                %                x=DT.X(geo_tri_temp,1)';
%                %                y=DT.X(geo_tri_temp,2)';
%                %                patch(x,y,'magenta','FaceAlpha',.7,'EdgeColor','magenta');
%             end
%          else
%             % Here we create ALL geodesic triangs.  This is used with alpha
%             % complexes.
%             geo_tri_temp = Geodesic_Tris{a,3}{b,1};
%             temp = [geo_tri_temp', [(geo_tri_temp(2:end))';geo_tri_temp(1)]];
%             X=[DT.X(temp(:,1),1)';DT.X(temp(:,2),1)'];
%             Y=[DT.X(temp(:,1),2)';DT.X(temp(:,2),2)'];
%             Z=[DT.X(temp(:,1),3)';DT.X(temp(:,2),3)'];
%             plot3(X,Y,Z,'b-');
%             %             x=DT.X(geo_tri_temp,1)';
%             %             y=DT.X(geo_tri_temp,2)';
%             %             patch(x,y,'magenta','FaceAlpha',.7,'EdgeColor','magenta');
%          end
%       end
%    end
% end
end% main function

function color_style = color_select(maxclass,i,j)
%
if maxclass(i).geodesics{j,5} == 1
   color_style = 'm-';
elseif maxclass(i).geodesics{j,5} == 2
   color_style = 'm-';
elseif maxclass(i).geodesics{j,5} == 3
   color_style = 'm-';
end
end %color_select

