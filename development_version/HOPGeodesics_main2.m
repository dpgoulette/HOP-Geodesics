function [hop,maxclass,minclass,maxconnect,one_cells] =...
   HOPGeodesics_main2(DataFileString)
% input - the full name of the file you want to analyze as a string.  Make
%         sure to include the file extension.  e.g. 'FlatDataExample.txt'
%
% output - HOP results in our data structures.  
%
%     This function is for testing and visualizing the results for
%     different choices.  You will be asked various questions at the
%     command window.  Even when plots pop up, the program is still
%     running.  So return to the command window for more options.

if ~ischar(DataFileString)
   error('Your input is not a string. It needs to be a file string')
end
Data = load(DataFileString);

% while true
%    fprintf('\nEnter the name of the file you wish to study in single quotes\n')
%    NameString = input('including the file extension: ');
%    Data = load(NameString);
%    clear NameString
%    break
% end

% Prepare the raw data for HOP. Calculate the Delaunay triangulation, the
% Voroinoi diagram, throw away "bad" data on the boundary of the data
% space, etc.
% NOTE!!  The "GoodEdges" matrix that is returned here may be altered if
% the user chooses to HOP on an alpha complex instead of the full delaunay.
[DT,VV,VC,BadDataID,DTedges] = HOPDataPrepare(Data);

quit_option = 0;
redo_complex = 1;
replot_option = 1;
while quit_option == 0
   while redo_complex == 1
      %now remove bad edges.  They are "bad" if they contain a bad data point.
      [R,~]=find(ismember(DTedges,BadDataID));
      R=unique(R);
      GoodEdges=DTedges;
      GoodEdges(R,:)=[];%delete the bad edges.
      
      GoodIndex=unique(GoodEdges);
      
      % User selects whether they want to hop on the complete Delaunay or an
      % alpha complex (subset of Delaunay).
      fprintf('Would you like to HOP on the full Delaunay 1-skeleton, or the\n')
      fprintf('1-skeleton of an alpha complex?\n')
      while true
         fprintf('   1) HOP on full Delaunay.\n')
         fprintf('   2) HOP on alpha complex.\n')
         alpha_option = input('Choose one of the above: ');
         if alpha_option == 2 || alpha_option == 1
            break
         else
            fprintf('ERROR! You must enter 1 or 2.\n\n')
         end
      end
      fprintf('\n')
      if alpha_option == 1
         % We are done prepping.
         fprintf('Finished the initial prep of the data.\n')
         one_cells = DTedges;
      else % HOP on an alpha complex 1-skeleton
         if size(Data,2) == 2
            % then Data is 2D
            [GoodEdges,one_cells] =...
                     AlphaCellsSelect2D(DT,GoodEdges,VV,VC,GoodIndex);
            fprintf('\nFinished selecting the alpha complex\n')
         else
            % Data is 3D
            GoodEdges = AlphaCellsSelect3D(DT.GoodEdges,VV,VC,GoodIndex);
            fprintf('\nFinished selecting the alpha complex\n')
         end
      end
      
      % Now create the hop data structure and create the path pointer that will
      % be used to do hop. Locate all of the hop maxima and minima. etc.
      [hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);
      
      % Now create the maxclass and minclass structure which holds all of the
      % information for each hop max/min class.
      [hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT);
      
      % Geodesic triangle function should go here.  We need maxclass to do it.
      % The function should locate all maxclass triangles and create a geodesic
      % triangle structure.
      
      % Calculate all of the geodesics connecting neighboring max classes.  Store
      % the results in the maxclass cell.
      maxconnect = hopmaxconnect(DT.X,maxclass,hop);
      
      % Select the geodesics that we will keep and which we will throw out.
      [GoodMaxGeodesics, maxconnect] = SelectGeodesics(maxconnect,maxclass);
      
      Geodesic_Tris = GeoTris(maxclass,maxindex);
      
      fprintf('\nDo you want to plot the results? \n');
      while true
         fprintf('   1) Yes.\n')
         fprintf('   2) No.\n')
         plot_option = input('Choose one of the above: ');
         if plot_option == 2 || plot_option == 1
            break
         else
            fprintf('ERROR! You must enter 1 or 2.\n\n')
         end
      end
      if plot_option == 1
         while replot_option == 1
            run('GeodesicPlot')
            fprintf('\nDo you want to replot with different options?\n')
            while true
               fprintf('   1) Yes.\n')
               fprintf('   2) No.\n')
               replot_option = input('Choose one of the above: ');
               if replot_option == 2 || replot_option == 1
                  break
               else
                  fprintf('ERROR! You must enter 1 or 2.\n\n')
               end
            end
            if replot_option == 2
               break
            end
         end
      end
      fprintf('\nWhat do you want to do?\n')
      while true
         fprintf('   1) Redo the complex.\n')
         fprintf('   2) Quit and return the final results.\n')
         Choice = input('Choose one of the above: ');
         if Choice == 2 || Choice == 1
            break
         else
            fprintf('ERROR! You must enter 1 or 2.\n\n')
         end
      end
      if Choice == 2
         redo_complex = 0;
         quit_option = 1;
      end
   end
end
fprintf('\n\nIf you want to close all of the figures you just plotted ')
fprintf('type this at command line:\n')
fprintf('\n   >> close all\n\n')
end% main
