% Change the first line in the script so that the "Data" variable holds
% whatever the name of your data matrix is.  The matrix must be 2D or 3D
% data.  

% Data = load('FlatDataExample.txt')

while true
   fprintf('\nWhat data file would you like to study?\n')
   NameString = input('Enter the name of the file in single quotes including the file extension: ');
   Data = load(NameString);
   break
end

% Prepare the raw data for HOP. Calculate the Delaunay triangulation, the
% Voroinoi diagram, throw away "bad" data on the boundary of the data
% space, etc.
[DT,VV,VC,BadDataID,DTedges,GoodIndex,GoodEdges] = HOPDataPrepare(Data);
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

fprintf('\nDo you want to plot the results? \n');
while true
   fprintf('   1) Yes.\n')
   fprintf('   2) No.\n')
   plot_option = input('Choose one of the above: ');
   if plot_option == 2
      clear plot_option
      break
   elseif plot_option == 1
      clear plot_option
      run('GeodesicPlot')
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
