% Change the first line in the script so that the "Data" variable holds
% whatever the name of your data matrix is.  The matrix must be 2D or 3D
% data.  
load('FlatDataExample.txt')
Data = FlatDataExample;

% The following block prepares the raw data and does the hop algorithm.  It
% also creates the max and min class structures.
[DT,VV,VC,BadDataID,DTedges,GoodIndex,GoodEdges] = HOPDataPrepare(Data);
[hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);
[hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT);
maxconnect = hopmaxconnect(DT.X,maxclass,hop);
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
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
