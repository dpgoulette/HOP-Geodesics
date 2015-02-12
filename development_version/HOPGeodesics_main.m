% This is the main script to run the HOPGeodesics data model.  

%%%%%  MAIN SCRIPT %%%%%%

% Change the first line in the script so that the "Data" variable holds
% whatever the name of your file holding your data array is.  The matrix
% must be 2D or 3D data.

Data = load('3D_example.txt');
% Data = load('FlatDataExample.txt');

% Prepare the raw data for HOP. Calculate the Delaunay triangulation and
% the Voroinoi diagram.  Then throw away "bad" data on the boundary of the
% data space.  Also, the user has an option to calculate an alpha complex
% in various ways. 
%
% NOTE! The three possibilities for this function are to allow for the user
% to rerun the model on a large data set without having to recalculate DT,
% VV and VC.  See the help section of this file and HOPDataPrepare.

if exist('last_run_data','var') == 1
   [DT, VV, VC, GoodEdges, GoodIndex, edge_alpha] =...
      HOPDataPrepare(Data, last_run_data);
else
   [DT, VV, VC, GoodEdges, GoodIndex, edge_alpha] = HOPDataPrepare(Data);
end


%%%%%%%%%%%%%% HOP and related data structures %%%%%%%%%%%%%%%%%
% The next two functions collectively do the HOP algorithm and catalog a
% variety of related data that result from doing HOP.  The key results are
% saved in the hop, maxclass and minclass structs.

% Now create the hop data structure and do the HOP algorithm based on our
% denisity function. Locate all of the hop maxima and minima. etc.
[hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);
% Now create the maxclass and minclass structure which holds all of the
% information for each hop max/min class.
[hop, maxclass, minclass] = HOPClasses(hop,maxindex,minindex,DT);

%%%%%% Finished with the HOP section and struct creation. %%%%%%



%%%%%%%%%%%%%%% HOP Geodesics section %%%%%%%%%%%%%%%
% In this section we find the geodesic paths which connect neighboring
% maxima.  The results in this section are all stored in the maxclass
% struct.

% Calculate all of the geodesics connecting neighboring max classes.
maxclass = hopmaxconnect(DT.X,maxclass,hop);

% Select the geodesics that we will keep and which we will throw out (if
% any).
maxclass = SelectGeodesics(maxclass, hop);

% The geodesics are already the 1D connections between maxima.  Now we want
% to find where geodesics form triangles (2D structure) and tetras (3D
% structure).  Of course there will only be tetras if the data is 3D to
% begin with.
maxclass = GeoTris(maxclass);
if size(DT.X,2) == 3
   maxclass = GeoTetras(maxclass);
end

% It classifies each geodesic as either being in a tetra, triangle, or only
% an isolated geodesic.  This roughly means that the geodesic is a part of
% a 3D, 2D or 1D structure.
maxclass = geodesic_classify(maxclass);

% Finished finding hop geodesics and related structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save the important data structures that were calculated to speed up
% rerunning the model.
last_run_data = struct('DT', {DT}, 'VV', {VV}, 'VC', {VC},...
   'GoodIndex',{GoodIndex}, 'GoodEdges',{GoodEdges},...
   'edge_alpha',{edge_alpha});

% Plotting section
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
      GeodesicPlot(DT, GoodIndex, maxclass, maxindex)
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
