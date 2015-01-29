% This is the main script to run the HOPGeodesics data model.  We first
% list the current issues that need to be resolved.

% Issues:
%
%     - restructure the main and dependencies for clarity.  Push the user
%     questions into the functions where possible.  Have the main only save
%     the main data needed for data analysis. 
%
%     - hopmaxconnect still uses
%     maxconnect internally.  It works but the code is not as simple as it
%     could be.  We possibly want to redo it to use only maxclass.
%
%     - GeoTris still creates a Geodesic_Tris cell that
%     we are no longer using.  (GeoTris is not returning this cell
%     currently.)  Consider removing this from the script.  Do everything
%     with the maxclass struct.
%
%     - GeoTetras still creates a cell for the tetras.  Same issue as the
%     previous issue.
%
%     - Need comments and cleanup throughout
%           hopmaxconnect
%           GeodesicPlot (and dependencies)
%
%     - Fix the inner for loop in HOPStructCreate.  Too many hidden
%     breakpoints. Should be a while loop.  Easier to read.
%
%     - Fix struct preallocation: hop, maxclass, minclass (where needed).
%
%     - Related to last issue: fix the assignment to the fields in
%     HOPClasses.  Unneeded preallocation (preassignment to a fields in not
%     needed if the array of structs is already initialized).
%
%     - HOPClasses still finds the boundary and interior edges but this is
%     still not well defined.  Consider revisiting and clarifying.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%%%%%  MAIN SCRIPT %%%%%%

% Change the first line in the script so that the "Data" variable holds
% whatever the name of your file holding your data array is.  The matrix
% must be 2D or 3D data.

% Data = load('3D_example.txt');
Data = load('FlatDataExample.txt');

% Prepare the raw data for HOP. Calculate the Delaunay triangulation and the
% Voroinoi diagram.  Then throw away "bad" data on the boundary of the data
% space.  Also, the user has an option to calculate an alpha complex in
% various ways.
[DT, VV, VC, GoodEdges, GoodIndex] = HOPDataPrepare(Data);

% Now create the hop data structure and do the HOP algorithm based on our
% denisity function. Locate all of the hop maxima and minima. etc.
[hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);

% Now create the maxclass and minclass structure which holds all of the
% information for each hop max/min class.
[hop, maxclass, minclass] = HOPClasses(hop,maxindex,minindex,DT);

% Calculate all of the geodesics connecting neighboring max classes.  Store
% the results in the maxclass cell.
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

% The following is a temporary function.  It classifies each geodesic as
% either being in a tetra, triangle, or only an isolated geodesic.  This
% roughly means that the geodesic is a part of a 3D, 2D or 1D structure.
maxclass = geodesic_classify(maxclass);

% Plotting the results
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
