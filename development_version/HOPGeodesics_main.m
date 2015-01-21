% Issues:
%     1) Remove maxconnect completely.  Redo everything to use maxclass.
%     2) Fix 2D plotting so it doesn't use maxconnect, GoodMaxGeodesics, or
%        Geodesic_Tris.  Use maxclass only.
%     2a) Once 2 is complet. Remove all usage of maxconnect,
%         GoodMaxGeodesics and Geodesic_Tris.
%     3) Add an option to keep all geodesics instead of doing selection.
%     4) GeoTris is inefficient in the way it finds ALL possible triangles.
%     5) Remove false triangles in 2d plot (or remove triangles all
%        together.
%     6) Remove all usage of Geodesic_Tetras

% Change the first line in the script so that the "Data" variable holds
% whatever the name of your data matrix is.  The matrix must be 2D or 3D
% data.

% Data = load('3D_example.txt');
Data = load('FlatDataExample.txt');

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
[DT, VV, VC, BadDataID, DTedges] = HOPDataPrepare(Data);

%now remove bad edges.  They are "bad" if they contain a bad data point.
% NOTE!!  The "GoodEdges" matrix that is created here may be altered if
% the user chooses to HOP on an alpha complex instead of the full delaunay.
[R,~]=find(ismember(DTedges,BadDataID));
R=unique(R);
GoodEdges=DTedges;
GoodEdges(R,:)=[];%delete the bad edges.
GoodIndex=unique(GoodEdges);
clear R


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
   % We are done.
   clear alpha_option
else % We will HOP on an alpha complex 1-skeleton
      clear alpha_option
      [GoodEdges,one_cells] = AlphaCellsSelect(DT,GoodEdges,VV,VC,GoodIndex);
      fprintf('\nFinished selecting the alpha complex\n')
end

fprintf('Finished the initial prep of the data.\n')

% Now create the hop data structure and create the path pointer that will
% be used to do hop. Locate all of the hop maxima and minima. etc.
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
[Geodesic_Tris, maxclass] = GeoTris(maxclass);
if size(DT.X,2) == 3
   [Geodesic_Tetras, maxclass] = GeoTetras(maxclass);
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
      run('GeodesicPlot')
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
