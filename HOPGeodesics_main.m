% Change the first line in the script so that the "Data" variable holds
% whatever the name of your data matrix is.  The matrix must be 2D or 3D
% data.
Data = FlatDataExample;

% The following block prepares the raw data and does the hop algorithm.  It
% also creates the max and min class structures.
[DT,VV,VC,BadDataID,DTedges,GoodIndex,GoodEdges] = HOPDataPrepare(Data);
[hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);
[hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT);
maxconnect = hopmaxconnect(DT.X,maxclass,hop);
GoodMaxGeodesics = SelectGeodesics(maxconnect,maxclass);
