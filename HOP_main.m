% Change this line to whatever data matrix that you want to do HOP on.  The
% matrix must be 2D or 3D data.
Data = flatpresent;

% The following block prepares the raw data and does the hop algorithm.  It
% also creates the max and min class structures.
[DT,VV,VC,BadDataID,DTedges,GoodIndex,GoodEdges] = HOPDataPrepare(Data);
[hop, maxindex, minindex] = HOPStructCreate(VV,VC,GoodEdges,GoodIndex,DT);
[hop,maxclass,minclass] = HOPClasses(hop,maxindex,minindex,DT);