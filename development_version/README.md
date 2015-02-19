HOP-Geodesics
=============
This is code is functional currently but it is not complete and is under active development (Feb 18th, 2015)
The code in this repository is the MATLAB code explained in the parent directory dpgoulette/HOP-Geodesics.
 
This code requires raw 2d or 3d data spatial data.  Two sample data sets are included in this directory: "3D_example.txt" and "FlatDataExample.txt."  You need to download all of the files in this repository and save them in a single directory (or in a directory that is in your MATLAB path).  Open the HOP_Geodesics_main.m file and read the comments on how to run that script.  The first line of code in that file loads one of the two example files and runs the model on that data set.  You can provide your own data of course.  HOP_Geodesics_main.m will run the model in an interactive session in the MATLAB command window.  Various options will be presented to you and the results will be plotted based on your choices.  The data structures that are created will be saved in your workspace.  The key results of the model are hop, maxclass and minclass (7-9 below).  Items 1-6 below are saved so that if you want to rerun the model with different settings, the code will use the saved data in the workspace.  So reruning the model is much faster than the first run.  The key saved data structures are:

| |Data variable name| explanation|
|----:|----------------------|------------------------------------------|
|1|DT |The Delaunay Triangulation object|
|2| VV | The Voronoi cell vertices |
|3| VC | The indices of the vertices surrounding each Voronoi cell |
|4| GoodIndex | The indices of the "good points" (We throw out data points on the boundary of the data space, the "good points" are the rest.) |
|5| GoodEdges | The edges in the graph that we apply HOP to. |
|6| edge_alpha | If you choose to apply HOP to the 1-skeleton of an alpha complex (which is a certain subset of the Delaunay triangulation), then the edges along with their corresponding alpha value are stored in this matrix. |
|7| hop | An array of structs that is as long as the raw data set.  This stores key information for each data point. |
|8| maxclass | An array of structs that is as long as the number of maxclasses.  The maxclasses are a partition of the data space that results from appling the HOP algorithm.  The fields in this array of structs hold key information about each maxclass. |
|9| minclass | Similar to maxclass. |

This code requires a version of MATLAB that has the TriRep class.  Specifically, you need the DelaunayTri subclass.  This is needed to create the Delaunay triangulation of the data as well as the Voronoi cells.




 
 

