HOP-Geodesics
=============
 
MATLAB code to automatically locate clusters and filaments in SDSS data using a modified version of the HOP algorithm.
 
This code requires raw 2d or 3d data.  You need to download all of the files in this repository.  Open the HOP_main.m file in matlab and change the variable "YourData" to be whatever matrix of points you wish to use (see the comments in HOP_main.m for more).
 
This code requires a version of MATLAB that has the TriRep class.  Specifically, you need the DelaunayTri subclass.  This is needed to create the Delaunay triangulation of the data as well as the Voronoi cells.
 
 

