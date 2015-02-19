# HOP-Geodesics
Code to automatically locate geometric/topolocial structures (such as clusters and filaments) in SDSS data using a modified version of the HOP algorithm.

## Brief explanation
Welcome!  You have found the code repository for HOP-Geodesics, which is code to analyze two and three dimensional point cloud data.  The code is under development and there may be substantial changes coming.  In a nutshell, this code takes raw point data and locates points of max density (locally).  Then it connects these local maxima along geodesic paths that trace out filimentary structures in the data (thereby tracing out the so-called "cosmic web").  We are developing approaches to detect 2D and 3D structures using these geodesics.  Our focus was on applying this algorithm to astronomical data.  We are developing and testing with a subset of galaxy positions from the Sloan Digital Sky Survey (SDSS), but any point cloud data set could be used (e.g. Millenium Simulation).  We are currently preparing a paper for publication which will explain the theory underlying the code in this repo.

## Somewhat longer explanation (woefully incomplete at the moment...)
The data that we are using in our study is a subset of the SDSS.  Our data consists of ~140,000 galaxy positions.  We only use galaxy position in 3D; we don't consider mass, velocity, etc.  For the sake of our work we think of these galaxies as a sample from the distribution of mass in the universe.  We are studying how to take a discrete sample (of galaxies) and reconstruct the geometric/topological structures that they were sampled from (the mass in the universe).  If you remove the parenthetical additions to the previous sentence you have a more general statement of our mathematical interest.  It just so happens that we were inspired by a problem posed to our research group by Dr. Jeff Scargle at NASA-Ames.  Although I can't explain the details here, it turns out that our results hold for any discrete n-dimensional point set paired with any metric on the space that the sample is from.  So our algorithms may in the future be used to understand the structure of the underlying sample space for observed data.

## How to run the code
As of today (Feb 18th, 2015) the current development version of the code is in the "development_version" directory.  See the readme in that directory.


