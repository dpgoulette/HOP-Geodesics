# HOP-Geodesics
Code to automatically locate geometric/topolocial structures (such as clusters and filaments) in SDSS data using a modified version of the HOP algorithm.

## Brief explanation
Welcome!  You have found the code repository for HOP-Geodesics, which is code to analyze two and three dimensional point cloud data.  The code is under development and there may be substantial changes coming.  In a nutshell, this code takes raw point data and locates points of max density (locally).  Then it connects these local maxima along geodesic paths that trace out filimentary structures in the data (thereby tracing out the so-called "cosmic web").  We are developing approaches to detect 2D and 3D structures as well.  Our focus was on applying this algorithm to astronomical data like data from the Sloan Digital Sky Survey (SDSS).  We are currently preparing a paper for publication which will explain the theory.

## Somewhat longer version (incomplete)
The data that we are using consists of ~140,000 galaxies from the SDSS.  We think of these galaxies as a sample from the distribution of mass in the universe.  We are studying how to take this discrete sample and reconstruct the structures that they were sampled from.  We accomplish this by assigning a density function to the set of points and locating local maxima.  Then we connect these maxima by certain geodesic paths in the simplicial complex that outline local 1-D substructures.  Our results generalize to any coordinate data with a metric on the data space. So they could be applied in various applications where you want to understand more about the structure of a sample space.

## How to run the code
(Coming soon)


