# Supersonic-Flow-over-a-Flat-Plate
Solution of supersonic flow over a flat plate based on the full Navier-Stokes equations. Following Anderson's CFD chapter 10.

The setup of the problem can be found in the setup.py file. You can change some things, such as the time steps and tolerance
for convergence in main. 

The meat of the code is in the MacCormack.py file, whcih contains all the functions for the program to integrate the PDEs. 

There are some minor discrepancies in the graphs generated w.r.t to those in Anderson, but I think the discrepancy is small enough to be neglected.