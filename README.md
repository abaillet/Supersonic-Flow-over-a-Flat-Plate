# Supersonic-Flow-over-a-Flat-Plate
Solution (in Python) of supersonic flow over a flat plate based on the full Navier-Stokes equations. Following Anderson's CFD chapter 10.

The setup of the problem can be found in the `setup.py` file. You can change some things, such as the time steps and tolerance
for convergence in `main`. 

The meat of the code is in the `MacCormack.py` file, which contains all the functions for the program to integrate the PDEs. Almost the entire code (except for the time loop) is vectorized and I think the speed is ok, the solution for the problem in Anderson takes about 10 s for me.

There are some minor discrepancies in the graphs generated w.r.t to those in Anderson, but I think the discrepancy is small enough to be neglected. Also, I didn't code up the M=25 at 200,000 ft example, but you should be able to implement this easily in `setup.py` (or any other freestream conditions you want).

Note that I used `tqdm` to display a progress bar. If you don't want to install this package I would suggest commenting out line 6 in `main.py` and changing line 38 to:

for i in (range(timeSteps)):
  ...
