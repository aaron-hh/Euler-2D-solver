# Euler-2D-solver

This directory contains all the code required to perform a 2D simulation of Euler equations. 

* **main.cpp** - This file contains the input data required for running the Euler 2-D simulation and the main function.
* **array.cpp** - This file contains the necessary functions for the data structure used in the simulation. 
* **array.H** - This file contains the function headers for `array.cpp`.
* **eos_Euler.cpp** - This file contains the ideal gas equation of state functions. 
* **eos_Euler.H** -  This file contains the function headers for `eos_Euler.cpp`.
* **system_Euler.cpp** - This file contains the functions for computing time step, defining initial conditions, defining boundary condition and iteration loops.    
* **system_Euler.H** - This file contains the function headers for `system_Euler.cpp`.
* **numerical_method.cpp** - This file contains the functions for SLIC and MUSCL Hancock solver (FORCE & HLLC).  
* **numerical_method.H** - This file contains the function headers for `numerical_method.cpp`.

