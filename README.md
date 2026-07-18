# Lid-Driven Cavity Flow LBM Solver (v2)

## Overview
This repository contains a 2D Lattice Boltzmann Method (LBM) solver written in standard C++ to simulate the classic Lid-Driven Cavity flow problem. The solver utilizes the D2Q9 velocity model and incorporates an $L_2$ error norm criterion to automatically evaluate steady-state convergence. 

This is Version 2, which includes corrected Zou/He boundary formulations for the domain corners, ensuring accurate mass conservation and stable vorticity generation.

## Features
* **Lattice Model:** D2Q9 (2 Dimensions, 9 discrete velocities).
* **Boundary Conditions:** Zou/He boundary conditions for stationary walls, the moving top lid, and all four corner nodes.
* **Convergence:** Implements an $L_2$ relative error check between macroscopic velocity fields at intervals, terminating the simulation when the error drops below $10^{-6}$.
* **Data Extraction:** Automatically pads graveyard nodes to perfectly extract the fluid domain center.

## Physical Parameters
* Reynolds Number ($Re$): 400
* Mach Number ($Ma$): 0.1
* Grid Resolution: $529 \times 529$
* Relaxation Time ($\tau$): Derived dynamically based on lattice viscosity.

## Outputs
The solver exports data in `.csv` format for easy post-processing (e.g., using Python/Matplotlib or ParaView):
* `centerline_x.csv`: Horizontal velocity profile across the geometric center.
* `centerline_y.csv`: Vertical velocity profile across the geometric center.
* `velocity_field.csv`: Full 2D domain extraction of $X$ and $Y$ macroscopic velocities.
* `Psi_lid.csv`: Computed stream function values.

## Compilation and Execution
The code is self-contained and only relies on the standard C++ library `<iostream>`, `<vector>`, `<cmath>`, and `<fstream>`.

Compile with GCC:
```bash
g++ -O3 -march=native src/main.cpp -o ldc_solver
