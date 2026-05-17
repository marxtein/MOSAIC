# MOSAIC (Magnetized Orbit Simulation and Analysis Integration Code)

## 1. Overview
MOSAIC is a comprehensive phase-space analysis and simulation framework designed for investigating energetic particle (EP) dynamics and wave-particle interactions in tokamak plasmas. 

This standalone package provides robust guiding-center / gyrokinetic test-particle orbit tracking and advanced phase-space diagnostics. Key capabilities include:
* High-Fidelity Orbit Integration: Evaluates guiding-center equations in Boozer magnetic flux coordinates using a 1st-order continuous (C^1) B-spline interpolation scheme for exact equilibrium field gradients.
* Characteristic Frequencies: Numerical extraction of poloidal bounce/transit and toroidal precession frequencies.
* Phase-Space Mapping: Construction of resonance lines and orbit accessibility boundaries in 3D Constants-of-Motion (CoM) coordinates (P_zeta, lambda, E) or alternative formats.
* Nonlinear Wave-Particle Interaction: Generation of Kinetic Poincare maps under single-mode time-dependent Alfvenic perturbations (e.g., RSAEs).
* Verification Tools: Direct comparative tools to overlay analytical resonance lines with perturbed distribution data (delta f) from gyrokinetic simulations (e.g., GTC).

This repository contains all the necessary scripts to reproduce the figures presented in the associated MOSAIC publication.

------------------------------------------------------------
## 2. Quick Start & How to Run

To reproduce the results and figures from the paper, open the project in MATLAB and navigate to the 'test_particle' directory. The primary interactive driver is:

    runme.m

The 'runme.m' script is divided into four main computational sections. Please enable and run only one section at a time by adjusting the control parameters at the top of that specific section. Each section automatically calls the core workflow script ('run_paper.m') to execute the task.

The Four Computational Tasks:
1. Test-Particle Orbit Calculation: Extracts and plots fundamental guiding-center trajectories.
2. Resonance-Line Computation: Calculates characteristic frequencies and constructs analytical phase-space resonance lines.
3. Distribution Comparison (Verification): Overlays computed resonance lines with delta f distribution data.
4. Kinetic Poincare Plots: Generates stroboscopic phase-space maps to reveal resonance island topologies and chaotic regions.

------------------------------------------------------------
## 3. Key Configuration Parameters

Inside 'runme.m', you can modify the following key parameters to customize the simulation or switch between different figure generation modes:

* ps_option : Phase-space / plotting mode selector
  - 'PLam_traj' : Trajectory projection
  - 'PLam' : P_zeta-lambda resonance mapping
  - 'MG' : Main resonance mode
  - 'GB' : Alternative mapping
  - 'Poincare' : Kinetic Poincare plots

* passing_option : Particle motion direction relative to plasma current
  - 'co-passing' or 'counter-passing'

* ant_option : Antenna / time-dependent perturbation switch
  - 'on' or 'off'

* Eperp_in : mu*B_a value for particle initialization (e.g., 20 in keV)
* PoinE : Energy bounds for Poincare plot initialization [E_down, E_up, E_prime] (e.g., [26 34.5 30] in keV)

Execution Examples:
* To generate Resonance Lines (Single Branch): Go to Section 2 of 'runme.m', set ps_option = 'PLam', passing_option = 'co-passing', ant_option = 'off', and run the section.
* To generate Kinetic Poincare Plots: Go to Section 4 of 'runme.m', set ps_option = 'Poincare', ant_option = 'on' (required for time-dependent maps), define Eperp_in and PoinE, and run the section.

------------------------------------------------------------
## 4. Output Files and Results

Upon execution, the code will automatically generate and save data and figures into the following directories:

* /save/*.mat : Stores computed numerical data arrays, including bounce/precession frequencies and phase-space grid matrices for post-processing.
* /output/*.eps : High-quality vector graphics ready for publication. 

------------------------------------------------------------
## 5. Dependencies and Requirements

* Environment: MATLAB R2023b or later.
* Toolboxes: No specialized toolboxes are required for core integration.
* External Data: The repository includes sample equilibrium data necessary to run the basic tests. Extended datasets can be supplied or generated via their respective external codes.

------------------------------------------------------------
## 6. License

This project is open-sourced under the GNU General Public License v3.0 (GPLv3). 
You are free to use, modify, and distribute this software under the terms of the GPLv3 license. See the LICENSE file in the root directory for the full legal text.

------------------------------------------------------------
## 7. Citation & Contact

If you use MOSAIC in your research or find this code helpful, please cite our paper:
Wang, Z. K., Bao, J., Zhang, W. L., et al. (2026). MOSAIC: An orbit-based phase space analysis tool for energetic particle physics in tokamak. Computer Physics Communications (Submitted/In Press).

Contact Information:
Jian Bao (Corresponding Author)
Institute of Physics, Chinese Academy of Sciences
Email: jbao@iphy.ac.cn
