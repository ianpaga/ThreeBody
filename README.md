Three-Body Problem: Ballistic trajectories
====

Description: The code solves the Newtonian equations for the 3-body problem. The setup consists of two stars, the donor and the accretor, and a test particle of almost vanishing mass.
Three possible outcomes are possible: direct impact (DI), disk formation (DF), or self-accretion (SA). If the code needs more time to find either of these outcomes, a fourth case is printed as NA.

## Figures:

![trajectory](https://github.com/ianpaga/ThreeBody/assets/57350668/caffeea7-c79d-4666-b18f-0ea6ca81faa4)
<img width="814" alt="scan_mass_J" src="https://github.com/ianpaga/ThreeBody/assets/57350668/61be3abb-f87e-416f-82bc-730d30516870">

## Requirements:
- gfortran compiler
- [yt project](https://yt-project.org/doc/index.html) for plots
- Python, Matplotlib, NumPy

## Compiling and running:
- make clean
- make
- ./ThreeBody

  Running this Fortran code will generate files in /data. 

  One can adjust the parameters num_of_case_lim and plot_interval in Global.f90 to increase the number of systems (possible trajectories) and the number of points for plotting each trajectory, respectively. To generate the scatter plots as shown above one needs a very large number of points to scan, roughly 1 million. The default values are plot_interval = 1 and num_of_case_lim = 10 to run in a few seconds in a local machine.

# Plotting:
- Run Three Body.ipynb (jupyter notebook)
- This notebook will grab the data files in /data to make plots

