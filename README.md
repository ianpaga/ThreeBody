Three-Body Problem: Ballistic trajectories
====

Description: The code solves the Newtonian equations for the 3-body problem. The setup consists of two stars, the donor and the accretor, and a test particle of almost vanishing mass.
Three possible outcomes are possible: direct impact (DI), disk formation (DF), or self-accretion (SA). If the code needs more time to find either of these outcomes, a fourth case is printed as NA.

## Figures:

![trajectory](https://github.com/ianpaga/ThreeBody/assets/57350668/caffeea7-c79d-4666-b18f-0ea6ca81faa4)
<img width="814" alt="scan_mass_J" src="https://github.com/ianpaga/ThreeBody/assets/57350668/61be3abb-f87e-416f-82bc-730d30516870">

## Requirements:
- C++ compiler
- [Boost Library](https://www.boost.org/)
- OpenMP
- Python, Matplotlib, NumPy

## Compiling and running:
- Run ./chpc to compile
- Run executable *.out
- Outputs *.raw and *.sum are large files. Use the bash script datfiles.sh to slice the data into smaller *.dat files
- Plot results with python potentials.py (see examples in /plots)
