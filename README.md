# compressible_acoustic

Modelling one-dimensional acoustic waves with a Hamiltonian finite volume scheme.
Boundary conditions concerned are all boundaries are solid walls.

There are two scripts:
compressible_acoustic_midpoint.m , compressible_acoustic_stormer_verlet.m

which refer to the different timesteppers that have been noted in the documentation 
compressible_acoustic.pdf

# Instructions

Make the folder 'MATLAB' the working directory in MATLAB, typing the name of the script will then run it

# NB
To be added: separate function for the initial condition to ensure consistency between scripts and errors.

# Firedrake folder
Adding Firedrake script to solve one- to three-dimensional compressible waves.
Testing will be done for one-dimensional waves.
Documentation to follow.

Currently Firedrake has been tested for constant basis and linear basis with 16 elements for theta = 0.5. 

To be added error analysis, tests with different theta.