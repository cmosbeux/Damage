The goal of this test, based on the existing "Advection_Reaction_TestCase" is to test the Semi-Lagrangian solver. 


# Simulation Setup

# Case description

The case consists in a 3 disc-shapes that are turning with a constant orbit given by a circular velocity field. We consider two metrics:

	1) the pure advection of theses shapes: there we look for methods that minimize numerical diffusion 
	2) a source term that is a function of these evolving shape, i.e. a source term that depends on the evolution of the particle along the path integral over time. We define it in the body force as
	```f90
	Particle Time Integral Source = Variable Hadv
	  Real matc "1.1*tx"
	``` 


## Variations

We test 3 different variations of the solver:

	1) Serial vs. Parallel (2 partitions)
	2) Nodal vs. Elemental Advection
	3) Particle reinitialisation vs No particle reinitialisation


This gives us a setup with 2x2x2=8 simulations. We then look at potential problems linked to each solution. 

The ".sif" file used for each simulation can be found in this repostiroy.

