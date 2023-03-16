# I. A few hint on the ParticleAdvectorSolver

## Type of particle:

1) `Nodal`: particles sometimes have problems in 3D, because at the surface velocity may point of of mesh. 
2) `Elemental`: more robust in this way but looks bad to the eye
3) `Discontinous Galerkin`: Not exactly DG since the initial location is actually on the gaussian points.
4) `Integration Point`: No additional benefits

## A few additional keywords

`Source Time Correction = Logical True`: when estimating the given source and facing dependencies automatically chooses old time from transient simulation in correct ratio.


# II.  Simulation Setup

The goal of this test, based on the existing "Advection_Reaction_TestCase" is to test the Semi-Lagrangian solver.

## Case description

The case consists in a 3 disc-shapes that are turning with a constant orbit given by a circular velocity field. We consider two metrics:

1) the pure advection of theses shapes: there we look for methods that minimize numerical diffusion. The variable is named `Hadv`
2) an advection that is a function of these evolving shape, i.e. a source term that depends on the evolution of the particle along the path integral over time. The Variable is named `Source`. 

We define the source term in the body force as

```f90
Particle Time Integral Source = Variable Hadv
  Real matc "1.1*tx"
``` 

## Variations

We test 3 different variations of the solver:

1) Serial vs. Parallel (2 partitions) (no change in the sif but we have a "serial" mesh and a "2-partition" mesh)
2) Nodal vs. Elemental Advection (sif files labelled " Nodal/Elem")
3) Particle reinitialisation vs No particle reinitialisation (sif files labelled "Reinit/NoReinit")

This gives us a setup with 2x2x2=8 simulations. We then look at potential problems linked to each solution. Here is a recapitulative table of the simulations:

<figure>
<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/Recap_table.png" width=75% height=75%>
</center>
</figure>

 

The folder containing the ".sif" files, USFs, and the mesh  used for each simulation can be found in this repostiroy under the name "AdvReac_TestCase.REF".

# Results

## 1. No Particle Reinitialization

In this case, we set `Reinitialize particles = Logical False`.

Here, the input variable is `HPart` and the solver outputs the advected value `Hadv`. 
We also compute the `Particle Time Integral` that we define as `Source`.

### 1.1 Nodal

The input and output variables are Nodal.
 
#### 1.1.1 Serial

In Serial, `Hadv` is well advected with no loss. `ParticlePathIntegral` evolves upfront the position of the source term `Real matc "1.1*tx"`. Partially linked to the fact that the `Particle Time Integral` advects the `Source` term that depends itself on `Hadv`; but the upfront position looks larger than a1 timestep advection.

<figure>
<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/NoReinit_Nodal_Serial.gif" width=50% height=50%>
<figcaption>Fig. Advected field Hadv (left) and ParticlePathIntegral (right). The grey circles on the right show the Hadv field which is used as a source for ParticlePathIntegral.</figcaption>
</center>
</figure>

#### 1.1.2 Parallel

No problem linked to the parallelisation.

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/NoReinit_Nodal_Parallel.gif" width=50% height=50%>
</center>

### 1.2 Elemental

The input variable is nodal, the ouput variable is elemental.

#### 1.2.1 Serial

In Serial, `Hadv` is well advected with no loss. `Source` remains to zero over the entire simulation. Notice that if we change `Hadv` for `Hpart` then we have an evolving source term. 

```f90
Particle Time Integral Source = Variable Hpart
  Real matc "1.1*tx"
```

Does `Particle Time Integral Source` needs a Nodal variable as an input?

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/NoReinit_Elemental_Serial.gif" width=50% height=50%>
</center>

#### 1.2.2 Parallel

No problem lined to the parallelisation.

<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/NoReinit_Elemental_Parallel.gif" width=50% height=50%>

## 2. Particle Reinitialization


> **⚠ WARNING** : diffusion related to reinitialization is inherent to the method. This is because reinitialization requires to reinterploate the input field to the particles. 
> Diffusion should decrease with mesh size.
> Diffusion can be decreases when following particles with very long timesteps. 
> The reinitialization should be avoided if possible (no reinitialization every timestep).

 
In this case, we set `Reinitialize particles = Logical True`.
Here, the input variable is `HPart` and the solver outputs the advected value `Hadv`. We then udpate `Hpart= Equals Hadv` in the bodyforce using the `UpdateExported` solver. 

We also compute the `Particle Time Integral` that we define as `Source`. 

### 2.1 Nodal 

### 2.1.1  Serial

Here we notice a clear numerical diffusion. Each timesetp uses the new interpolated `Hpart` as a restart and we loose information. The ParticlePathIntegral seems to work great.

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Reinit_Nodal_Serial.gif" width=50% height=50%>
</center>

### 2.1.2 Parallel

Parallelization of the ParticlePathIntegral seems to have indexation problems, the field directly evolves in a "speckled" pattern.

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Reinit_Nodal_Parallel.gif" width=50% height=50%>
</center>

### 2.2 Elemental


<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Reinit_Elemental_Serial.gif" width=50% height=50%>
</center>

### 2.2.1 Serial

Important numerical diffusion linked to the reinitilization and problem with the ParticlePathIntegral that remains to zero. 

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Reinit_Elemental_Parallel.gif" width=50% height=50%>
</center>

### 2.2.2 Parallel

No particular issue linked to the parallelisation. 

<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Reinit_Elemental_Parallel.gif" width=50% height=50%>
</center>

# Restarts

Here we check how the different implementations handle restart files. We show one example of a problem linked to the restart, i.e. the `ParticlePathIntegral` that does not account for initial conditions. For that, we restart the simulation at half-time (after a rotation of 180 degree). 

The restart is executed from the cass 1.1.1 (i.e. No Reinitialisation / Nodal / Serial). Problem: the "input" variable is intialized with `Particle Time Integral = Equals Source`. `Source` is the result of the previous simulation but is not accounted for.
<center>
<img src="https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/animations/Restart_NoReinit_Nodal_Serial.gif" width=50% height=50%>
</center>

Currently, the `SUBROUTINE ParticlePathIntegral` does not initialize the particle to an inital condition, integrating this to the subroutine would be very useful.



