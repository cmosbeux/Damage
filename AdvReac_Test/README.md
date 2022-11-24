

# Simulation Setup

The goal of this test, based on the existing "Advection_Reaction_TestCase" is to test the Semi-Lagrangian solver. 

## Case description

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

# Results

## 1. No Particle Reinitialization

In this case, we set `Reinitialize particles = Logical False`.

Here, the input variable is `HPart` and the solver outputs the advected value `Hadv`. 
We also compute the `Particle Time Integral` that we define as `Source`.

### 1.1 Nodal

The input and output variables are Nodal.
 
#### 1.1.1 Serial

In Serial, `Hadv` is well advected with no loss. `Source` evolves but with a surprising "diffusion" upfront the position of the source term `Real matc "1.1*tx"`. Linked to the fact that the `Particle Time Integral` advects the `Source` term that depends itself on `Hadv`? 

[![Watch the video](https://github.com/cmosbeux/Damage/blob/main/AdvReac_Test/NoReinit_Nodal_Serial.avi)]

#### 1.1.2 Parallel

No problem linked to the parallelisation.


### 1.2 Elemental

The input variable is nodal, the ouput variable is elemental.

#### 1.2.1 Serial

In Serial, `Hadv` is well advected with no loss. `Source` remains to zero over the entire simulation. Notice that if we change `Hadv` for `Hpart` then we have an evolving source term. 

```f90
Particle Time Integral Source = Variable Hpart
  Real matc "1.1*tx"
```

Does `Particle Time Integral Source` needs a Nodal variable as an input?

#### 1.2.2 Parallel

No problem lined to the parallelisation.


## 2. Particle Reinitialization

In this case, we set `Reinitialize particles = Logical True`.
Here, the input variable is `HPart` and the solver outputs the advected value `Hadv`. We then udpate `Hpart= Equals Hadv` in the bodyforce using the `UpdateExported` solver

We also compute the `Particle Time Integral` that we define as `Source`.


### 2.1 Nodal 

### 2.1.1  Serial

Here we notice a clear numerical diffusion. Each timesetp uses the new interpolated `Hpart` as a restart and we loose information.

### 2.1.2 Parallel

No problem linked to the parallelisation.

###2.2 Elemental

### 2.2.1 Serial








