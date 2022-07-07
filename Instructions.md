

# 1. Advection-Reaction 

From Krug et al. (2014), damage can describe in Elmer/Ice as a scalar value which is advected through the media according to the following advection-reaction equation (with a reaction term set to zero):

$\dfrac{\partial D}{\partial t} + \boldsymbol{u} \nabla D= f(\chi)$,

where $\boldsymbol{u}$ is the velocity vector from the Stokes system, and $f(\chi)$ the source term for damage.

# 2.	Semi-Lagrangian Advection

Another solution is to use a semi-Lagrangian particle advection scheme. By definition, a Lagrangian description of a system focuses on following individual particles along their trajectories (evolving over time) as opposed to the classical Eulerian description used in Elmer, which focuses on the variation of system variables at fixed locations (the grid points). The semi-Lagrangian scheme uses the Eulerian framework but with an Lagrangian discretization of the equations.

For the particle advection we assume that the fields are transported diffusion-free carried by a velocity field $\boldsymbol{u}$. The particles are initialized at the nodes of the mesh and followed backward in time. When their position backward in time is recovered, it is communicated to the nodes as the value of the convected field. This is done by evaluating the following integral:

$\vec{r} = \vec{r_0}+\int_{0}^{-\delta t} \boldsymbol{u} \ dt$.

When the particles have been advected, the field may be evaluated from 

$D = D (\vec{r_0},t) = D(\vec{r},t-\delta t)$,

where the damage D in each node depends on the initial position and the earlier timestep. 

A source term I(χ) that depends on the evolution of the particle along the path integral can be evaluated over space or time. In our case, that relates to the damage creation over δt.

$I_t (\chi)= \int_{0}^{-\delta t} f(\chi) \ dt$

We end up with

$D=D(\vec{r_0},t)+ I_t (\chi)$

An interpolation scheme is then utilized to estimate the value of the variable at the grid points surrounding the point where the particle originated from. 

# 3. Solver details

The `ParticleAdvector` solver relies on the `ParticleUtils` where a source term can be accounted for by checking the bodyforce `Particle Distance Integral Source` be `Particle Time Integral Source`. A variable that evolves over distance or time, is then created.

We can call the `ParticleAdvanceTimestep`, which allows particle to move with the velocity field, and a given timestep (defining the path and distance parkoured by the particle). 

The `ParticlePathIntegral` subroutine integrates variable over the path. This is activated by the source term in the `bodyforce`. The source term is evaluated at each node and added to the `TimeIntegVar`. The relevant code for an integral over time is (`line 5597-5611`):

```f90
      ! Path integral over time
      IF( TimeInteg ) THEN
        Source(1:n) = ListGetReal( BodyForce,'Particle Time Integral Source', &
            n, Indexes, Found )
        IF( Found ) THEN
          SourceAtPath = SUM( Basis(1:n) * Source(1:n) )
          IF( UseGradSource ) THEN
            DO i=1,dim
              SourceAtPath = SourceAtPath + 0.5*SUM( dBasisdx(1:n,i) * Source(1:n) ) * &
                  ( PrevCoord(i) - Coord(i) )
            END DO
          END IF
          TimeIntegVar % Values(No) = TimeIntegVar % Values(No) + dtime * SourceAtPath
        END IF
      END IF
```

The code also allows for integral over distance but we do not use it for damage advection. 

## Application to damage

We want to account for the damage creation by using a source term. It can be described in the bodyforce section using the `damage USF` as a source for the Particle Time Integrale variable:

```f90
Body Force 1
  Particle Time Integral Source = Variable "Damage"
    Real Procedure "USF_Damage" "SourceDamage"
End
```

The integral as to be defined as a time integral as the development of the damage is calculated at each timestep, for a given χ (Chi, which can be exported at the dime of the ComputeDevStress solver). 

The solver can be executed at the end of each timestep, when Navier-Stokes and the stresses have been computed. `Particle Time integral` is one of the internal variable of the `ParticleAdvector` and can be used for the damage advection. The keyword `Result Variable 3 = String "Damage"` can be used to export the integrated damage (instead of Particle Time Integral). 

Here we describe line by line the options we use in the solver:

* The particules can be initialized at the center of each element (`Logical True`) or at the node (`Logical False`)
```f90
Advect Elemental = Logical True
```

* We can decide to either reinitialize the particlue at each timestep or not. 
```f90
 Reinitialize Particles = Logical False
 Particle Dt Constant = Logical False
```

* We can describe the "internal" Timestepping strategy (the maximum number of timestep to backtrace the particule along its path)
```f90
Simulation Timestep Sizes = Logical True
Max Timestep Intervals = Integer 5
```

* Time in average 4 steps in each element
```f90
Timestep Unisotropic Courant Number = Real 0.25
Max Timestep Size = Real 1.0e5
```

* Give up integration if particles are tool old:
```f90
Max Integration Time = Real 1.0e5
```

* The velocity variable name can be given so that the particule can follow the flow. In our case, we compute the flow with Stokes and the solution is the `Flow Solution`:
```f90
Velocity Variable Name = String "Flow Solution"
```

* Runge Kutta can be used for a forward integration in time. In our case we set the parameter to False, as we only backtrace the particules. Gradient and source gradient corrections can be added. This is an alternative way of increasing the accuracy of the integral with respect to the forward Runge Kutta scheme. Here the gradient of the velocity field and the gradient of the source are evaluated at the point of the particle to account for the curvature of the flow. 
```f90
Runge Kutta = Logical False
Velocity Gradient Correction = Logical True
Source Gradient Correction = Logical True

```

* We can decide to display some info at the end of the solver (optional)
```f90
Show some info in the end
Particle Info = Logical True
Particle Time = Logical True
```

* Some internal variable of the solver can be calculated and given as an output. In our case, the damage is associated to the `particle time integral` with its source term provided in the bodyforce. The `particle time integral` is then copied to the result variable `Damage` so that it can be used to update the source term at the next timestep.
```f90
! The internal variables for this solver
  Variable 1 = String "Particle distance"
  Variable 2 = String "Particle time"
  Variable 3 = String "Particle time integral"
  Variable 4 = String "Particle distance integral"
  Variable 5 = String "Particle Velocity_Abs"
  
  Result Variable 3 = String "Damage"
```

The complete Solver section can be seen here after:

```f90
Solver 10
  Equation = ParticleAdvector
  Procedure = "ParticleAdvector" "ParticleAdvector"

! Initialize particles at center of elements (as opposed to nodes)
  Advect Elemental = Logical True

  Reinitialize Particles = Logical False
  Particle Dt Constant = Logical False

! Timestepping strategy
  Simulation Timestep Sizes = Logical True
  Max Timestep Intervals = Integer 5

! Time in average 4 steps in each element
  Timestep Unisotropic Courant Number = Real 0.25
  Max Timestep Size = Real 1.0e5

! Give up integration if particles are tool old
  Max Integration Time = Real 1.0e5

  Velocity Variable Name = String "Flow Solution"

! Integration forward in time (we do not use Runge Kutta)
  Runge Kutta = Logical False
! But we add gradient corrections
  Velocity Gradient Correction = Logical True
  Source Gradient Correction = Logical True

! Show some info in the end
  Particle Info = Logical True
  Particle Time = Logical True

! The internal variables for this solver
  Variable 1 = String "Particle distance"
  Variable 2 = String "Particle time"
  Variable 3 = String "Particle time integral"
  Variable 4 = String "Particle distance integral"
  Variable 5 = String "Particle Velocity_Abs"
  Result Variable 3 = String "Damage"
End

```

For boundary conditions, a particle wall can be applied where the Particle should not leave the body (e.g., at the bedrock interface). 
