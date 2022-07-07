

# 1. Advection-Reaction 

From Krug et al. (2014), damage can describe in Elmer/Ice as a scalar value which is advected through the media according to the following advection-reaction equation (with a reaction term set to zero):

∂D/∂t+u ∇D=f(χ),

where u is the velocity vector from the Stokes system, and f(χ) the source term for damage.

# 2.	Semi-Lagrangian Advection

Another solution is to use a semi-Lagrangian particle advection scheme. By definition, a Lagrangian description of a system focuses on following individual particles along their trajectories (evolving over time) as opposed to the classical Eulerian description used in Elmer, which focuses on the variation of system variables at fixed locations (the grid points). The semi-Lagrangian scheme uses the Eulerian framework but with an Lagrangian discretization of the equations.

For the particle advection we assume that the fields are transported diffusion-free carried by a velocity field u. The particles are initialized at the nodes of the mesh and followed backward in time. When their position backward in time is recovered, it is communicated to the nodes as the value of the convected field. This is done by evaluating the following integral:

r ⃗=(r_0 ) ⃗+∫_0^(-δt)▒u dt.

When the particles have been advected, the field may be evaluated from 

D=D((r_0 ) ⃗,t)=D(r ⃗,t-δt),

where the damage D in each node depends on the initial position and the earlier timestep. 

A source term I(χ) that depends on the evolution of the particle along the path integral can be evaluated over space or time. In our case, that relates to the damage creation over δt.

I_t (χ)= ∫_0^(-δt)▒〖f(χ)〗 dt

We end up with

D=D((r_0 ) ⃗,t)+ I_t (χ)

An interpolation scheme is then utilized to estimate the value of the variable at the grid points surrounding the point where the particle originated from. 

# 3. Solver details

The `ParticleAdvector` solver relies on the `ParticleUtils` where a source term can be accounted for by checking the bodyforce `Particle Distance Integral Source` be `Particle Time Integral Source`. A variable that evolves over distance or time, is then created.

We can call the `ParticleAdvanceTimestep`, which allows particle to move with the velocity field, and a given timestep (defining the path and distance parkoured by the particle). 

The `ParticlePathIntegral` subroutine integrates variable over the path. This is activated by the source term in the `bodyforce`. The source term is evaluated at each node and added to the `TimeIntegVar`. The relevant code for an integral over time is:


