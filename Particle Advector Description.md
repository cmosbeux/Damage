# Particle Advector 

## 1. Main routine 
The `ParticleAdvector` calculates the advection of fields using particle tracking. Given a field of some quantity (damage for example) and a transport field (velocity), this subroutine calculates how that quantity moves over time based on the velocityfield. 


At first, the code initializes the particles. The subroutines then loops over a specified number of timesteps and calculates the timestep size at each iteration.

The subroutine involves updating the particle positions based on the field being advected. This is done using an interpolation method to estimate the field values at each particle position, and then updating the particle position based on the velocity field.

# 2. Set the particle velocities

The subroutine is used to setup the velocities of each particle. 

```fortran90
SUBROUTINE SetParticleVelocities( FirstStep )
    LOGICAL :: FirstStep
    
    TYPE(Element_t), POINTER :: BulkElement
    ...
    ...
```

The subroutine is only called the first time. It then loops over the particle and check their status to skip lost, initiated, fixed, or on a wall boundary particles. If the particle is not in an element, then is is considered lost (`PARTICLE_LOSTS`).

 If the has moved and is still on an element, the subroutine retrieves the velocity at the particle's current position from the mesh.


```fortran90
 IF( .NOT. Visited ) THEN
      Mesh => GetMesh()
      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), dBasisdx(n, 3) )
      
      Params => GetSolverParams()
      
      VariableName = ListGetString(Params,'Velocity Variable Name',Found)
      IF(.NOT. Found) VariableName = 'flow solution'
      VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
        CALL Fatal('ParticleFieldInteraction','Velocity field variable does not exist: '//TRIM(VariableName))           
      END IF
      UseGradVelo = GetLogical( Params,'Velocity Gradient Correction',Found)

      IF(UseGradVelo .AND. Particles % RK2 ) THEN
        CALL Warn('ParticlePathIntegral','Quadratic source correction incompatibe with Runge-Kutta')
        UseGradVelo = .FALSE.
      END IF

      IF( .NOT. Particles % DtConstant ) THEN
        DtVar => ParticleVariableGet( Particles,'particle dt')
        IF( .NOT. ASSOCIATED( DtVar ) ) THEN
          CALL Fatal('SetParticleVelocities','Required field > particle dt < not present!')
        END IF
      END IF      

      SpeedMin = ListGetConstReal( Params,'Particle Min Speed',Found)
      IF(.NOT. Found) SpeedMin = EPSILON( SpeedMin )

      NewLost =  0

      Visited = .TRUE.
    END IF 
      
```

If the particle has moved and is still on an element, the subroutine retrieves the velocity at the particle's current position from the mesh. 

The velocity gradient correction term accounts for the variation of the velocity field across a particle's path and is used to improve the accuracy of the particle tracking simulation. However, in the context of the RK2 method, this term introduces additional numerical errors because the velocity gradient correction term is quadratic in the velocity field. The RK2 method already has inherent quadratic terms, and the addition of the quadratic velocity gradient correction term can cause numerical instability and inaccuracies. 

```fortran90
    Coordinate => Particles % Coordinate
    Velocity => Particles % Velocity
    Coord = 0.0_dp
    Velo = 0.0_dp

    OldLost = 0
    FixedLost = 0
  

    IF( Particles % DtConstant ) THEN
      dtime = Particles % DtSign * Particles % dtime
    END IF

    SkipZeroTime = .NOT. ( Particles % DtConstant .OR.  FirstStep ) 
    
  ```
  Now we will loop over the particles and check it is not lost.
  
  ```fortran90
  
    DO No = 1, Particles % NumberOfParticles
      Status = GetParticleStatus( Particles, No )
      IF( Status >= PARTICLE_LOST .OR. &
          Status <= PARTICLE_INITIATED .OR. &
          Status == PARTICLE_FIXEDCOORD .OR. &
          Status == PARTICLE_WALLBOUNDARY ) THEN
        OldLost = OldLost + 1
	CYCLE
      END IF

```

The element index containing the particle is retrieved and its velocity is only computed if the particle has moved.

```fortran90
      ElementIndex = GetParticleElement( Particles, No )
      IF( ElementIndex == 0 ) THEN
        Particles % Status(No) = PARTICLE_LOST
        NewLost(1) = NewLost(1) + 1
        CYCLE       
      END IF

      ! If the particle has not moved then it cannot have
      ! any change in the velocity.
      IF( SkipZeroTime ) THEN
        IF( ABS( DtVar % Values(No) ) < TINY( dtime ) ) CYCLE
      END IF
```

Element info and particle coordinates are retrieved.

```fortran90      
      BulkElement => Mesh % Elements( ElementIndex )
      Coord(1:dim) = Coordinate( No, 1:dim )
```

```fortran90
      !-------------------------------------------------------------------------
      ! Get velocity from mesh
      !-------------------------------------------------------------------------
      IF( UseGradVelo ) THEN
        stat = ParticleElementInfo( BulkElement, Coord, &
            SqrtElementMetric, Basis, dBasisdx )
      ELSE
        stat = ParticleElementInfo( BulkElement, Coord, &
            SqrtElementMetric, Basis )
      END IF

      IF(.NOT. Stat ) THEN
        Particles % Status(No) = PARTICLE_LOST
        NewLost(2) = NewLost(2) + 1
        CYCLE
      END IF

If `USeGradVelo` is true, the velocity is approximaterd by a first-order taylor series of the velocity field, accounting that velocity field may vary spatially. This approximation assumes that the veliocity field is smooth and continously differentiable.  

```fortran90
      IF( UseGradVelo ) THEN
        CALL GetVectorFieldInMesh(VeloVar,BulkElement, Basis, VeloAtPoint, &
            dBasisdx, GradVeloAtPoint )
	
  IF( .NOT. Particles % DtConstant ) THEN
          dtime = Particles % DtSign * DtVar % Values(No)
        END IF
        DO i=1,dim
          Velo(i) = VeloAtPoint(i) + &
              0.5_dp * SUM( GradVeloAtPoint(i,1:dim) * VeloAtPoint(1:dim) ) * dtime        
        END DO
      ELSE
        CALL GetVectorFieldInMesh(VeloVar, BulkElement, Basis, VeloAtPoint )
	Velo(1:dim) = VeloAtPoint(1:dim)
      END IF

      Speed = SQRT( SUM( Velo(1:dim) ** 2 ) )
      IF( Speed < SpeedMin ) THEN
 	Particles % Status(No) = PARTICLE_FIXEDCOORD
        Velocity( No, 1:dim ) = 0.0_dp
        FixedLost = FixedLost + 1
      ELSE
        Velocity( No, 1:dim ) = Velo(1:dim)
      END IF
    END DO
  ```
  
  If a particle is lost, then we call for some warnings:
  ```fortran90 
    IF( .FALSE. ) THEN
      IF( NewLost(1) > 0 ) CALL Warn('SetParticleVelocities','Some particles could not be located')
      PRINT *,'Total number of particles:',Particles % NumberOfParticles
      PRINT *,'Passive particles:',OldLost
      PRINT *,'New lost particles:',NewLost
      PRINT *,'New fixed velo particles:',FixedLost
    END IF
```

