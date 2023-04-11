FC = elmerf90
FFLAGS = ''
# Object files
OBJ = ParticleAdvector.o MyParticleUtils.o

# Main target
all: ParticleAdvector

# Compile object files
ParticleAdvector.o: ParticleAdvector.F90
	$(FC) $(FFLAGS) -c ParticleAdvector.F90

MyParticleUtils.o: MyParticleUtils.F90
	$(FC) $(FFLAGS) -c MyParticleUtils.F90

# Link object files to create executable
ParticleAdvector: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o ParticleAdvector
  
# Clean up object files
clean:
	rm -f *.o
