all: mesh compile mpi 
mesh:
	ElmerGrid 1 2 square.grd
	ElmerGrid 1 2 square.grd -metis 2 0 0
compile:
	elmerf90 USFs.F90 -o USFs

serial:
	ElmerSolver 

mpi:
	mpirun -np 2 ElmerSolver_mpi
