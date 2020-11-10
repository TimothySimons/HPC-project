particle: particle.c
	 gcc -I -L particle.c -o particle

particle_omp: particle.c
	gcc -I -L particle.c -fopenmp -o particle_omp

particle_mpi: particle_mpi.c
	mpicc particle_mpi.c -o particle_mpi

debug: 
	gcc -g -I -L particle.c -o debug 

debug_omp: 
	gcc -g -I -L particle.c -fopenmp -o debug_omp

run: particle
	./particle 2000 100 100 10 5
	cat solution.txt
	python plot_solution.py solution.txt

run_omp: particle_omp
	./particle_omp 2000 100 100 10 5
	cat solution.txt
	python plot_solution.py solution.txt

run_mpi: particle_mpi
	mpiexec -n 4 particle_mpi 100 20 20 10 5

clean: 
	rm *.exe *.out solution.txt solution.pdf

