particle: particle.c
	 gcc -I -L particle.c -o particle

particle_omp: particle.c
	gcc -I -L particle.c -fopenmp -o particle_omp

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

clean: 
	rm *.exe *.out solution.txt solution.pdf

