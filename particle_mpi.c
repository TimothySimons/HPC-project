/*
 * Genetic algorithm for 2D Lennard Jones particle simulation
 * M. Kuttel October 2020
 */
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define DEFAULT_POP_SIZE 300 //bigger population is more costly
#define DEFAULT_NUM_PARTICLES 30 //more PARTICLES is more costly

// consts
static const int X_DEFAULT=20; //width of box
static const int Y_DEFAULT=20; //length of box
static const double MUTATION_RATE=0.10; //how often random mutations occur
static const double MAX_STAGNATION=100; // maximum number of generations
static const double MAX_GEN=1000; // maximum number of generations
static const double ITERATIONS=5; //number of times the whole process is run
static const double TOLERANCE=50; //not used... yet


//each person has x and y location in box
typedef struct{
    int x_pos;
    int y_pos;
} position;


// box pattern
typedef struct{
    position *person;
    double fitness;
} box_pattern;


typedef struct {
    int population_index;
    double fitness;
} population_best;


//display the box pattern
void printbox(box_pattern box,int num_particles){
    int i;
    for (i=0; i<num_particles-1; i++){
        printf("%d,%d\t", box.person[i].x_pos, box.person[i].y_pos);
    }
    printf("%d,%d\t:fitness %f\n",box.person[i].x_pos,box.person[i].y_pos,box.fitness);
}


//print the box pattern to file
void printboxFile(box_pattern box,FILE *f,int num_particles ){
    int i;
    for (i=0; i<num_particles-1; i++){
        fprintf(f,"%d,%d\t",box.person[i].x_pos,box.person[i].y_pos);
    }
    fprintf(f,"%d,%d\n",box.person[i].x_pos,box.person[i].y_pos,box.fitness);
}


/* FITNESS FUNCTION  - this is key*/
// NOTE: this is the most cost intensive operation of the program (by far) 
// NOTE: check out gprof - https://www.maketecheasier.com/profile-c-program-linux-using-gprof/
double calcFitness(box_pattern box, int num_particles){
    double fitness=0.0;
    int i,j;
    double x,y,r,tmp;
    for (i =0; i<num_particles-1; i++) {
        for (j=i+1; j<num_particles; j++) { //cycle through all pairs to calc distances
            x = (double)box.person[i].x_pos - (double)box.person[j].x_pos;
            y = (double) box.person[i].y_pos - (double)box.person[j].y_pos;
            r = sqrt((x*x)+(y*y));
            tmp = 2.0/r;
            //fitness-= 1.0/r; // electric repulsion
            //fitness-= pow(tmp,6); //purely repulsive function
            fitness+= (pow(tmp,12) - pow(tmp,6)); //Lennard-Jones function
        }
    }
    return fitness;
}


/* Creates initial random population */
void initPopulation(box_pattern *box, int population_size, int xmax, int ymax, int num_particles){
    int i,p;
    for (p=0; p<population_size; p++) {
        for (i=0; i<num_particles; i++){
            box[p].person[i].x_pos = (rand() % (xmax + 1));
            box[p].person[i].y_pos = (rand() % (ymax + 1));
        }
        box[p].fitness = calcFitness(box[p], num_particles);
    }
}


/* create child from parents */
box_pattern crossover(box_pattern child, box_pattern parentOne, box_pattern parentTwo, int splitPoint,int num_particles){
    int i=0;
    for (i=0; i<splitPoint; i++){ //copy over parentOne up to splitPoint
        child.person[i].x_pos=parentOne.person[i].x_pos;
        child.person[i].y_pos=parentOne.person[i].y_pos;
    }
    i--;
    if((rand()%(2) == 1) && (i<num_particles) && (i>=0)) //50% of time split in middle of person, more mixing
        child.person[i].y_pos=parentTwo.person[i].y_pos;

    for (i=splitPoint; i<num_particles; i++){ //copy over parentTwo from splitPoint to end
        child.person[i].x_pos=parentTwo.person[i].x_pos;
        child.person[i].y_pos=parentTwo.person[i].y_pos;
    }
    child.fitness = calcFitness(child, num_particles); //calculate fitness
    return child;
}


/* deep copy b into a [does a=b] */
void copybox(box_pattern *a, box_pattern *b,int num_particles){
    int i;
    for (i=0; i<num_particles; i++){
        (*a).person[i].x_pos=(*b).person[i].x_pos;
        (*a).person[i].y_pos=(*b).person[i].y_pos;
    }
    (*a).fitness=(*b).fitness;
}


/* main GA function - does selection, breeding, crossover and mutation */
population_best breeding(box_pattern *box, int population_size, int x_max, int y_max, int num_particles) {
    
    box_pattern max_parent; //keep track of highest from previous generation
    max_parent.person = malloc(num_particles * sizeof(position));
    copybox(&max_parent, &box[0], num_particles); //set max to first one
    
    int i;

    // allocate memory  
    box_pattern *new_generation = (box_pattern*) malloc(sizeof(box_pattern) * (population_size));
    for(i=0; i<population_size; i++) {
        new_generation[i].person = malloc(num_particles * sizeof(position));
    }

    for (i=0; i<population_size; i+=2) { //two children
        // determine breeding pair, with tournament of 2 (joust)
        int one = rand() % (population_size);
        int two = rand() % (population_size);
        int parentOne = two;
        if (box[one].fitness > box[two].fitness) parentOne = one;

        one = rand() % (population_size);
        two = rand() % (population_size);
        int parentTwo = two;
        if (box[one].fitness > box[two].fitness) parentTwo = one;
        
        int splitPoint = rand() % num_particles; 
        new_generation[i] = crossover(new_generation[i], box[parentOne], box[parentTwo], splitPoint, num_particles); //first child
        if (i+1 < population_size) {
            new_generation[i+1] = crossover(new_generation[i+1], box[parentTwo], box[parentOne], splitPoint, num_particles); //second child
        }

        // mutation first child
        double mutation = rand()/(double)RAND_MAX;
        if (mutation <= MUTATION_RATE ){
            int mutated = rand() % num_particles;
            new_generation[i].person[mutated].x_pos=(rand()%(x_max + 1));
            new_generation[i].person[mutated].y_pos=(rand()%(y_max + 1));
        }

        // mutation second child
        if (i+1 < population_size) {
            mutation = rand()/(double)RAND_MAX; //mutation second child
            if (mutation <= MUTATION_RATE ){
                int mutated = rand() % num_particles;
                new_generation[i+1].person[mutated].x_pos=(rand()%(x_max + 1));
                new_generation[i+1].person[mutated].y_pos=(rand()%(y_max + 1));
            }
        }
    }

    // find maximum parent fitness to keep and minimum new generation to throw away.
    // NOTE: if new_generation[0].fitness isnan then the fittest individual is never updated.
    // NOTE: this is because (x > NaN) always evaluates to false.
    new_generation[0].fitness = calcFitness(new_generation[0], num_particles);
    double min_fitness, max_fitness;
    if isnan(new_generation[0].fitness) {
        min_fitness = 0;
        max_fitness = 0;
    } else {
        min_fitness = new_generation[0].fitness;
        max_fitness = new_generation[0].fitness;
    }
    
    int min_box = 0;
    int highest = 0;
    for (i=1; i<population_size; i++){
        if (box[i].fitness > max_parent.fitness) {
            copybox(&max_parent, &box[i], num_particles); //replace lowest fitness with highest parent
        }
        new_generation[i].fitness = calcFitness(new_generation[i], num_particles);
        if (new_generation[i].fitness < min_fitness) {
            min_fitness=new_generation[i].fitness;
            min_box=i;
        }
        if (new_generation[i].fitness > max_fitness) {
            max_fitness=new_generation[i].fitness;
            highest=i;
        }
    }

    //copies
    for (i=0; i<population_size; i++){
        if (i == min_box) {
            copybox(&box[i], &max_parent, num_particles);
        } else {
            copybox(&box[i], &new_generation[i], num_particles);
        }
    }

    if (max_parent.fitness > max_fitness) {
        max_fitness = max_parent.fitness;
        highest = min_box;
    }

    population_best best_box;
    best_box.population_index = highest;
    best_box.fitness = max_fitness;

    //release memory
    for(i=0;i<population_size;i++)
        free(new_generation[i].person);
    free(new_generation);
    free(max_parent.person);

    return best_box;
}


void send_receive_best(int world_rank, int rank_count, population_best best, box_pattern *sub_population, int sub_pop_size, int num_particles) {
    if (world_rank == 0 || world_rank == rank_count) {
        int i;

        double best_fitness;
        int *xs = (int *) malloc(sizeof(int) * num_particles);
        int *ys = (int *) malloc(sizeof(int) * num_particles);

        if (world_rank == rank_count) {
            int best_index = best.population_index;
            best_fitness = best.fitness;
            for (i=0; i<num_particles; i++){
                xs[i] = sub_population[best_index].person[i].x_pos;
                ys[i] = sub_population[best_index].person[i].y_pos;
            }
            MPI_Ssend(&best_fitness, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Ssend(xs, num_particles, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Ssend(ys, num_particles, MPI_INT, 0, 2, MPI_COMM_WORLD);
        }
        
        if (world_rank == 0) {
            MPI_Recv(&best_fitness, 1, MPI_DOUBLE, rank_count, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(xs, num_particles, MPI_INT, rank_count, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(ys, num_particles, MPI_INT, rank_count, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int index = rand() % (sub_pop_size);
            if (best_fitness > sub_population[index].fitness) {
                sub_population[index].fitness = best_fitness;
                for (i=0; i<num_particles; i++){
                    sub_population[index].person[i].x_pos = xs[i]; 
                    sub_population[index].person[i].y_pos = ys[i]; 
                }
            }
        }

        free(xs);
        free(ys);
    }
}


bool is_finished(int world_rank, int stag, int max_stag){
    bool finished;
    if (world_rank == 0 && stag >= max_stag) {
        if (stag >= max_stag) {
            finished = true;
        } else {
            finished = false;
        }
    }
    MPI_Bcast(&finished, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return finished;
}


int main(int argc, char *argv[]) {
    int population_size = DEFAULT_POP_SIZE;
    int x_max = X_DEFAULT;
    int y_max = Y_DEFAULT;
    int num_particles = DEFAULT_NUM_PARTICLES;
    int iter = ITERATIONS;
    int k,i;
    if (argc >= 2) {
        population_size = atoi(argv[1]); //size population first command line argument
        if (argc >= 4) {
            x_max = atoi(argv[2]); //x dimension
            y_max = atoi(argv[3]); //y dimension
        }
        if (argc >= 5) num_particles = atoi(argv[4]);
        if (argc == 6) iter = atoi(argv[5]);
    }

    int gen_count = 0;
    
    // init
    int world_rank, world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    srand(world_rank);

    FILE *f;
    if (world_rank == 0) { //
        // only root process has access to the file
        f = fopen("solution.txt", "w"); 
        printf("Starting optimization with particles=%d, population=%d, width=%d,length=%d for %d iterations\n", num_particles, population_size, x_max, y_max, iter);
        printf("Writing dimensions to file\n");
        fprintf(f, "%d,%d\n", x_max, y_max); //write box dimensions as first line of file
    }

    // split population amongst processes
    int sub_pop_size = population_size / world_size;
    if (world_rank == 0) {
        sub_pop_size += population_size % world_size;
    }

    // allocate memory
    box_pattern *sub_population;
    sub_population = (box_pattern*) malloc(sizeof(box_pattern) * sub_pop_size);
    for(i=0;i<sub_pop_size;i++) {
        sub_population[i].person = malloc(num_particles * sizeof(position));
    }

    // run simulation k times
    for (k=0; k<iter; k++) {
        printf("Rank %d: initializing sub-population=%d iter=%d\n", world_rank, sub_pop_size, k);
        initPopulation(sub_population, sub_pop_size, x_max, y_max, num_particles);
        
        double start_time;
        if (world_rank == 0) {
            start_time = MPI_Wtime();
        }

        double max_fitness = 0;
        int gen = 0, highest = 0, stop = 0, rank_count = 1;
        while (gen < MAX_GEN) {
            population_best best_box = breeding(sub_population, sub_pop_size, 
                    x_max, y_max, num_particles);

            highest = best_box.population_index;
            double fitness = best_box.fitness;
            if (fitness <= max_fitness) {
                ++stop;
            } else {
                max_fitness = fitness;
                stop = 0;
            }

            send_receive_best(world_rank, rank_count, best_box, sub_population, 
                    sub_pop_size, num_particles);
            gen += 1;

            if (++rank_count % world_size == 0) {
                rank_count = 1;
            }

            if (is_finished(world_rank, stop, MAX_STAGNATION)) {
                break;
            }
        }

        double end_time;
        if (world_rank == 0) {
            end_time = MPI_Wtime();
            double total_time = end_time - start_time;
            printf("Root finished GA in %.4f seconds\n", total_time);
        }

        printf("Rank %d fitness: %f\n", world_rank, max_fitness);
        if (world_rank == 0) {
            if (f == NULL) {
                printf("Error opening file!\n");
                exit(1);
            }
            printboxFile(sub_population[highest], f, num_particles);
        }
        
        gen_count += gen;
    }

    if (world_rank == 0) {
        fclose(f);
    }

    // release memory
    for(i=0;i<sub_pop_size;i++) 
        free(sub_population[i].person); 
    free(sub_population);

    printf("Rank %d average generations: %f\n", world_rank, (double)gen_count/(double)k);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
