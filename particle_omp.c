/*
 * Genetic algorithm for 2D Lennard Jones particle simulation
 * M. Kuttel October 2020
 */

#include <math.h>
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
static const double ITERATIONS=10; //number of times the whole process is run
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
    #pragma omp parallel private(p,i)
    {
        for (p=0; p<population_size; p++) {
            for (i=0; i<num_particles; i++){
                // if (omp_get_thread_num()==0) {
                //     printf("threadNum=%d,i=%d,p=%d\n",omp_get_thread_num(),i,p);
                // }
                box[p].person[i].x_pos = (rand() % (xmax + 1));
                box[p].person[i].y_pos = (rand() % (ymax + 1));
            }
            box[p].fitness = calcFitness(box[p], num_particles);
        }
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


/* Main GA function - does selection, breeding, crossover and mutation */
population_best breeding(box_pattern *box, int population_size, int x_max, int y_max, int num_particles){
    box_pattern max_parent; //keep track of highest from previous generation
    max_parent.person = malloc(num_particles * sizeof(position));
    copybox(&max_parent, &box[0], num_particles); //set max to first one
    int i;
    box_pattern *new_generation = (box_pattern*) malloc(sizeof(box_pattern) * (population_size));
    // no parallel here since malloc isn't an expensive function?

    for(i=0;i<population_size;i++) {
        new_generation[i].person=malloc(num_particles*sizeof(position));
    }

    int one,two,parentOne,parentTwo,splitPoint,mutated;
    double mutation;

    // TODO: add #pragma omp parrallel for directive to this for loop.
    // NOTE: reason being it contains calcFitness (bottleneck)
    #pragma omp parallel private(i,one, two, parentOne,parentTwo,splitPoint,mutation,mutated) shared(box,new_generation,num_particles)
//  this parallel has too much overhead and slows down the program with small params
    {
        for (i=0; i<population_size; i+=2) { //two children
            // printf("(#threads=%d) iteration: %i\n",omp_get_thread_num(), i);
            // Determine breeding pair, with tournament of 2 (joust)
            one = rand() % (population_size);
            two = rand() % (population_size);
            parentOne = two;
            if (box[one].fitness > box[two].fitness) parentOne = one; //joust

            one = rand() % (population_size);
            two = rand() % (population_size);
            parentTwo = two;
            if (box[one].fitness > box[two].fitness) parentTwo=one; //joust

            splitPoint = rand() % num_particles; //split chromosome at point
            new_generation[i] = crossover(new_generation[i], box[parentOne], box[parentTwo], splitPoint, num_particles); //first child
            new_generation[i+1] = crossover(new_generation[i+1], box[parentTwo], box[parentOne], splitPoint, num_particles); //second child

            // Mutation first child
            mutation = rand()/(double)RAND_MAX;
            if (mutation <= MUTATION_RATE ){
                mutated = rand() % num_particles;
                new_generation[i].person[mutated].x_pos=(rand()%(x_max + 1));
                new_generation[i].person[mutated].y_pos=(rand()%(y_max + 1));
            }
            mutation = rand()/(double)RAND_MAX; //mutation second child
            if (mutation <= MUTATION_RATE ){
                mutated = rand() % num_particles;
                new_generation[i+1].person[mutated].x_pos=(rand()%(x_max + 1));
                new_generation[i+1].person[mutated].y_pos=(rand()%(y_max + 1));
            }
        }
    }
    // find maximum parent fitness to keep and minimum new generation to throw away.
    // NOTE: if new_generation[0].fitness isnan then the fittest individual is never updated.
    // NOTE: this is because (x > -nan) always evaluates to false.
    // #pragma omp barrier // reason for barrier is that the access to new_generation below might cause a data race
    new_generation[0].fitness = calcFitness(new_generation[0], num_particles);
    double min_fitness, max_fitness;
    if (isnan(new_generation[0].fitness)) {
        min_fitness = 0;
        max_fitness = 0;
    } else {
        min_fitness = new_generation[0].fitness;
        max_fitness = new_generation[0].fitness;
    }
    int min_box = 0;
    int highest = 0;

    // TODO: add #pragma omp parrallel for directive to this for loop.
    // TODO: ensure comparison against min and max are in a #pragma omp critical region.
    // NOTE: use named critical regions (name for min; diff. name for max)
    // NOTE: reason for parallelise - contains calcFitness (bottleneck)
    #pragma omp parallel private(i) shared(box, max_parent,num_particles, new_generation,min_fitness,min_box,max_fitness,highest)
    {   
        // printf("%d\n", omp_get_thread_num());
        for (i=1; i<population_size; i++){

            if (box[i].fitness > max_parent.fitness) {
                copybox(&max_parent, &box[i], num_particles); //replace lowest fitness with highest parent
            }
            new_generation[i].fitness = calcFitness(new_generation[i], num_particles);
            #pragma omp critical
            if (new_generation[i].fitness < min_fitness) {
                min_fitness=new_generation[i].fitness;
                min_box=i;
            }
            #pragma omp critical
            if (new_generation[i].fitness > max_fitness) {
                max_fitness=new_generation[i].fitness;
                highest=i;
            }
        }
    }

    //copies
    #pragma omp parallel private(i) shared(box, max_parent,num_particles, new_generation,min_box)
    {   
        for (i=0; i<population_size; i++){
            //printbox(new_generation[i]);
            if (i == min_box) {
                copybox(&box[i], &max_parent, num_particles);
            } else {
                copybox(&box[i], &new_generation[i], num_particles);
            }
            // printbox(box[i]);
        }
    }
    if (max_parent.fitness > max_fitness) { //previous generation has the best
        max_fitness = max_parent.fitness;
        highest = min_box;
    }

    population_best best_box;
    best_box.population_index = highest;
    best_box.fitness = max_fitness;


    // #pragma omp parallel private(i) shared(new_generation)
    // {
        for(i=0;i<population_size;i++) {
            free(new_generation[i].person); //release memory
        }
    // }
    free(new_generation); //release memory
    free(max_parent.person);
    return best_box;
}


int main(int argc, char *argv[] ) {
    int population_size=DEFAULT_POP_SIZE;
    int x_max =X_DEFAULT;
    int y_max=Y_DEFAULT;
    int num_particles=DEFAULT_NUM_PARTICLES;
    int iter=ITERATIONS;
    int k,i;
    if (argc >=2) {
        population_size = atoi(argv[1]); //size population first command line argument
        if (argc>=4) {
            x_max=atoi(argv[2]); //x dimension
            y_max=atoi(argv[3]); //x dimension
        }
        if (argc>=5) num_particles=atoi(argv[4]);
        if (argc==6) iter =atoi(argv[5]);
    }

    printf("Starting optimization with particles = %d, population=%d, width=%d,length=%d for %d iterations\n",num_particles,population_size,x_max,y_max,iter);
    int gen_count = 0;
    FILE *f = fopen("solution.txt", "w");
    printf("Writing dimensions to file\n");
    fprintf(f, "%d,%d\n", x_max, y_max); //write box dimensions as first line of file

    box_pattern * population;
    population = (box_pattern*) malloc(sizeof(box_pattern) * population_size); //allocate memory
    for(i=0;i<population_size;i++)
        population[i].person = malloc(num_particles * sizeof(position)); //allocate memory

    double start; 
    double end; 
    double averageTimeTaken = 0;
    start = omp_get_wtime(); 

    for (k=0; k<iter; k++) { //k is number of times whole simulation is run
        // populate with initial population
        printf("initializing population\n");
        initPopulation(population,population_size,x_max,y_max,num_particles);
        printf("=========%d\n", k);

        double max_fitness = 0;
        int stop = 0;
        int gen = 0, highest = 0;
        while (stop < MAX_STAGNATION && gen < MAX_GEN) {
            population_best best_box = breeding(population, population_size, x_max, y_max, num_particles);
            highest = best_box.population_index;
            double fitness = best_box.fitness;
            if (fitness <= max_fitness) {
                ++stop;
            } else {
                max_fitness = fitness;
                stop = 0;
            }
            gen += 1;
        }
        end = omp_get_wtime(); 
        double timeTaken = end - start;
        printf("Work took %f seconds\n", timeTaken);
        averageTimeTaken+=timeTaken;

        printf("# generations= %d \n", gen);
        printf("Best solution:\n");
        printbox(population[highest], num_particles);
        if (f == NULL) {
            printf("Error opening file!\n");
            exit(1);
        }
        printboxFile(population[highest], f, num_particles);
        printf("---------");
        gen_count += gen;
        start = omp_get_wtime(); 
    }
    fclose(f);

    for(i=0;i<population_size;i++) 
        free(population[i].person); //release memory
    free(population); //release memory

    printf("Average generations: %f\n", (double)gen_count/(double)k);

    averageTimeTaken=averageTimeTaken/ (double)iter;
    printf("Average time taken: %f with %i particles" , averageTimeTaken, num_particles);
    return 0;
}
