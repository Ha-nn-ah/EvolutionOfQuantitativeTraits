// Simulation-Code for haploid WF population
    // One Locus Model: Adaptation from New Mutations

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <sstream> 
#include <iostream> 
#include <fstream>
#include <vector> 
#include <math.h>
#include <string>
#include <cstdlib>
// #include <omp.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <time.h>
#include <algorithm>  

using namespace std;

class Generation{
    public:
        int N; // population size: often used 10000
        double s; // selection rate
        double mu_locus; // mutation rate
        double BackMut; // 1 = backmutation same rate as mutation; 0 = no backmutation
        double Stop; // stopping condition (number of generations the population is evolving)
        int Dtau; // after how many generations an output is generated
        int t;
        int seed; // seed for random number generator- Pseudo-random, as we use a fixed seed for the same run
        vector<int> popArray; // the number of mutant individual
        int popMut; // the number of mutant individual
        double popWeights[2]; // the weights for multinomial sampling
        double popWeightsBeforeMut[2];
        unsigned int popLocus[2]; // an array to use in evolve
        
  void Initialization (double pop_size,double sel,double mutation_rate,int initial_mut,double tau,int delta_tau,int selected_seed); // Initialise the random number generator from the GNU-Library. We choose the MT19937 generator.
        bool Evolve (void);
        void Mutation (void);
        void Reproduction (gsl_rng *rng); 
};

void Generation::Initialization(double pop_size,double sel,double mutation_rate,int initial_mut, double tau,int delta_tau,int selected_seed){
    N = pop_size;
    s = sel;
    mu_locus = mutation_rate;
    BackMut = 0; // 1 = backmutation same rate as muation; 0 = no backmutation
    Stop = tau;
    Dtau = delta_tau;
    t = 0;
    seed = selected_seed;
    popMut = initial_mut;
    for(int i=0;i<2;i++){
      popLocus[i]=0;
      popWeights[i]=0;
    }
} 

void Generation::Reproduction(gsl_rng *rng){
    
    //***selection***
      double fitness = exp(s); // exponential directional selection model
      popWeightsBeforeMut[0] = (N-popMut)/(double)N; // frequency of the ancestral state (fitness of the ancestral state is normalized to 1)
      popWeightsBeforeMut[1] = popMut/(double)N*fitness; // frequency of the derived state times the fitness of the derived state
                     
    //***mutation***
      popWeights[0] = popWeightsBeforeMut[0]*(1-mu_locus) + popWeightsBeforeMut[1]*mu_locus*BackMut;
      popWeights[1] = popWeightsBeforeMut[1]*(1-mu_locus*BackMut) + popWeightsBeforeMut[0]*mu_locus;

    //***reproduction**
      gsl_ran_multinomial(rng, 2, N, popWeights, popLocus);
     
    popMut=popLocus[1]; // how many individuals are now from mutant type
    
    // save the number of mutants in the population every Dtau generations
    if ((t+1) % Dtau == 0) {
        popArray.push_back(popMut);
    }
}

bool Generation::Evolve(void){
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937); // Generate a random number generator called rng, of type gsl_rng
    gsl_rng_set (rng, (unsigned long)seed);
    
    // loop: each loop is one generation, continues until the population reaches stop condition (certain generation has been reached)
	while (t<Stop){
      Reproduction(rng);
      t++;
      seed++;
    }
    gsl_rng_free(rng);
    return 1;
}

int main(){

    // the following needs to be adjusted:
    int NumR = 20000; // number of simulation runs (typically 20000)
    int PopSize = 10000; // population size (typically 10000)
    double SelCoeff = 0.1; // selection coefficient
    double MutRate = 0; // population scaled mutation rate Theta = N*mu
    int InitialMut = 0; // set to 0 (all individuals have the ancestral state in the beginning); or 1 (start with one mutant)
    double StopCond = 100; // stopping condition (number of generations)
    int DeltaT = 25; // after how many generations an output is generated
    
    ofstream myfile;
    myfile.open ("simulation.tsv");
    
    for (int j=0;j<NumR;j++){
        Generation gen;
        gen.Initialization (PopSize,SelCoeff,MutRate/((double)PopSize),InitialMut,StopCond,DeltaT,j);

        if(!gen.Evolve()){
            NumR++;
        }
        else{
            for(int i=0;i<gen.popArray.size();i++){ // number of derived alleles in the total population for each run
                myfile << gen.popArray[i] <<"\t";
            }
            myfile << endl;
        }
    }
    myfile.close();
    return 0;
}
