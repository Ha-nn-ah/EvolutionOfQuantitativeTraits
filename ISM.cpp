// Simulation-Code for haploid WF population
    // Infinite Sites Model: Adaptation from New Mutations of an Additive quantitative Trait with Unlinked Loci

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
        int numloci; // number of loci where a mutation happened (but maybe lost again)
        int N; // population size: often used 10000
        double s; // selection rate
        double mu_locus; // mutation rate
        int f; // distribution mutation effects (0 = equal, 1 = exponential, 2 = truncated normal)
        double a; // mean mutation effect
        vector<double> selArray; // the mutation effects for each locus
        double Stop; // stopping condition (number of generations the population is evolving)
        int Dtau; // after how many generations an output is generated
        int t;
        int seed; // seed for random number generator- Pseudo-random, as we use a fixed seed for the same run
        vector<int> popArray; // the number of mutant individual at numloci in the population
        double popArrayWeights[2]; // the weights for multinomial sampling. holds only for one locus at a time.
        unsigned int popArrayLocus[2]; // an array to use in evolve, this will hold information for a single locus, and is just used to copy and paste from popArray
    
        int numlocifix; // number of fixed mutated loci
        int numlociloss; // number of extinct mutated loci
        double ZMean; // phenotypic mean in the population
        double ZVar; // phenotypic variance in the population
        vector<int> numlociArray; // number of loci where a mutation happened every Dtau generations
        vector<int> numlocifixArray; // number of fixed mutated loci every Dtau generations
        vector<int> numlocilossArray; // number of extinct mutated loci every Dtau generations
        vector<double> ZMeanArray; // phenotypic mean in the population every Dtau generation
        vector<double> ZVarArray; // phenotypic variance in the population every Dtau generation
    
    void Initialization (double pop_size,double sel,double mutation_rate,int mut_effect,double mean_effect,double tau,int delta_tau,int selected_seed); // Initialise the random number generator from the GNU-Library. We choose the MT19937 generator. Intitialise arrays for population, age and fitness.
        bool Evolve (void);
        void Mutation (void);
        void Reproduction (gsl_rng *rng); 
};

void Generation::Initialization(double pop_size,double sel,double mutation_rate,int mut_effect,double mean_effect,double tau,int delta_tau,int selected_seed){
    numloci = 0; // ancestral population in the beginning
    N = pop_size;
    s = sel;
    mu_locus = mutation_rate;
    f = mut_effect;
    a = mean_effect;
    Stop = tau;
    Dtau = delta_tau;
    t = 0;
    seed = selected_seed;
    for(int i=0;i<2;i++){
      popArrayLocus[i]=0;
      popArrayWeights[i]=0;
    }
    numlocifix=0;
    numlociloss=0;
    ZMean=0;
    ZVar=0;
}

void Generation::Reproduction(gsl_rng *rng){
                       
    for (int l = 0; l<numloci; l++){
        
        if (popArray[l] > 0) {
            
            //***selection***
                double fitness = exp(s*selArray[l]); // exponential directional selection model
                popArrayWeights[0] = (N-popArray[l])/(double)N; // frequency of the ancestral state at this locus (fitness of the ancestral state is normalized to 1)
                popArrayWeights[1] = popArray[l]/(double)N*fitness; // frequency of the derived state at this locus times the fitness of the derived state
     
            //***reproduction**
                gsl_ran_multinomial(rng, 2, N, popArrayWeights, popArrayLocus);
     
            popArray[l]=popArrayLocus[1]; // how many individuals are now from mutant type at locus l
        }
    }
    
    //***mutation***
        popArrayWeights[0] = 1-mu_locus; // probability of no mutation
        popArrayWeights[1] = mu_locus; // probability of a mutation at a new locus

        gsl_ran_multinomial(rng, 2, N, popArrayWeights, popArrayLocus);
    
    // add loci which mutated in this generation
    for (int i=0; i<popArrayLocus[1];i++){
        popArray.push_back(1);
        
        if (f == 0) { // equal mutation effects
            selArray.push_back(a);
        } else if (f == 1){ // exponential distributed mutation effects
            selArray.push_back(gsl_ran_exponential(rng, a));
        } else if (f == 2){ // truncated normal distributed mutation effects
            selArray.push_back(abs(gsl_ran_gaussian(rng, a*1.2533141373155001))); // normal distributed N[0,v] truncated \pi a^2 = 2 v^2
        } else {return;}
               
        numloci++;
    }
    
    // save output every Dtau generations
    if ((t+1) % Dtau == 0) {
        numlocifix = 0; // fixed mutants
        numlociloss = 0; // extinct mutants
        ZMean = 0; // mean
        ZVar = 0; // variance
        for (int i=0; i<numloci; i++){
            ZMean = ZMean + selArray[i]*(popArray[i]/((double)N));
            ZVar = ZVar + pow(selArray[i],2)*(popArray[i]/((double)N))*(1-(popArray[i]/((double)N)));
            if (popArray[i] == N) {
                numlocifix++;
            }
            if (popArray[i] == 0) {
                numlociloss++;
            }
        }
        ZMeanArray.push_back(ZMean);
        ZVarArray.push_back(ZVar);
        numlociArray.push_back(numloci);
        numlocifixArray.push_back(numlocifix);
        numlocilossArray.push_back(numlociloss);
    }
}

bool Generation::Evolve(void){
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937); // Generate a random number generator called rng, of type gsl_rng
    gsl_rng_set (rng, (unsigned long)seed);
    
    // loop: each loop is one generation, continues until the population reaches stopping condition (certain generation has been reached)
    
    Reproduction(rng); // generation 0
	
    while (t<Stop){ // generation 1 to StopCond
        t++;
        seed++;
        Reproduction(rng);
    }
    
    gsl_rng_free(rng);
    return 1;
}

int main(){

    // the following needs to be adjusted:
    int NumR = 20000; // number of simulation runs (typically 20000)
    int PopSize = 10000; // population size (typically 10000)
    double SelCoeff = 0.001; // selection coefficient
    double MutRate = 0.5; // population scaled mutation rate Theta = N*mu
    int MutEffect = 0; // distribution mutation effects (0 = equal, 1 = exponential, 2 = truncated normal)
    double MeanEffect = 1; // mean mutation effect
    double StopCond = 10; // stopping condition (number of generations)
    int DeltaT = 1; // after how many generations an output is generated
    
    ofstream myfileMean;
    myfileMean.open ("mean.tsv");
    ofstream myfileVar;
    myfileVar.open ("variance.tsv");
    ofstream myfileMut;
    myfileMut.open ("mutants.tsv");
    ofstream myfileFix;
    myfileFix.open ("fixed.tsv");
    ofstream myfileLoss;
    myfileLoss.open ("extinct.tsv");
    
    for (int j=0;j<NumR;j++){
        Generation gen;
        gen.Initialization (PopSize,SelCoeff,MutRate/((double)PopSize),MutEffect,MeanEffect,StopCond,DeltaT,j);

        if(!gen.Evolve()){
            NumR++;
        }
        else{
            
            for(int i=0;i<gen.numlociArray.size();i++){ // number of mutated loci for each run
              myfileMut << gen.numlociArray[i] <<"\t";
            }
            myfileMut << endl;
            
            for(int i=0;i<gen.numlocifixArray.size();i++){ // number of fixed mutated loci for each run
              myfileFix << gen.numlocifixArray[i] <<"\t";
            }
            myfileFix << endl;
            
            for(int i=0;i<gen.numlocilossArray.size();i++){ // number of extinct mutated loci for each run
              myfileLoss << gen.numlocilossArray[i] <<"\t";
            }
            myfileLoss << endl;
            
            for(int i=0;i<gen.ZMeanArray.size();i++){ // phenotypic mean of the population for each run
              myfileMean << gen.ZMeanArray[i] <<"\t";
            }
            myfileMean << endl;
            
            for(int i=0;i<gen.ZVarArray.size();i++){ // phenotypic variance of the population for each run
              myfileVar << gen.ZVarArray[i] <<"\t";
            }
            myfileVar << endl;
        }
    }
    
    myfileMean.close();
    myfileVar.close();
    myfileMut.close();
    myfileFix.close();
    myfileLoss.close();
    return 0;
}
