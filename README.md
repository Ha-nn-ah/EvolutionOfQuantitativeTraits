--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
	README.txt for the Wright-Fisher simulations
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------

The simulation package contains four C++ files and one R file:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%	OneLocus.cpp	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This program tracks the number of mutants for one locus. The simulation code has been used for Sect. 3 and 4 and App. E.

- Input: parameters
	population size N, selection coefficient s, 
	population-scaled mutation rate \tilde{Theta} (0 for main part), 
	initial number of mutants in the population (1 for main part, 0 for Appendix C),
	stopping condition n (number of generations), 
	time step dt (output generated every dt generations)

- Output: simulation.tsv
	i-th column = number of mutants after i*dt generations; for i = 1,2,...,n/dt
	rows = runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%	ISM_Loci.cpp	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This program evaluates the number of mutants at each locus according to the Infinite Sites Model. It has been used for Figure 3.5.

- Input: parameters
	population size N, selection coefficient s, 
	population-scaled mutation rate \Theta, 
	mutation effect distribution f, mean mutation effect,	
	stopping condition \tau (number of generations)

- Output: simulation.tsv
	i-th column = number of mutants at locus i in generation \tau
	rows = runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ISM_SegSites.cpp	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This program evaluates the number of mutants at each locus according to the Infinite Sites Model. It has been used for Figure 6.1.

- Input: parameters
	population size N, selection coefficient s, 
	population-scaled mutation rate \Theta, 
	mutation effect distribution f, mean mutation effect,	
	stopping condition \tau (number of generations) or G (mean phenotype reached in the population)

- Output: simulation.tsv
	1st column = number of simulated generations (\tau)
	2nd column = number of mutated loci, where ancestral state get fixed again during \tau generations
	3rd column = mean phenotype in the population after \tau generations
	i-th column = number of mutants arising from the (i-3)-th successful mutation event in generation \tau; i = 4,5,...
	rows = runs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%	ISM.cpp		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This program keeps track of different quantities (see output) in the Infinite Sites Model. After summarizing with the help of ISM_summary.R, the results have been used for Sect. 4 and 5 and App. D and E.

- Input: parameters
	population size N, selection coefficient s, 
	population-scaled mutation rate \Theta, 
	mutation effect distribution f, mean mutation effect,	
	stopping condition \tau (number of generations),
	time step dt (output generated every dt generations)

- Output: 
	i-th column = xxx after i*dt generations; for i = 1,2,...,tau/dt
	rows = runs
	
	mean.tsv	xxx = phenotypic mean
	variance.tsv	xxx = phenotypic variance
	mutants.tsv	xxx = number of mutation events
	fixed.tsv	xxx = number of mutated loci, where derived state get fixed
	extinct.tsv	xxx = number of mutated loci, where ancestral state get fixed again

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%	ISM_summary.R	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script merges the five output files of ISM.cpp to one tsv-file.

- Input: output files of ISM.cpp
	parameter dt used in ISM.cpp

- Output: simulation.tsv
	1st column: generation
	2nd column: expected phenotypic mean
	3rd column: expected phenotypic variance
	4th column: expected number of mutation events
	5th column: expected number of mutated loci, where derived state get fixed
	6th column: expected number of mutated loci, where ancestral state get fixed again
	7th column: expected number of segregating sites
	8th column: variance of the phenotypic mean
	9th column: variance of the phenotypic variance
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
--------------------------------------------------------------------------------------
	OXS-TERMINAL
--------------------------------------------------------------------------------------

The C++ scripts are executed by...

g++ $(gsl-config --cflags) OneLocus.cpp $(gsl-config --libs) && ./a.out

g++ $(gsl-config --cflags) ISM_Loci.cpp $(gsl-config --libs) && ./a.out

g++ $(gsl-config --cflags) ISM.cpp $(gsl-config --libs) && ./a.out
