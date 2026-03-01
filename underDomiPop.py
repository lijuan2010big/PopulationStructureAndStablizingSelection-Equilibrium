import numpy as np
import argparse

def build_population(pop_size, p):
    # Generate individuals with initial genotypes for one locus with specified allele frequencies
    p00 = p**2
    p11 = (1 - p)**2
    p01 = 1-p00-p11
    population = np.random.choice(3, pop_size,p=[p00, p01, p11])
    return population.astype(int)

def wright_fisher_simulation(num_demes, pop_size_per_deme, max_generations, Nm, Ns, Nu, init_allele_freq, convergence_threshold, outfile):
    total_pop_size = num_demes * pop_size_per_deme
    mutation_rate = Nu/pop_size_per_deme
    s=Ns/pop_size_per_deme
    migration_rate=Nm/pop_size_per_deme

    # Initialize populations for each deme
    populations = [build_population(pop_size_per_deme, init_allele_freq) for i in range(num_demes)]
    allele_freqs_history = []
    out_fmt = ['%d'] + ['%1.6f'] * num_demes
    out_file = open(outfile, 'w')

    for gen in range(max_generations):
        ## Record allele frequency for all generations
        if gen%500==0:
            print(gen)
            # Calculate allele frequencies for each deme
            output_info=[gen]
            allele_freqs = [1-np.mean(population)/2 for population in populations]
            output_info.extend(allele_freqs)
            np.savetxt(out_file, np.column_stack(output_info), fmt=out_fmt)
 
        # Migration
        # Determine the number of migrants based on the migration rate
        num_migrants = np.random.binomial(total_pop_size, migration_rate)
        if num_migrants>0:
            # Create a migration pool by concatenating populations from all demes
            migration_pool = np.random.choice(np.concatenate(populations), num_migrants)
            # Randomly select individuals for substitution by individuals in migrant pool
            dest_deme=np.random.randint(0, num_demes, num_migrants)
            dest_ind=np.random.randint(0, pop_size_per_deme, num_migrants)
            for i in range(num_migrants):
                # Replace individual in the deme with the migrant
                populations[dest_deme[i]][dest_ind[i]] = migration_pool[i]

        # Mutation
        # Introduce mutations in the population
        num_mutations = np.random.binomial(2*total_pop_size, mutation_rate)
        dest_deme=np.random.randint(0, num_demes, num_mutations)
        dest_ind=np.random.randint(0, pop_size_per_deme, num_mutations)

        for i in range(num_mutations):
        # Mutate the allele: change it to a random allele (0, 1, or 2) different from the current one
            current_allele = populations[dest_deme[i]][dest_ind[i]]
            new_allele = 1 if current_allele != 1 else np.random.choice([0,2])
            populations[dest_deme[i]][dest_ind[i]] = new_allele

        # Selection and Reproduction (random mating within each deme)
        for i in range(num_demes):
            population = populations[i]
            genotype_freq = np.array([np.sum(population == 0) , np.sum(population == 1), np.sum(population == 2)])/pop_size_per_deme
            new_genotype_freq = genotype_freq*[1, 1-s, 1]
            new_genotype_freq = new_genotype_freq / sum(new_genotype_freq)
            p = new_genotype_freq[0]+new_genotype_freq[1]/2
            populations[i] = build_population(pop_size_per_deme, p)


    out_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wright-Fisher simulation with underdominance selection, migration between demes, and mutations at one locus")
    parser.add_argument("--num_demes", type=int, default=100, help="Number of demes")
    parser.add_argument("--pop_size", type=int, default=100, help="Population size per deme")
    parser.add_argument("--max_generations", type=int, default=100002, help="Maximum number of generations")
    parser.add_argument("--Nm", type=float, default=0.1, help="Migration rate (proportion of individuals involved in migration)")
    parser.add_argument("--Ns", type=float, default=0.1, help="Selection coefficient (against heterozygotes)")
    parser.add_argument("--Nu", type=float, default=1e-2, help="Mutation rate (per allele per generation)")
    parser.add_argument("--init_allele_freq", type=float, default=0, help="Frequency of allele 0 in the initial population")
    parser.add_argument("--convergence_threshold", type=float, default=1e-8, help="Threshold for convergence in allele frequencies")
    parser.add_argument("--outfile", type=str, help="Name of the output text file")
    
    args = parser.parse_args()

    # Run simulation
    wright_fisher_simulation(args.num_demes, args.pop_size, args.max_generations, args.Nm, args.Ns, args.Nu, args.init_allele_freq, args.convergence_threshold, args.outfile)

