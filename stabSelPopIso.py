import numpy as np
import argparse

class Population:
    def __init__(self, pop_size, num_loci, mutation_rate, recombination_rate, selection_strength,allele_freq):
        self.pop_size = pop_size # Population size
        self.num_loci = num_loci
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate
        self.selection_strength = selection_strength
        self.genotypes = np.zeros((pop_size, num_loci, 2), dtype=int) # Diploid: 2 alleles per locus

        for locus in range(num_loci):
            self.genotypes[:, locus, 0] = np.random.binomial(1, allele_freq[locus], pop_size)
            self.genotypes[:, locus, 1] = np.random.binomial(1, allele_freq[locus], pop_size)

    ### Mutation
    def mutate(self):
        mutation_events = np.random.rand(self.pop_size, self.num_loci, 2) < self.mutation_rate
        self.genotypes ^= mutation_events

    #### Recombination
    def recombine(self, parent):
        select_alleles = np.random.rand(self.pop_size, self.num_loci) < self.recombination_rate
        gamete = np.where(select_alleles[:,:], parent[:, :,0], parent[:, :,1])
        return gamete

    #### Trait features
    def traits(self):
        # assume allele 0 has effect -0.05, allele 1 has effect 0.05. The effect change of a substitution from A to a is 0.1. The optimal is in the middle of two extremes [-10, 10].
        return np.sum(self.genotypes, axis=(1,2))*0.1

    def traits_mean(self):
        return np.mean(self.traits())

    def traits_variance(self):
        return np.var(self.traits())

    #### Selection
    def reproduce(self, optimum):
        fitness = np.exp(-self.selection_strength *(self.traits()-self.traits_mean())*(self.traits()+self.traits_mean()-2*optimum))
        probabilities = fitness / np.sum(fitness)


        parent1_genotypes=self.genotypes[np.random.choice(self.pop_size, size=self.pop_size, p=probabilities)]
        parent2_genotypes=self.genotypes[np.random.choice(self.pop_size, size=self.pop_size, p=probabilities)]

        self.genotypes[:,:,0]=self.recombine(parent1_genotypes)
        self.genotypes[:,:,1]=self.recombine(parent2_genotypes)

    def heterozygosity(self):
        frequencies = self.allele_frequencies()
        heterozygosity = 2 * frequencies * (1-frequencies)
        return heterozygosity

    def allele_frequencies(self):
        allele_counts = np.sum(self.genotypes == 0, axis=(0, 2))
        frequencies = allele_counts / (2 * self.pop_size)
        return frequencies

class Simulation:
    def __init__(self, num_generations, pop_size, num_loci,
                 Nu, recombination_rate, selection_strength, initial_allele_freq, convergence_threshold, output_file):
        self.num_generations = num_generations
        self.pop_size=pop_size
        self.num_loci=num_loci
        self.convergence_threshold = convergence_threshold
        self.traits_variance_history = []
        self.output_file = output_file

        mutation_rate=Nu/pop_size

        if initial_allele_freq >= 0 and initial_allele_freq <= 1 :
            allele_freq = np.random.choice([1-initial_allele_freq, initial_allele_freq], num_loci, p=[0.5,0.5])
            print(allele_freq)
        else:
            print("Beta")
            allele_freq = np.random.beta(4*Nu, 4*Nu, num_loci)

        self.population = Population(pop_size, num_loci, mutation_rate, recombination_rate, selection_strength, allele_freq)


    def run(self):
        # Save to file
        out_file = open(self.output_file, 'w')
        out_fmt = ['%d'] + ['%1.6f'] * 3
        header= ['Gen'] + ["het"] + ["traitMean"] +  ["traitVar"] 
        np.savetxt(out_file, np.column_stack(header), fmt='%s')

        # Simulation
        optimum0=self.num_loci*0.1*2*0.5
        print(optimum0)
        for gen in range(self.num_generations):
            if gen%100 == 0:
                print(gen)
                # collect output info
                output_info=[gen]
                output_info.append(np.mean(self.population.heterozygosity()))
                output_info.append(self.population.traits_mean())
                output_info.append(self.population.traits_variance())
                # save to a file
                np.savetxt(out_file, np.column_stack(output_info), fmt=out_fmt)

            self.population.mutate()
            self.population.reproduce(optimum0)

        out_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wright-Fisher forward simulation with polygenic adaptation and population structure.")
    parser.add_argument('--num_generations', type=int, default=40002, help='Number of generations to simulate.')
    parser.add_argument('--num_demes', type=int, default=100, help='Number of subpopulations.')
    parser.add_argument('--deme_size', type=int, default=100, help='Sizes of subpopulations.')
    parser.add_argument('--num_loci', type=int, default=100, help='Number of loci per subpopulation.')
    parser.add_argument('--Nu', type=float, default=1e-2, help='Mutation rates for subpopulations.')
    parser.add_argument('--recombination_rate', type=float, default= 0.5, help='Recombination rates for subpopulations.')
    parser.add_argument('--Ns', type=float, default= 0.05, help='Selection strengths for subpopulations.')
    parser.add_argument('--Nm', type=float, default=1, help='Migration rate between subpopulations.')
    parser.add_argument('--initial_allele_freq', type=float, default=12, help='The initial allele frequencies of allele 0. If initial_allele_freq is not a value in [0,1], draw allele frequencies from a beta(0.5,0.5).')
    parser.add_argument('--convergence_threshold', type=float, default=1e-4, help='The convergence_threshold indicates whether a population is at equilibrium.')
    parser.add_argument("--output_file", type=str, help="Name of the output text file")

    args = parser.parse_args()

    FST=1/(1+4*args.Nm+8*args.Nu)
    pop_size=round(args.num_demes*args.deme_size/(1-FST))

    Vs=pop_size*0.1*0.1/(args.Ns*2)
    # Create and run the simulation
    simulation = Simulation(
        num_generations=args.num_generations,
        pop_size=pop_size,
        num_loci=args.num_loci,
        Nu=args.Nu,
        recombination_rate=args.recombination_rate,
        selection_strength=1/(Vs*2),
        initial_allele_freq=args.initial_allele_freq,
        convergence_threshold=args.convergence_threshold,
        output_file=args.output_file
    )

    simulation.run()
