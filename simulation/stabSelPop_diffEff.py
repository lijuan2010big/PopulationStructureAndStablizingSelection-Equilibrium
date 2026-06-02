import numpy as np
import os.path
import argparse

class Population:
    def __init__(self, size, num_loci, mutation_rate, recombination_rate, selection_strength, effect_sizes, allele_freq):
        self.size = size # Population size
        self.num_loci = num_loci
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate
        self.selection_strength = selection_strength
        self.effect_sizes = np.array(effect_sizes)
        self.genotypes = np.zeros((size, num_loci, 2), dtype=int) # Diploid: 2 alleles per locus
        self.optimum = np.sum(effect_sizes)

        for locus in range(num_loci):
            self.genotypes[:, locus, 0] = np.random.binomial(1, allele_freq[locus], size)
            self.genotypes[:, locus, 1] = np.random.binomial(1, allele_freq[locus], size)

    ### Mutation
    def mutate(self):
        mutation_events = np.random.rand(self.size, self.num_loci, 2) < self.mutation_rate
        self.genotypes ^= mutation_events

    #### Recombination
    def recombine(self, parent):
        select_alleles = np.random.rand(self.size, self.num_loci) < self.recombination_rate
        gamete = np.where(select_alleles[:,:], parent[:, :,0], parent[:, :,1])
        return gamete


    #### Trait features
    def traits(self):
        # Variable effect. 
        return (np.sum(np.sum(self.genotypes, axis=2)*self.effect_sizes, axis=1))

    def traits_mean(self):
        return np.mean(self.traits())

    def traits_variance(self):
        return np.var(self.traits())

    #### Selection
    def reproduce(self):
        # The trait range is 0, 2 * num_loci * locus_effect. We assume the middle as the optimum num_loci * locus_effect
        fitness = np.exp(-self.selection_strength *(self.traits()-self.traits_mean())*(self.traits()+self.traits_mean()-2*self.optimum))
        probabilities = fitness / np.sum(fitness)

        parent1_genotypes=self.genotypes[np.random.choice(self.size, size=self.size, p=probabilities)]
        parent2_genotypes=self.genotypes[np.random.choice(self.size, size=self.size, p=probabilities)]

        self.genotypes[:,:,0]=self.recombine(parent1_genotypes)
        self.genotypes[:,:,1]=self.recombine(parent2_genotypes)


    def heterozygosity(self):
        frequencies = self.allele_frequencies()
        heterozygosity = 2 * frequencies * (1-frequencies)
        return heterozygosity

    def allele_frequencies(self):
        allele_counts = np.sum(self.genotypes == 0, axis=(0, 2))
        frequencies = allele_counts / (2 * self.size)
        return frequencies

    def genetic_variance(self):
        return np.sum(self.heterozygosity() * self.effect_sizes**2)

class Simulation:
    def __init__(self, num_generations, num_demes, subpopulation_size, num_loci,
                 Nu, recombination_rate, selection_strength, Nm, initial_allele_freq, convergence_threshold, effect_sizes, initFrq_file, output_file):
        self.num_generations = num_generations
        self.num_demes = num_demes
        self.subpopulation_size=subpopulation_size
        self.total_pop_size = num_demes * subpopulation_size
        self.migration_rate = Nm/subpopulation_size
        self.convergence_threshold = convergence_threshold
        self.traits_variance_history = []
        self.output_file = output_file
        self.effect_sizes = effect_sizes

        self.c1= np.where(effect_sizes == 0.05, 1, 0)
        self.c1_num=np.sum(self.c1)
        self.c2= np.where(effect_sizes == 0.1, 1, 0)
        self.c2_num=np.sum(self.c2)
        self.c3= np.where(effect_sizes == 0.2, 1, 0)
        self.c3_num=np.sum(self.c3)

        mutation_rate=Nu/subpopulation_size

        #if os.path.isfile(initFrq_file):
        #    allele_freq=np.loadtxt(initFrq_file)
        #    allele_freq=np.transpose(allele_freq)
        #    for i in range(num_loci):
        #        np.random.shuffle(allele_freq[i])
        #    np.random.shuffle(allele_freq)
        #    allele_freq=np.transpose(allele_freq)
        #    print(allele_freq)
        #    self.populations = [Population(subpopulation_size, num_loci, mutation_rate, recombination_rate, selection_strength, effect_sizes, allele_freq[i]) for i in range(num_demes)]
        #else:
        #if initial_allele_freq >= 0 and initial_allele_freq <= 1 :
        #    allele_freq = np.full(num_loci, 1-initial_allele_freq)
        #else:
            # suppose to draw from beta distribution with alpha=beta=4*N*mutation_rate with expected heterozygosity 0.04.
        #    allele_freq = np.random.beta(4*Nu, 4*Nu, num_loci)
        #self.populations = [Population(subpopulation_size, num_loci, mutation_rate, recombination_rate, selection_strength, effect_sizes, allele_freq) for _ in range(num_demes)]

        self.populations = [Population(subpopulation_size, num_loci, mutation_rate, recombination_rate, selection_strength, effect_sizes, np.random.beta(4*Nu, 4*Nu, num_loci)) for _ in range(num_demes)]


    def check_convergence(self):
        if len(self.traits_variance_history) < 20:
            return False  # Insufficient history for convergence check
        recent_freqs1 = np.array(self.traits_variance_history[-10:-5])
        recent_freqs2 = np.array(self.traits_variance_history[-5:])
        max_diff = np.abs(np.mean(recent_freqs1) - np.mean(recent_freqs2))
        return max_diff < self.convergence_threshold

    def metapop_heterozygosity(self):
        frequencies = []
        for i, pop in enumerate(self.populations):
            frequencies.append(pop.allele_frequencies())
        frequencies=np.mean(frequencies, axis=0)
        heterozygosity = 2 * frequencies * (1-frequencies)
        return heterozygosity

    def metapop_trait(self):
        traits_in_metapop = []
        for i, pop in enumerate(self.populations):
            traits_in_metapop.append(pop.traits())
        return np.concatenate(traits_in_metapop)

    def metapop_trait_mean(self):
        return np.mean(self.metapop_trait())

    def metapop_trait_variance(self):
        return np.var(self.metapop_trait())

    def run(self):
        # Save to file
        out_file = open(self.output_file, 'w')
        out_fmt = ['%d'] + ['%1.6f'] * 5 +['%1.6f'] * ((self.num_demes)*6)
        header= ['Gen'] 
        for i in range(self.num_demes+1):
            if i>0:
                header = header + ["c1_het"+str(i)] + ["c2_het"+str(i)] + ["c3_het"+str(i)] + ["traitMean"+str(i)] +  ["traitVar"+str(i)] + ["Vg"+str(i)]
            else:
                header = header + ["c1_het"+str(i)] + ["c2_het"+str(i)] + ["c3_het"+str(i)] + ["traitMean"+str(i)] +  ["traitVar"+str(i)]
        np.savetxt(out_file, np.column_stack(header), fmt='%s')

        # Simulation
        for gen in range(self.num_generations):
            # output trait mean, trait variance, heterozygosity in each population
            if gen%200 == 0:
                print(gen)
                # collect output info
                output_info=[gen]
                output_info.append(np.dot(self.metapop_heterozygosity(), self.c1)/self.c1_num)
                output_info.append(np.dot(self.metapop_heterozygosity(), self.c2)/self.c2_num)
                output_info.append(np.dot(self.metapop_heterozygosity(), self.c3)/self.c3_num)
                output_info.append(self.metapop_trait_mean())
                output_info.append(self.metapop_trait_variance())
                for pop in self.populations:
                    output_info.append(np.dot(pop.heterozygosity(), self.c1)/self.c1_num)
                    output_info.append(np.dot(pop.heterozygosity(), self.c2)/self.c2_num)
                    output_info.append(np.dot(pop.heterozygosity(), self.c3)/self.c3_num)
                    output_info.append(pop.traits_mean())
                    output_info.append(pop.traits_variance())
                    output_info.append(pop.genetic_variance())
                # save to a file
                np.savetxt(out_file, np.column_stack(output_info), fmt=out_fmt)

            for pop in self.populations:
                pop.mutate()
                pop.reproduce()

            # Migration
            # Determine the number of migrants based on the migration rate
            #num_migrants=int(self.mean_num_migrants+np.random.rand())
            num_migrants=np.random.binomial(self.total_pop_size, self.migration_rate)

            if num_migrants>0:
                # Randomly select individuals to build migrant pool
                inds = np.random.randint(0, self.total_pop_size, num_migrants)
                source_deme,source_ind = divmod(inds, self.subpopulation_size)
                migration_pool =  [self.populations[source_deme[i]].genotypes[source_ind[i]] for i in range(num_migrants)]
                # Randomly select individuals for substitution by individuals in migrant pool
                inds = np.random.randint(0, self.total_pop_size, num_migrants)
                dest_deme,dest_ind = divmod(inds, self.subpopulation_size)
                # Copy individuals from the source deme to the destine deme.
                for i in range(num_migrants):
                    while source_deme[i] == dest_deme[i]:
                        dest_deme[i] = np.random.randint(0,self.num_demes)
                    self.populations[dest_deme[i]].genotypes[dest_ind[i]]=migration_pool[i]
                    
           # self.traits_variance_history.append(self.metapop_trait_variance())


        out_file.close()

import numpy as np
import argparse
from stabSelPop_diffEff import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wright-Fisher forward simulation with polygenic adaptation and population structure.")
    parser.add_argument('--num_generations', type=int, default=40002, help='Number of generations to simulate.')
    parser.add_argument('--num_demes', type=int, default=100, help='Number of subpopulations.')
    parser.add_argument('--subpopulation_size', type=int, default=100, help='Sizes of subpopulations.')
    parser.add_argument('--num_loci', type=int, default=100, help='Number of loci per subpopulation.')
    parser.add_argument('--Nu', type=float, default=1e-2, help='Mutation rates for subpopulations.')
    parser.add_argument('--recombination_rate', type=float, default= 0.5, help='Recombination rates for subpopulations.')
    parser.add_argument('--Vs', type=float, default= 0.5, help='Selection strengths for subpopulations. Vs: the variance of the fitness function.')
    parser.add_argument('--Nm', type=float, default=1, help='Migration rate between subpopulations.')
    parser.add_argument('--initial_allele_freq', type=float, default=12, help='The initial allele frequencies of allele 0. If initial_allele_freq is not a value in [0,1], draw allele frequencies from a beta(0.5,0.5).')
    parser.add_argument('--convergence_threshold', type=float, default=1e-6, help='The convergence_threshold indicates whether a population is at equilibrium.')
    parser.add_argument("--initFrq_file", type=str, help="Name of the initial allele frquency file")
    parser.add_argument("--output_file", type=str, help="Name of the output text file")

    args = parser.parse_args()

    effect_sizes = np.array([0.05]*int(args.num_loci*0.6) + [0.1]*int(args.num_loci*0.3) + [0.2]*int(args.num_loci*0.1))
    np.random.shuffle(effect_sizes)

    #Vs default is 0.5, so that Ns is 0.25, 1, 4

    # Create and run the simulation
    simulation = Simulation(
        num_generations=args.num_generations,
        num_demes=args.num_demes,
        subpopulation_size=args.subpopulation_size,
        num_loci=args.num_loci,
        Nu=args.Nu,
        recombination_rate=args.recombination_rate,
        selection_strength=1/(args.Vs*2),
        Nm=args.Nm,
        initial_allele_freq=args.initial_allele_freq,
        convergence_threshold=args.convergence_threshold,
        effect_sizes=effect_sizes,
        initFrq_file=args.initFrq_file,
        output_file=args.output_file
    )

    simulation.run()


