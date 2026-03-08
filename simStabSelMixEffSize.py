import numpy as np
import argparse
from stabSelPopMixEffSize import *
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wright-Fisher forward simulation with polygenic adaptation and population structure.")
    parser.add_argument('--num_generations', type=int, default=50002, help='Number of generations to simulate.')
    parser.add_argument('--num_demes', type=int, default=100, help='Number of subpopulations.')
    parser.add_argument('--deme_size', type=int, default=200, help='Sizes of subpopulations.')
    parser.add_argument('--num_loci', type=int, default=200, help='Number of loci per subpopulation.')
    parser.add_argument('--Nu', type=float, default=1e-2, help='Mutation rates for subpopulations.')
    parser.add_argument('--recombination_rate', type=float, default= 0.5, help='Recombination rates for subpopulations.')
    parser.add_argument('--Vs', type=float, default= 1, help='Selection strengths for subpopulations.')
    parser.add_argument('--Nm', type=float, default=1, help='Migration rate between subpopulations.')
    parser.add_argument('--initial_allele_freq', type=float, default=12, help='The initial allele frequencies of allele 0. If initial_allele_freq is not a value in [0,1], draw allele frequencies from a beta(0.5,0.5).')
    parser.add_argument('--convergence_threshold', type=float, default=1e-4, help='The convergence_threshold indicates whether a population is at equilibrium.')
    parser.add_argument("--output_file", type=str, help="Name of the output text file")
    parser.add_argument("--prop", type=float, default=1, help="Proportion of NuL")
    parser.add_argument("--effSize_file", type=str, default="/none")


    args = parser.parse_args()

    if os.path.exists(args.effSize_file):
        r1=np.loadtxt(args.effSize_file)
        effectSize=np.delete(r1[0],0)
    else:
        Ns=np.random.exponential(scale=args.Nu*args.num_loci*args.prop, size=args.num_loci)  
        effectSize=np.round(np.sqrt((Ns*2*args.Vs/args.deme_size)),4)

    print(effectSize)

    simulation = Simulation(
        num_generations=args.num_generations,
        num_demes=args.num_demes,
        deme_size=args.deme_size,
        num_loci=args.num_loci,
        Nu=args.Nu,
        recombination_rate=args.recombination_rate,
        selection_strength=1/(args.Vs*2),
        Nm=args.Nm,
        effectSize=effectSize,
        initial_allele_freq=args.initial_allele_freq,
        convergence_threshold=args.convergence_threshold,
        output_file=args.output_file
    )

    simulation.run()

### As variable effect sizes