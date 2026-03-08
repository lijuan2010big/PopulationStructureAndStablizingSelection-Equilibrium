import numpy as np
import argparse
from underDomiPop import *

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

