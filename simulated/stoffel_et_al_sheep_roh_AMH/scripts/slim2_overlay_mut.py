# Keywords: Python, tree-sequence recording, tree sequence recording

from sys import argv
import inspect
import msprime, pyslim, gzip
import tskit
import numpy as np
#import matplotlib.pyplot as plt
import os
import random

# get args
script, run_name, demo_mod, pop_size = argv
infile = "../output/" + demo_mod + "/trees/" + run_name + ".trees"
outfile = "../output/" + demo_mod + "/vcfs/" + run_name + ".vcf"

pop_size = int(pop_size)
# tree sequence
# ts = pyslim.load(infile) # "slim_sim/sheep.trees"
ts = tskit.load(infile) # "slim_sim/sheep.trees"

# alive = ts.individuals_alive_at(0)
# print(f"There are {len(alive)} individuals alive from the final generation.")
#ts.individual(201).id

# why was each individuals retained in tree sequence?
# indiv_types = {"first_gen" : 0,
#                "remembered" : 0,
#                "alive" : 0}
# for ind in mutated.individuals():
#    if ind.flags & pyslim.INDIVIDUAL_FIRST_GEN:
#       indiv_types['first_gen'] += 1
#    if ind.flags & pyslim.INDIVIDUAL_REMEMBERED:
#       indiv_types['remembered'] += 1
#    if ind.flags & pyslim.INDIVIDUAL_ALIVE:
#       indiv_types['alive'] += 1
# 
# for k in indiv_types:
#    print(f"Number of individuals that are {k}: {indiv_types[k]}")
   
   
# is recapitation necessary?
# recapitation adds diversity present in the initial generation;
# will it make a difference? In fact, most segments of the genome have
# already coalesced, i.e. 13379 out of 15471 have one root
sum([t.num_roots == 1 for t in ts.trees()])
sum([t.num_roots > 0 for t in ts.trees()])

# recapitate anyway.
recap = pyslim.recapitate(ts, recombination_rate=9.15e-9, ancestral_Ne=pop_size, random_seed=np.random.randint(1,30000000))

# verify that it worked
# orig_max_roots = max(t.num_roots for t in ts.trees())
# recap_max_roots = max(t.num_roots for t in recap.trees())
# print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
#       f"and after recapitation, it was {recap_max_roots}.")


# some individuals might not be simplified away because they have nodes that
# are required to describe the genealogies of the sample
# keep_indivs = np.random.choice(recap.individuals_alive_at(0), 200, replace=False)
# keep_nodes = []
# for i in keep_indivs:
#    keep_nodes.extend(recap.individual(i).nodes)
# recap_simp = recap.simplify(keep_nodes)

# add mutations. Also wrap in SlimTreeSequence so that pyslim can still work with it.
mutated = msprime.mutate(recap, rate=4.17e-9, 
                                  random_seed=np.random.randint(1,30000000), 
                                  keep=True)
# print(f"The tree sequence now has {mutated.num_mutations} mutations, "
#       f"and mean pairwise nucleotide diversity is {mutated.diversity()}.")

# only keep individuals that are alive today
# mutated=mutated.simplify()

n_dip_indv = int(mutated.num_samples / 2)
ind_names=[0] * n_dip_indv
for i in range(n_dip_indv):
      ind_names[i] = mutated.individual(i).metadata["pedigree_id"]

#ind_names_sub = random.sample(ind_names, 200)

#indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
indv_names = [f"tsk_{str(i)}indv" for i in ind_names]

alive = pyslim.individuals_alive_at(mutated, 0).tolist()

# subset 100 individuals
ind_sub = random.sample(range(len(alive)), 100)
alive_sub = [alive[i] for i in ind_sub]
indv_names_sub = [indv_names[i] for i in ind_sub]

with open(outfile, "w") as vcf_file:
    mutated.write_vcf(vcf_file, individuals = alive_sub, individual_names=indv_names_sub)
    








  
