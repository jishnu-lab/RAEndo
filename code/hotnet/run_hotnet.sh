#!/bin/bash

#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=hotnet_RA_pgr_permute

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --ntasks-per-node=4

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

#SBATCH --error=/ix/djishnu/Priyamvada/Immport/Network_analysis/hotnet.err

# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
#hotnet2=..
#num_cores=-1
#num_network_permutations=50
#num_heat_permutations=1000

module purge 
module load python/anaconda2.7-4.2.0


# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
#hotnet2=..
num_cores=-1

python /ix/djishnu/TFH/hotnet2-master/HotNet2.py \
	-nf  /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5 \
	-pnp /ix/djishnu/Priyamvada/Immport/Network_analysis/Networks/Networks_to_run_hotnet/permuted_networks/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4_##NUM##.h5 \
	-hf  /ix/djishnu/Priyamvada/Immport/forJavad/LDAK_score_RA_prg65.csv \
	-np 500 \
	-hp 500 \
	-o /ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute \
	-c  $num_cores
