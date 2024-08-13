#!/bin/bash

#BSUB -P acc_kennylab
#BSUB -n 1
#BSUB -W 12:00
#BSUB -q premium
#BSUB -R rusage[mem=8000]
#BSUB -J "fst[1]"
#BSUB -o /sc/arion/projects/igh/kennylab/christa/fst/logs/out.%J.%I 
#BSUB -e /sc/arion/projects/igh/kennylab/christa/fst/logs/err.%J.%I 
 

. /sc/arion/projects/igh/kennylab/christa/miniconda/etc/profile.d/conda.sh
conda activate ibd

num=$LSB_JOBINDEX 
genotype_file="../../data/biome/geno/gght_v2_topmed_allchr"
fam_file="louvain_redone_clusters.txt" 
output="louvain_redone/fst_maf"

# run plink fst calculation over final louvain cluster assignments 
# assumes a file where column 1 is IID, column 2 is FID, and column 3 is the final cluster assignment for an individual 
# ex: 7241534	7241534	cluster1_1_1_1
# uses the PCA sites, not sure if this matters too much 

../plink2 --bfile $genotype_file \
    --maf 0.05 \
    --extract "../mapping_ibd/pca/pruned_ids.prune.in" \
    --pheno $fam_file \
    --fst PHENO1 'method'='hudson' \
    --out $output