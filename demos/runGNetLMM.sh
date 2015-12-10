#!/bin/bash 

BFILE=./../data/1000G_chr22/chrom22_subsample20_maf0.10 #specify here bed basename
FFILE=./../data/1000G_chr22/ones.txt
PFILE=./out/pheno
CFILE=./out/chrom22
ASSOC0FILE=./out/lmm
GFILE=./out/genes

ANCHOR_THRESH=1e-6
ANCHORFILE=./out/cisanchor_thresh1e-6_wnd2000.txt
WINDOW=2000

VFILE=./out/vstructures_thresh1e-6_wnd2000
ASSOCFILE=./out/gnetlmm_thresh1e-6_wnd2000

PLOTFILE=./out/power.pdf

# Generate phenotypes
#./../GNetLMM/bin/gNetLMM_simPheno --bfile $BFILE --pfile $PFILE

# Compute covariance matrix
#./../GNetLMM/bin/gNetLMM_preprocess --compute_covariance --bfile $BFILE --cfile $CFILE

# Run initial association scan
#for i in $(seq 0 10000 40000)
#do
#    ./../GNetLMM/bin/gNetLMM_analyse --initial_scan --bfile $BFILE --pfile $PFILE --cfile $CFILE.cov --assoc0file $ASSOC0FILE.startSnp_$i --startSnpIdx $i --nSnps 10000 --ffile $FFILE
#done

# Merging results
#./../GNetLMM/bin/gNetLMM_analyse --merge_assoc0_scan  --assoc0file $ASSOC0FILE --nSnps 10000 --bfile $BFILE
# Compute marginal gene-gene correlations
#./../GNetLMM/bin/gNetLMM_analyse --gene_corr --pfile $PFILE --gfile $GFILE 

# Compute anchors 
#./../GNetLMM/bin/gNetLMM_analyse --compute_anchors  --bfile $BFILE --pfile $PFILE --assoc0file $ASSOC0FILE --anchorfile $ANCHORFILE --anchor_thresh=$ANCHOR_THRESH  --window=$WINDOW --cis


for i in $(seq 0 10 90)
do
    ./../GNetLMM/bin/gNetLMM_analyse --find_vstructures  --pfile $PFILE  --gfile $GFILE --anchorfile $ANCHORFILE  --assoc0file $ASSOC0FILE --window $WINDOW --vfile $VFILE.startTrait_$i --bfile $BFILE --startTraitIdx $i --nTraits 10
    
     ./../GNetLMM/bin/gNetLMM_analyse --update_assoc --bfile $BFILE --pfile $PFILE --cfile $CFILE.cov --ffile $FFILE --vfile $VFILE.startTrait_$i --assocfile $ASSOCFILE.startTrait_$i --startTraitIdx $i --nTraits 10
done

# Merge csv files
#./../GNetLMM/bin/gNetLMM_analyse --concatenate --files $VFILE
#./../GNetLMM/bin/gNetLMM_analyse --concatenate --files $ASSOCFILE

# Write to matrix
#./../GNetLMM/bin/gNetLMM_analyse --merge_assoc --assoc0file $ASSOC0FILE --assocfile $ASSOCFILE

# Plot results
#./../GNetLMM/bin/gNetLMM_postprocess --assocfile $ASSOCFILE --assoc0file $ASSOC0FILE --plotfile $PLOTFILE --pfile $PFILE --bfile $BFILE --window $WINDOW



