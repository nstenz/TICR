#!/bin/bash

# Input data settings
MB_BLOCK="mb-block.txt"
DATASET_ROOT="chr4-subset"
DATASET="$DATASET_ROOT.nex"

# Fiddle around with these and see how things change!
MDL_BLOCK_SIZE=10
CONVERGENCE_THRESHOLD=0.02

START=$(date +%s)

# In case something goes wrong we want to stop immediately,
# in which case please contact us with your error
set -e

# Reset number of generations and sample frequency in MrBayes block
sed -i 's/ngen=[0-9]\+/ngen=100000/' $MB_BLOCK
sed -i 's/samplefreq=[0-9]\+/samplefreq=100/' $MB_BLOCK

# Run mdl using a block size of 10 (= 3240 parsimony calculations)
../scripts/mdl.pl $DATASET -b $MDL_BLOCK_SIZE -o mdl-example

# Run MrBayes on the 19 resulting blocks
../scripts/mb.pl mdl-example/$DATASET_ROOT.tar.gz -m $MB_BLOCK -o mb-example

# Check that the MCMC chains are below the desired threshold (should be all except one)
../scripts/mb.pl mb-example -c $CONVERGENCE_THRESHOLD

# Remove runs below the desired threshold (should be just the one)
../scripts/mb.pl mb-example -r $CONVERGENCE_THRESHOLD

# Rerun the block which fell below the threshold with more generations
sed -i 's/ngen=[0-9]\+/ngen=200000/' $MB_BLOCK
sed -i 's/samplefreq=[0-9]\+/samplefreq=200/' $MB_BLOCK
../scripts/mb.pl mdl-example/$DATASET_ROOT.tar.gz -m $MB_BLOCK -o mb-example

# Check that the MCMC chains are below the desired threshold (should be all now)
../scripts/mb.pl mb-example -c $CONVERGENCE_THRESHOLD

# Remove runs below the desired threshold (shouldn't be any)
../scripts/mb.pl mb-example -r $CONVERGENCE_THRESHOLD

# Run BUCKy on all 210 possible 4-taxon sets
../scripts/bucky.pl mb-example/$DATASET_ROOT.mb.tar -o bucky-example

# Create population tree using Quartet Max Cut
../scripts/get-pop-tree.pl bucky-example/$DATASET_ROOT.CFs.csv

# Determine branch lengths for population tree and draw it
cp bucky-example/$DATASET_ROOT.CFs.csv .
Rscript ../scripts/getTreeBranchLengths.r $DATASET_ROOT

# Run TICR
Rscript ../scripts/TICR.r $DATASET_ROOT

END=$(date +%s)

TIME=`echo "scale=2; ($END - $START) / 60" | bc | sed 's/^\./0./'`

echo -e "\nTotal execution time: $TIME minute(s)."
