#!/bin/bash

DATASET_ROOT="chr4-subset"
DATASET="$DATASET_ROOT.nex"

# Rename example datafile so we don't delete it
mv $DATASET tmp.nex

# Delete script output
rm -rf *-example/
rm -f $DATASET_ROOT*

# Rename example datafile back to its original name
mv tmp.nex $DATASET
