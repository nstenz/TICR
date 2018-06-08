Goal: Want to modify `mb.pl` and `bucky.pl` for SLURM.

## mb.pl
Instead of the `mb.pl` script that parallelizes the MrBayes runs per gene, we want to use SLURM to parallelize the work. We have all the nexus files in the same folder, along with the following scripts:
- The script `paste-mb-block.jl` will read all the nexus files in the directory, and will read a textfile with the MrBayes block to paste on all nexus files (option `‑m, ‑‑mb‑block` in the original `mb.pl` script). All the nexus files will be renamed: `1.nex, 2.nex,....`
- The `mb-slurm-submit.sh` script will parallelize each gene run with SLURM. Change `--array` to the correct number of genes


## bucky.pl
Original `bucky.pl` runs `mbsum` and then `bucky` for all quartets. We want to have two separate scripts: 

### mbsum.pl
Not done yet, but suppose that you have a directory where each file is of the form *.t and is a MrBayes output file.
Use mbsum to summarize each file.  Remove the first 1000 trees of each for burnin.
```
  for X in *.t; do mbsum -n 1000 $X; done
```
This will create a file named <filename>.in for each file named <filename>.t .
Warning!  It will overwrite files with the name <filename>.in if they exist.

This can also be done with `mbsum-t-files.jl`. Warning: this script is hard-coded to 3 independent runs in MrBayes (but can be easily changed).

### bucky-slurm.pl 
This script will run **just one** bucky for a given quartet and it will take as input:
  - mbsum folder name
  - output name
  - bucky arguments: alpha,...
  - integer for the given quartet

The bucky script will provide as output a file for that quartet with the .concordance and a file with the parsed output in the form of the CF table (to append later).

It is run in SLURM with `bucky-slurm-submit.sh`. Change `--array` to the number of quartets.

**Note** Need to put the executable of bucky in `/workspace/software/bin`



