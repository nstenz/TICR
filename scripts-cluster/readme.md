Goal: run the MrBayes and BUCKy steps (done by `mb.pl` and `bucky.pl`)
on a cluster that uses the job scheduler SLURM.

## mb.pl

Instead of the `mb.pl` script that parallelizes the MrBayes runs per gene, we want to use SLURM to parallelize the work. We have all the nexus files in the same folder, along with the following scripts:

- The julia script [`paste-mb-block.jl`](paste-mb-block.jl) takes a directory as argument.
  It will read all the nexus files in that directory (that is, files ending with ".nex").
  It will also read a text file with the MrBayes block, to paste onto all nexus files.
  This is similar to the option `‑m, ‑‑mb‑block` in the original `mb.pl` script,
  but the name of this translate file is hard-coded into the script: modify it to your needs,
  near the top of the file.
  A new directory will be created, with the new nexus files in it, and with loci renamed:
  `1.nex`, `2.nex`, ... A translation table (gene ID - gene name) will also be created
  (the name for this translation table is also hard-coded: modify it to your needs).

- The submit script [`mb-slurm-submit.sh`](mb-slurm-submit.sh) will parallelize
  all the individual-gene MrBayes runs with SLURM. Change `--array` to the correct number of genes.


In the original pipeline, `bucky.pl` runs `mbsum` on each gene,
and then `bucky` (using all genes) for each quartets.
Here we have two separate steps, explained below. 

## mbsum.pl

There is no such script yet, because this step is fast and the shell can suffice.
Suppose that we have a directory where each file is of the form `*.t` and is a MrBayes output file.
We can use `mbsum` to summarize a single `.t` file that correspond to one gene,
removing the first 1000 trees of each for burnin.

```shell
for X in *.t; do mbsum -n 1000 $X; done
```
This would create a file named <filename>.in for each file named <filename>.t .
Warning! It would overwrite files with the name <filename>.in if they exist.

Alternatively, the mbsum step can also be done with
the julia script [`mbsum-t-files.jl`](mbsum-t-files.jl).
This script takes as argument the folder name where all the tree files are located.  
**Warning**: a burnin of 2500 generations is hard coded: this can easily be changed,
near the top of the file.

## bucky-slurm.pl 

The perl script [`bucky-slurm.pl`](bucky-slurm.pl) will run `bucky`
**just once**, for a single 4-taxon set (or quartet) and it will take as input:
  - mbsum folder name: containing all files created by `mbsum`,
    one per locus. These files need to be named `*.in` to be used.
  - output name: `-o` or `--out-dir` name of the directory to store output files in
  - bucky arguments: `-a` or `--alpha` for the prior alpha value,
    and `-n` or `--ngen` number of generations
  - integer for the given quartet

It will produce, as output, a file for the given 4-taxon set with the
`.concordance` file and a file with the parsed output `.cf`
in the form of the CF table (to append later).

The perl script [`bucky-slurm.pl`](scripts-cluster/bucky-slurm.pl) can be run by SLURM
with the submit script [`bucky-slurm-submit.sh`](scripts-cluster/bucky-slurm-submit.sh).
Change `--array` to the appropriate number of 4-taxon sets for your data.

**Note** the `bucky` executable needs to be placed in `/workspace/software/bin`.
Adjust the path as appropriate for your cluster.

For each 4-taxon set, the parsed output file `.cf` contains a
single line with data for this 4-taxon set, in this order:

```
"taxon1,taxon2,taxon3,taxon4,CF12_34,CF12_34_lo,CF12_34_hi,CF13_24,CF13_24_lo,CF13_24_hi,CF14_23,CF14_23_lo,CF14_23_hi,ngenes"
```

The line above should serve as header before concatenating
all the individual parsed output files across all 4-taxon sets,
with something like `cat bucky/*.cf >> CFtable.csv` .
