# run this from the command line with:
# Rscript --vanilla TICRtest.r file_name_root
# example:
# Rscript --vanilla TICRtest.r chr4-subset
# Requirements:
# - R package "ape"
# - TICR functions in file testingTreeWithQuartetCF.r, placed in ../ticr/
# - quartet CF data: csv file with these column names in the header:
#   "taxon1",...,"taxon4","CF12.34","CF13.24","CF14.23" 

library(ape)
source("../ticr/testingTreeWithQuartetCF.r")

filename.root <- "chr4-subset"
args <- commandArgs(TRUE)
if (length(args)>0)
  filename.root <- args[1]
cat("root for file names:",filename.root,"\n")

tree.filename <- paste0(filename.root, ".QMClengths.tre")
buckyCF.filename <- paste0(filename.root, ".CFs.csv")
output.filename <- paste0(filename.root, ".ticr.txt")
partialtree.pdf.filename <- paste0(filename.root, ".ticr.pdf")
cat("output will be redirected to",output.filename,"\n")
cat("                      and to",partialtree.pdf.filename,"\n")

zz <- file(output.filename, open="wt")
sink(zz)
sink(zz, type="message")

tmp <- c(paste0("taxon",1:4),"CF12.34","CF13.24","CF14.23")
# tmp = names of columns to extract
dat <- read.csv(buckyCF.filename)[,tmp]
if (ncol(dat) != 7)
  stop("column headers not found in CF data (",buckyCF.filename,")\n")
guidetree <- read.tree(tree.filename)

# prelimary calculations: may take a little while, like 19.5 seconds on the test example
prelim <- test.tree.preparation(dat,guidetree)
Ntaxa <- length(guidetree$tip.label)
Nquartets <- dim(dat)[1] # number of 4-taxon sets with CF data
# indices of internal edges: those that might be collapsed in a partially resolved tree
internal.edges <- which(guidetree$edge[,2] > Ntaxa)


#---------- Test of Panmixia ----------------------------------------#
cat("\n#----------------------------------------------#")
cat("\nTest of panmixia:\n\n")
panmixia <- test.one.species.tree(dat,guidetree,prelim,edge.keep=NULL)
panmixia[1:4] # gives alpha, negative log pseudo-likelihood,
              # X2 and p-value from chi-square test of panmixia
panmixia$outlier.table # table of number of quartets with low / large p-values
cat("\n(Warnings may appear from the chi-square test having fewer than
5 expected 4-taxon sets in some outlier categories (e.g. p<.01))\n")

#---------- Test of the fully-resolved Tree -------------------------#
cat("\n#----------------------------------------------#")
cat("\nTest of binary tree:\n\n")
fulltree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=internal.edges) # 0.186 seconds
fulltree[1:4] # gives alpha, negative log pseudo-likelihood,
              # X2 and p-value from chi-square test of panmixia
fulltree$outlier.table # table of number of quartets with low / large p-values

#---------- Search for partially-resolved Tree ----------------------#
# forward + backward search starting from panmixia:
cat("\n#----------------------------------------------#")
cat("\nSearch for partial tree, starting from panmixia:\n\n")
resP <- stepwise.test.tree(dat,guidetree,search="both", kbest=15,
                           maxiter=100, startT="panmixia")
resP

# forward + backward search starting from the full tree:
cat("\nSearch for partial tree, starting from binary tree:\n\n")
resF <- stepwise.test.tree(dat,guidetree,search="both", kbest=15, 
                           maxiter=100, startT="fulltree")
resF

cat("\nPlot of partially resolved tree, branch lengths in coalescent units:\n")
cat("\tSee output in TICRtest.pdf\n\n")
pdf(file=partialtree.pdf.filename, width=5,height=5)
plot.species.tree(guidetree,resF$edges)
add.scale.bar()
dev.off()

# to recover the test on the partially-resolved tree without doing the search all over:
cat("\nTest of partially resolved tree:\n\n")
partialTree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=resF$edges)
partialTree[1:4]
partialTree$outlier.table
