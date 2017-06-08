# run this from the command line with:
# Rscript --vanilla TICR.r file_name_root
# example:
# Rscript --vanilla TICR.r chr4-subset
# Requirements:
# - R package "ape"
# - TICR functions in file testingTreeWithQuartetCF.r, placed in ../ticr/
# - quartet CF data: csv file with these column names in the header:
#   "taxon1",...,"taxon4","CF12.34","CF13.24","CF14.23"
# Claudia: save at the end an Rda file with most important objects

library(ape)
source("testingTreeWithQuartetCF.r")

filename.root <- "chr4-subset"
args <- commandArgs(TRUE)
if (length(args)>0)
  filename.root <- args[1]
cat("root for file names:",filename.root,"\n")

tree.filename <- paste0(filename.root, ".QMClengths.tre")
buckyCF.filename <- paste0(filename.root, ".CFs.csv")
output.filename <- paste0(filename.root, ".ticr.txt")
rda.filename <- paste0(filename.root, ".Rda")
partialtree.pdf.filename <- paste0(filename.root, ".ticr.pdf")
cat("output will be redirected to",output.filename,"\n")
cat("                      and to",partialtree.pdf.filename,"\n")

zz <- file(output.filename, open="wt")
sink(zz)
sink(zz, type="message")

tmp <- c(paste0("taxon",1:4),"CF12_34","CF13_24","CF14_23")
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

cat("\nPlots of results for each test (panmixia, full tree, partial tree)")
#cat("\n\t(branch lengths in coalescent units in trees):")
cat("\n\tSee output in",partialtree.pdf.filename,"\n\n")
pdf(file=partialtree.pdf.filename, width=10.5,height=8)
par(mar=c(3.2,3.1,2.5,.5), tck=-0.01, mgp=c(2.2,.4,0), las=0, oma=c(0,0,2.5,0))


#---------- Test of Panmixia ----------------------------------------#
cat("\n#----------------------------------------------#")
cat("\nTest of panmixia:\n\n")
panmixia <- test.one.species.tree(dat,guidetree,prelim,edge.keep=NULL)
panmixia[1:4] # gives alpha, negative log pseudo-likelihood,
              # X2 and p-value from chi-square test of panmixia
panmixia$outlier.table # table of number of quartets with low / large p-values
cat("\n(Warnings may appear from the chi-square test having fewer than
5 expected 4-taxon sets in some outlier categories (e.g. p<.01))\n")
mtext(paste0("Test of panmixia\nalpha=",signif(panmixia$alpha,2),
             ", X2=",signif(panmixia$X2,2),", X2 pval=",signif(panmixia$chisq.pval,2)),
      side=3, adj=0.01, outer=T, line=-1)


#---------- Test of the fully-resolved Tree -------------------------#
cat("\n#----------------------------------------------#")
cat("\nTest of binary tree:\n\n")
fulltree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=internal.edges) # 0.186 seconds
fulltree[1:4] # gives alpha, negative log pseudo-likelihood,
              # X2 and p-value from chi-square test of panmixia
fulltree$outlier.table # table of number of quartets with low / large p-values
mtext(paste0("Test of full tree\nalpha=",round(fulltree$alpha,2),
             ", X2=",round(fulltree$X2,2),", X2 pval=",signif(fulltree$chisq.pval,3)),
      side=3, adj=0.01, outer=T, line=-1)

#---------- Search for partially-resolved Tree ----------------------#
# forward + backward search starting from panmixia:
cat("\n#----------------------------------------------#")
cat("\nSearch for partial tree, starting from panmixia:\n\n")
resP <- stepwise.test.tree(dat,guidetree,search="both", kbest=15,
                           maxiter=100, startT="panmixia")
resP[1:7]
resP$outlier.table
mtext(paste0("Test of partial tree (search started from panmixia)\nalpha=",
             round(resP$alpha,2),", X2=",round(resP$X2,2),", X2 pval=",signif(resP$chisq.pval,3)),
      side=3, adj=0.01, outer=T, line=-1)


# forward + backward search starting from the full tree:
cat("\nSearch for partial tree, starting from binary tree:\n\n")
resF <- stepwise.test.tree(dat,guidetree,search="both", kbest=15,
                           maxiter=100, startT="fulltree")
resF[1:7]
resF$outlier.table
mtext(paste0("Test of partial tree (search started from full tree)\nalpha=",
             round(resF$alpha,2),", X2=",round(resF$X2,2),", X2 pval=",signif(resF$chisq.pval,3)),
      side=3, adj=0.01, outer=T, line=-1)

save(panmixia, fulltree, resP, resF, guidetree, file=rda.filename)

dev.off()

# to recover the test on the partially-resolved tree without doing the search all over:
#cat("\nTest of partially resolved tree:\n\n")
#partialTree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=resF$edges)
#partialTree[1:4]
#partialTree$outlier.table
