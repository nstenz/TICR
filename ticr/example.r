# Cecile Ane, September 2014
# test of ILS on a population tree, and search for partially resolved tree
# needs:
# 1. concordance factors on many (or all) 4-taxon sets
# 2. guide tree

#------------- Preparation ------------------------------------------#
# upload ape library and necessary functions:
library(ape)
source("testingTreeWithQuartetCF.r")

# read the data and tree.
# First concordance factors: one line per quartet, 3 CFs per line. 
# columns 1-4 should have these headers:  taxon1 taxon2 taxon3 taxon4
# columns 5-7 should contain the CFs of quartets: 12|34, 13|24 and 14|23.
dat = read.csv("quartetCF.csv")
names(dat)[5:7] = c("cf12","cf13","cf14")
# Second: read the tree
guidetree = read.tree("tree.tre")

# prelimary calculations: may take a little while, like 19.5 seconds on the test example
prelim = test.tree.preparation(dat,guidetree)
Ntaxa = length(guidetree$tip.label)
# indices of internal edges: those that might be collapsed in a partially resolved tree
internal.edges = which(guidetree$edge[,2] > Ntaxa)


#---------- Test of Panmixia ----------------------------------------#

panmixia <- test.one.species.tree(dat,guidetree,prelim,edge.keep=NULL) # 0.10 seconds
panmixia$alpha        # 19.86
panmixia$minus.pll    # -51737.76 = negative log pseudo-likelihood
panmixia$X2           # 685.498 --> very large, so the full tree is rejected
unlist(panmixia[4:6]) # proportions used by the chi-square test.
# p.01..01 p.05..05 p.10..10 
#      424      401      956 

#---------- Test of the fully-resolved Tree -------------------------#
fulltree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=internal.edges) # 0.186 seconds
fulltree$alpha        # 165.9
fulltree$minus.pll    # -109566.8 = negative log pseudo-likelihood
fulltree$X2           # 149.9901 --> very large, so the full tree is rejected
unlist(fulltree[4:6]) # proportions used by the chi-square test.
# p.01..01 p.05..05 p.10..10
#      440      928     1195 

#---------- Search for partially-resolved Tree ----------------------#

# forward + backward search starting from the full tree:
resF <- stepwise.test.tree(dat,guidetree,search="both", kbest=15, 
                           maxiter=100, startT="fulltree") # 52.9 seconds
# At each step, the best choice was to removed an edge.
resF
# Nedge=21 edges kept:
# 1  2  4  6  7  8 11 14 20 21 23 24 31 34 35 36 38 39 44 47 53
# 6 edges not included: 3  5 19 22 37 48
# alpha = 173.7946
# negPseudoLoglik = -111005.3
# X2 = 298.7756 --> still large, so the partially-resolved tree is still rejected

# to recover the test on the partially-resolved tree without doing the search all over:
partialTree <- test.one.species.tree(dat,guidetree,prelim,
                    edge.keep=c(1,2,4,6,7,8,11,14,20,21,23,24,31,34,35,36,38,39,44,47,53))

# To start the search from panmixia:
resP <- stepwise.test.tree(dat,guidetree,search="both", kbest=15,
                           maxiter=100, startT="panmixia") # 96.2 seconds
# At each step, the best choice was to add an edge.
# The final tree is the same as when starting from the full tree:
resP$Nedge       # still 21
sort(resP$edges) # the same 21 edges are kept
sort(resP$notincluded) # the same 6 edges are collapsed.
