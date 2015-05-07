# Cecile Ane, September 2014 - May 2015
# test of ILS on a population tree, and search for partially resolved tree
# needs:
# 1. concordance factors on many (or all) 4-taxon sets
# 2. guide tree
# 3. ape package already installed in R

#------------- Preparation ------------------------------------------#
#              upload ape library and necessary functions            #

library(ape)
source("testingTreeWithQuartetCF.r")

# read the CF data:
# First concordance factors: one line per quartet, 3 CFs per line. 
# columns 1-4 should have these headers:  taxon1 taxon2 taxon3 taxon4
# columns 5-7 should contain the CFs of quartets: 12|34, 13|24 and 14|23.
dat = read.csv("quartetCF.csv")
# read the guide tree:
guidetree = read.tree("tree.tre")

# prelimary calculations: may take a little while, like 19 seconds on the test example
prelim = test.tree.preparation(dat,guidetree)
Ntaxa = length(guidetree$tip.label)
# indices of internal edges: those that might be collapsed in a partially resolved tree
internal.edges = which(guidetree$edge[,2] > Ntaxa)


#---------- Test of Panmixia ----------------------------------------#

panmixia <- test.one.species.tree(dat,guidetree,prelim,edge.keep=NULL) # 0.10 seconds
panmixia[1:4] # gives alpha=19.86, negative log pseudo-likelihood=-51737.76,
              # X2=685.498 and p-value=3e-148 from chi-square test of panmixia
panmixia$outlier.table # table of number of quartets with low / large p-values
#             .01    .05     .10   large
# observed 424.00  401.0  956.00 25624.0
# expected 274.05 1096.2 1370.25 24664.5
# example: we expect 274 4-taxon sets to have a p-value<0.01 just by chance,
#          but we actually observe 401 such 4-taxon sets

#---------- Test of the fully-resolved Tree -------------------------#

fulltree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=internal.edges) # 0.186 seconds
fulltree[1:4] # alpha=165.9, negative log pseudo-likelihood=-109566.8,
              # X2=149.9901 and p-value=2.6e-32 from chi-square test of binary tree
fulltree$outlier.table
#             .01    .05     .10   large
# observed 440.00  928.0 1195.00 24842.0
# expected 274.05 1096.2 1370.25 24664.5

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
# X2 = 298.7756, chisq.pval=1.8e-64 --> the partially-resolved tree is still rejected

# To start the search from panmixia:
resP <- stepwise.test.tree(dat,guidetree,search="both", kbest=15,
                           maxiter=100, startT="panmixia") # 96.2 seconds
# At each step, the best choice was to add an edge.
# Same partial tree with same 21 resolved edges:
resP
 
# to see and re-analyze the partially-resolved tree without doing the search all over:
edges2keep <- c(1,2,4,6,7,8,11,14,20,21,23,24,31,34,35,36,38,39,44,47,53)
# or just: edges2keep <- resF$edges
plot.species.tree(guidetree,edges2keep)
partialTree <- test.one.species.tree(dat,guidetree,prelim,edge.keep=edges2keep)
partialTree[1:4]
partialTree$outlier.table

#-------- Identify taxa most responsible for extra outlier quartets -----#

outlier.4taxa <- which(partialTree$outlier.pvalues < 0.01)
length(outlier.4taxa) # 483 40-taxon sets with outlier p-value below 0.01
q01 = as.matrix(dat[outlier.4taxa,1:4])
sort(table(as.vector(q01)),decreasing=T)
# Cnt_1 Vind_1  Hau_0 Ragl_1  Mnz_0   Nw_0 Qar_8a   Co_1  A_Lyr Bsch_0   Jm_0 
#   239    239    126    123    120     92     74     66     65     55     53 
#  Is_0  Pla_0 Rome_1  Kro_0  Tha_1   Tu_0   Wa_1 Pna_17   Uk_1   En_2   Et_0 
#    52     51     51     45     44     42     39     37     37     36     36 
# Uod_1  Hey_1   Wt_5 Da1_12  Dra_0   Ha_0   Yo_0  Van_0 
#    34     32     27     26     26     24     21     20 
# So: Cnt_1 appears in 239 of these 483 outlier 4-taxon sets.

sum(apply(q01,1,function(x){"Cnt_1" %in% x | "Vind_1" %in% x}))
# 266 outlier 4-taxon sets have either Cnt_1 or Vind_1
sum(apply(q01,1,function(x){"Cnt_1" %in% x & "Vind_1" %in% x}))
# 212 outlier 4-taxon sets have both Cnt_1 or Vind_1
