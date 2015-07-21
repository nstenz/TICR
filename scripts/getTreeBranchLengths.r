#!/usr/bin/Rscript

# input:  quartet CFs in csv file, and tree in nexus file.
# output: same tree topology, but branch lengths in coalescent units
#         u = -log(3/2 * (1 - mean CF of all quartets defining the branch))
#         u = 0.1 arbitrarily on external edges.
# Warning: the input tree needs to be binary, fully resolved.
# Depends on R package 'ape'
# Execute from terminal like this:
# Rscript --vanilla getTreeBranchLengths.r
# or like this to change the input file names, and outgroup for plotting:
# Rscript --vanilla getTreeBranchLengths.r new_file_name_root
# Rscript --vanilla getTreeBranchLengths.r new_file_name_root outgroup_name
# example:
# Rscript --vanilla getTreeBranchLengths.r chr4-subset A_Lyr

# 1. reads the tree, extract all quartets defining each branch,
#    writes them in file 'edge2quartets.txt'
# 2. reads in quartet CFs from data file. see format requirements below.
#    calculates average CF of quartets on each branch, the
#    associated coalescent units, writes them in csv file.
# 3. write these coalescent edge lengths into the tree, outputs
#    tree in newick format and, if desired,
#    plots the tree (pdf), edges annotated with branch lengths and mean CF.

#------------------------------------------------------------------#
# prelim: adapt file names as necessary                            #
#------------------------------------------------------------------#

filename.root <- "BD_walk_Q30_all.v1.reduced"
outgroup <- NULL
args <- commandArgs(TRUE)
if (length(args)>0)
  filename.root <- args[1]
if (length(args)>1)
  outgroup <- args[2]
cat("root for file names:",filename.root,"\n")
cat("outgroup:", ifelse(is.null(outgroup),"last taxon",outgroup),"\n")
tree.filename <- paste(filename.root, ".QMC.tre", sep="")
buckyCF.filename <- paste(filename.root, ".CFs.csv", sep="")
# required format looks like this: specific column names,
#
# taxon1 taxon2 taxon3 taxon4 CF12.34 CF13.24 CF14.23
#  A_Lyr Bsch_0 Da1_12  Dra_0   0.987  0.0065  0.0065
#  A_Lyr Bsch_0   Co_1  Hau_0   0.265  0.3430  0.3920
#   ...                                          ...
#  Uod_1 Vind_1   Wt_5   Yo_0   0.272   0.460   0.267
#  Van_0 Vind_1   Wt_5   Yo_0   0.425   0.293   0.282
#  Van_0 Vind_1   Wa_1   Wt_5   0.482   0.269   0.249
edgeLengths_meanCF.filename <- paste(filename.root,".CFbyEdge.csv",sep="")
tree.withbranchlengths.filename <- paste(filename.root,".QMClengths.tre",sep="")
drawTree.withBranchLengths <- TRUE
tree.pdf.filename <- paste(filename.root,".QMClengths.pdf",sep="")
cf.threshold <- 0.38 # edges with CF above this value will be drawn thicker

#------------------------------------------------------------------#
# step 1: read tree and extract all quartets defining each branch  #
#------------------------------------------------------------------#

library(ape)
tre <- read.tree(tree.filename)
ntax <- length(tre$tip.label) # number of taxa
cat("\ttree was read.",ntax,"taxa.\n")
# arbitrarily root with last taxon by default. For technical reasons
# only, because the tree is considered as unrooted.
if (is.null(outgroup)){  outg <- ntax
} else {
 outg <- which(tre$tip.label == outgroup)
 if (length(outg) == 0)
   stop("outgroup taxon ",outgroup," was not found in the tree\n")
 if (length(outg) > 1)
   stop("outgroup taxon ",outgroup," found more than once in the tree\n")
}
tre <- root(tre,outgroup=outg,resolve.root=FALSE)
nedg <- dim(tre$edge)[1] # total number of edges
intEdge <- which(tre$edge[,2]>ntax) # indices of internal edges
outEdge <- which(tre$edge[,2] == outg) # external edge to single outgroup

# Make a list of all descendants from each edge
edgeDescendants <- vector("list",nedg)
for (i in nedg:1){
 childnode <- tre$edge[i,2]
 if (childnode <= ntax){ # external edge, childnode = leaf
  edgeDescendants[[i]] <- tre$tip.label[childnode]
 } else { # internal edge
  tmp <- which(tre$edge[,1]==childnode) # indices of the 2 children edges
  edgeDescendants[[i]] <- c(edgeDescendants[[tmp[1]]],edgeDescendants[[tmp[2]]])
 }
}
# edgeDescendants[intEdge]

# Make a list of 4-way partitions, one for each edge
part1=vector("list",nedg) # descendants from child 1 
part2=vector("list",nedg) # descendants from child 2
part3=vector("list",nedg) # descendants from sibling
part4=vector("list",nedg) # all remaining tata = complement of descendants from parent

for (i in 1:nedg){
  childnode = tre$edge[i,2] # if child node = leaf: do nothing
  if (childnode > ntax){
    child = which(tre$edge[,1]==childnode) # extracted 2 children edges
    part1[[i]] = edgeDescendants[[child[1]]]
    part2[[i]] = edgeDescendants[[child[2]]]
    sibling = which(tre$edge[,1]==tre$edge[i,1])
    # vector of both the edge and its sibling(s). Only one sibling
    # unless the edge connects to the root, in which case there are 2 siblings:
    # the outgroup and the true sibling edge
    # (which might be internal or external).
    if (length(sibling) == 3){ # to identify root case and remove false sibling
      out <- which(sibling == outEdge)
      sibling <- sibling[-out]
    }
    if (sibling[1]==i){ sibling=sibling[2] } else {sibling=sibling[1]}
    part3[[i]] <- edgeDescendants[[sibling]]
    part4[[i]] <- setdiff(tre$tip.label, c(part1[[i]],part2[[i]],part3[[i]]))
  }
}
# i=1; part1[i]; part2[i]; part3[i]; part4[i]

# Make a list of quartets now, for internal edges only
Nquartets <- sapply(part1[intEdge],length)*sapply(part2[intEdge],length) *
             sapply(part3[intEdge],length)*sapply(part4[intEdge],length)
# Nquartets is the number of quartets per edge
# sum(Nquartets) # 8780 (cp) or 7216 (pda): much fewers than choose(30,4)=27405
dat <- data.frame(edge   = rep(intEdge,Nquartets),
                  taxon1 = rep(NA,sum(Nquartets)),
                  taxon2 = NA, taxon3=NA, taxon4=NA,
                  quartet= NA, CF=NA
                  )
j <- 0 # indexes the row
for (i in intEdge){
 for (t1 in part1[[i]]){ for (t2 in part2[[i]]){ for (t3 in part3[[i]]){ for (t4 in part4[[i]]){
  j <- j+1
  o <- order(c(t1,t2,t3,t4))
  r <- rank( c(t1,t2,t3,t4))
  quartet <- c(t1,t2,t3,t4)[o]
  dat$taxon1[j] <- quartet[1]
  dat$taxon2[j] <- quartet[2]
  dat$taxon3[j] <- quartet[3]
  dat$taxon4[j] <- quartet[4]
  # the quartet on this edge is t1 t2 | t3 t4, but written 
  # using indices when the 4 taxa are sorted alphabetically. 
  if (which(r==1) >=3 ){ # in this case r[3]==1 or r[4]==1
   left  <- sort(r[3:4])
   right <- sort(r[1:2])
  } else { # in this case r[1]==1 or r[2]==1
   left  <- sort(r[1:2])
   right <- sort(r[3:4])
  }
  dat$quartet[j] <- paste(left[1],left[2],"|",right[1],right[2],sep="")
 }}}}
}
# head(dat)
# subset(dat,edge==1)
# subset(dat,edge==31)
# write.table(dat,"edge2quartets.txt",row.names=F,quote=F)
Nqrt <- dim(dat)[1] # number of 4-taxon sets whose internal branch covers a single edge in the tree.
# This is much smaller than the total number of 4-taxon sets: choose(ntax,4)
# (8780 for stratified subsample, 7216 for PD, compared to choose(30,4)=27,405 4-taxon sets total)
cat("\tlisted the",Nqrt,"quartets associated with edges in tree.\n")

#------------------------------------------------------------------#
# step 2: get mean CF for internal edges, and coalescent lengths   #
#------------------------------------------------------------------#

#dat <- read.table(edge2quartet.filename, header=T)
#for (i in 2:5){ dat[,i] = as.character(dat[,i]) }
cf <- read.csv(buckyCF.filename, header=T)
if (dim(cf)[1] != choose(ntax,4)) # 27405 with 30 taxa
  cat("Warning: it looks like there are missing (or extra) quartets in\n\t",
      buckyCF.filename,".\n\tRead ",dim(cf)[1]," rows, expected ",
      choose(ntax,4)," quartets.\n",sep="")
# create an extra matrix and columns "t1" ... "t4" to identify
# the 4-taxon set uniquely: 4 taxon names sorted alphabetically
cf$taxon1 <- as.character(cf$taxon1); cf$t1 <- NA
cf$taxon2 <- as.character(cf$taxon2); cf$t2 <- NA
cf$taxon3 <- as.character(cf$taxon3); cf$t3 <- NA
cf$taxon4 <- as.character(cf$taxon4); cf$t4 <- NA
taxonOrder <- matrix(NA,dim(cf)[1],4)
for (i in 1:dim(cf)[1]) taxonOrder[i,] <- order(cf[i,1:4])
for (i in 1:4){
  for (j in 1:4){
    ind <- taxonOrder[,i]==j
    cf[ind,paste0("t",i)] <- cf[ind,paste0("taxon",j)]
  }
}
# extract CF value for the quartet resolution in the tree
for (i in 1:Nqrt){
 ind <- NA
 ind <- which(cf$t1==dat$taxon1[i] & cf$t2==dat$taxon2[i] &
              cf$t3==dat$taxon3[i] & cf$t4==dat$taxon4[i])
 if (length(ind)>1)
   stop(paste("found multiple rows in CF file for these taxa:",
              dat$taxon1[i],dat$taxon2[i],dat$taxon3[i],dat$taxon4[i]))
 if (length(ind)==0)
   next
 # resolution needs to be "CF12.34" or "CF13.24" or "CF14.23"
 #            with numbers corresponding to ordering of taxa in cf.
 # The line below would work turn 13|24 into CF13.24 or 12|34 into CF12.34,
 # but assumes that the taxa are alphabetically ordered in cf. 
 # resolution <- paste("CF",sub('\\|','.',dat$quartet[i]),sep="")
 
 # The more complicated solution below allows for any ordering of
 # the 4 taxa ('taxon1', ... , 'taxon4') in cf. 
 sister2one <- as.numeric(substr(dat$quartet[i], 2,2)) # taxoni sister to taxon1, in dat
 leftclade <- sort(taxonOrder[ind,c(1,sister2one)])    # in cf
 rightclade <- setdiff(1:4,leftclade) # will be sorted
 if (leftclade[1]==1){
   resolution <- paste0("CF1",leftclade[2],".",rightclade[1],rightclade[2])
 } else
   resolution <- paste0("CF1",rightclade[2],".",leftclade[1],leftclade[2])
 dat$CF[i] <- cf[ind,resolution]
}
ind <- which(is.na(dat$CF)) # to identify quartet with missing CF data
if (length(ind)>0)
  cat(length(ind),"quartets defining edges in the tree have missing CF data\n")

newdat <- data.frame(edge      = unique(dat$edge),
                     meanCF    = tapply(dat$CF, dat$edge, mean, na.rm=T),
                     sdCF      = tapply(dat$CF, dat$edge, sd,   na.rm=T),
                     Nquartets = as.numeric(table(dat$edge))
                     )
newdat$edge.length = -log(3/2 * (1-newdat$meanCF))
# newdat; hist(newdat$meanCF, breaks=30, col="tan")
write.csv(newdat,edgeLengths_meanCF.filename,quote=F,row.names=F)
cat("\tlisted internal edges with mean (and sd) of CFs of quartets\n\t  on each edge, and associated coalescent units.\n\t  See in",edgeLengths_meanCF.filename,"\n")

#------------------------------------------------------------------#
# step 3: write edge lengths into the tree, output tree, and plot  #
#------------------------------------------------------------------#

tre$edge.length <- rep(NA,nedg)
extEdge <- which(tre$edge[,2]<=ntax) # 27 internal edges, 30 external
tre$edge.length[extEdge] <- 0.1 # this is arbitrary: no meaning
tre$edge.length[newdat$edge] <- newdat$edge.length
tre$edge.length[tre$edge.length<0] <- 0
# might happen if the edge is contradicted by data, i.e. if mean CF <1/3
tre$edge.length[tre$edge.length==Inf] <- 30 # arbitrary large
# happens if the edge has CF 1.0, corresponding to infinite # coalescent units.
# replacing by 30. exp(-30) is close enough to 0 discordance.

write.tree(tre,tree.withbranchlengths.filename)
cat("\twrote tree with branch lengths (coalescent units) to\n\t ",tree.withbranchlengths.filename,"\n")

if (drawTree.withBranchLengths){
  ew <- rep(1,nedg) # used for edge width
  ew[newdat$edge[newdat$meanCF >= cf.threshold]] <- 3
  pdf(tree.pdf.filename,height=10,width=14)
  plot(tre,no.margin=T,edge.width=ew,use.edge.length=F,show.tip.label=F)
  edgelabels(round(tre$edge.length[newdat$edge],2),newdat$edge,
             frame="n",adj=c(.5,-0.5),cex=.8)
  edgelabels(round(newdat$meanCF,2),newdat$edge,frame="n",adj=c(.5,1.3),font=3,cex=.8)
  tiplabels(tre$tip.label, frame=c("none"), bg=NULL, adj=c(0,.5))
  dev.off()
  cat("\tplotted tree in",tree.pdf.filename,"\n")
  cat("\t above branches: coalescent units\n\t below branches: average quartet concordance factors.\n")
}
