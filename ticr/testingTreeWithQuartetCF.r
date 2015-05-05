# TICR: Tree Incongruence Checking with R
# test of ILS on a population tree. 
# Cecile Ane, September 2014
# input: concordance factors on many (or all) 4-taxon set
#        guide tree

# See how these functions are used in file example.r

stepwise.test.tree = function(cf, guidetree, search="heuristic", method="PLL",
                              kbest=5, maxiter=100, startT="panmixia"){
# search: "heuristic" (method etc. ignored) or
#         "both" directions (method etc. are used)
# kbest: lower value for faster, less thorough search.
# startT: starting partial tree for stepwise search. One of:
#         "panmixia", "fulltree", or numeric vector of edge numbers.

 # change all external edges of the guide tree to 0
 external.edge = which(guidetree$edge[,2]<=length(guidetree$tip.label))
 internal.edge = which(guidetree$edge[,2]> length(guidetree$tip.label))
 guidetree$edge.length[external.edge]=0
 # this is modifying guidetree in the current environment only, not in gloval env.
 Ntip = length(guidetree$tip.label)
 tax2id = 1:Ntip; names(tax2id)=guidetree$tip.label
 M = nrow(cf) # 27405: number of quartets with data
 cf.obs = as.matrix(cf[,5:7])
 logCFsum = sum(log(cf.obs)) 
 dat.names = as.matrix(cf[,1:4]) # -> taxon names as characters, not factors
 dat.leaves = matrix(tax2id[dat.names],nrow=M,ncol=4)

 tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
 names(tab0)=c(".01",".05",".10","large") # tab0*M:  .01    .05     .10   large
                                          #       274.05 1096.2 1370.25 24664.5 
 cat("determining tree traversal post-order... ")
 nodes.postorder = rep(NA,Ntip+guidetree$Nnode)
 guidetree$node.order = rep(NA,Ntip+guidetree$Nnode)
 nextnode = 1
 set.node.postordertraversal = function(nodeid){
  # external: nextnode, nodes.postorder, guidetree$node.order
  if (nodeid > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==nodeid) # child edges
   for (ce in ces){
    set.node.postordertraversal(guidetree$edge[ce,2])
   }
  }
  guidetree$node.order[nodeid] <<- nextnode
  nodes.postorder[nextnode] <<- nodeid
  nextnode <<- nextnode+1
 }
 set.node.postordertraversal(Ntip+1) # ID for the root
 #plot(guidetree);tiplabels(adj=0);nodelabels(adj=1);edgelabels()
 cat("done.\n")

 cat("calculating matrix of descendant relationships... ")
 node2descendant = matrix(FALSE, Ntip+guidetree$Nnode, Ntip+guidetree$Nnode)
 # entry [i,j] is TRUE if j is a descendant of i.
 for (n in nodes.postorder){
  node2descendant[n,n] = T
  if (n > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==n) # child edges     
   for (cn in guidetree$edge[ces,2]){ # child nodes
    node2descendant[n,] = node2descendant[cn,] | node2descendant[n,]
   }
  }
 }
 cat("done.\n")

 get.mrca = function(nodes){
  n = max(guidetree$node.order[nodes]) # first candidate node
  for (n in nodes.postorder){
   if (all(node2descendant[n,nodes])){ return(n) }
  }
 }

 #guidetree.dist = round(cophenetic(guidetree),digits=12)
 #get.dominant = function(x){
 # quartet.dist = guidetree.dist[x,x]
 # threedist = c(quartet.dist[1,2]+quartet.dist[3,4],quartet.dist[1,3]+quartet.dist[2,4],
 #               quartet.dist[1,4]+quartet.dist[2,3])
 # return(order(threedist)[1]) #  index for which treedist is min
 #} #dominant = apply(dat.names,1,get.dominant)

 cat("calculating matrix of edges spanned by each quartet... ")
 dominant.table  = c(1:3,1:3) # ie tree displays 12|34 if dominant = 1
                              #                  13|24               2
                              #                  14|23               3
 dominant = rep(NA,M)
 quartet2edge = matrix(FALSE,M,length(guidetree$edge.length)) # big matrix!

 for (i in 1:M){
  x = dat.leaves[i,]
  ancestor = sapply(list(x[1:2],x[c(1,3)],x[c(1,4)],x[3:4],x[c(2,4)],x[2:3]),get.mrca)
  ind = which.min(guidetree$node.order[ancestor])
  a1 = ancestor[ind]
  dominant[i] = dominant.table[ind]
  mrca = ancestor[which.max(guidetree$node.order[ancestor])]
  ancestor[ind] = Ntip+1 # replacing by the root to find 2nd best, next
  a2 = ancestor[which.min(guidetree$node.order[ancestor])]
  mynode = a1
  while (mynode!=a2 & mynode!=mrca){
   parentE = which(guidetree$edge[,2]==mynode)
   quartet2edge[i,parentE]=TRUE
   mynode = guidetree$edge[parentE,1]
  }
  if (mynode == mrca){ # the mrca of all 4 taxa is on path between both ancestors
   mynode = a2
   while (mynode!=mrca){
    parentE = which(guidetree$edge[,2]==mynode)
    quartet2edge[i,parentE]=TRUE
    mynode = guidetree$edge[parentE,1]
   }
  }
 } # i=4288; dat.names[i,]; quartet2edge[i,]
 colnames(quartet2edge)=as.character(1:length(guidetree$edge.length))
 cat("done.\n")

 total.ew = rep(1,M) %*% quartet2edge[,internal.edge] # number quartet per edge
 qw =  1/quartet2edge %*% rep(1,dim(quartet2edge)[2])
    # quartet weight: 1/Ne for each edge, where Ne= # edges for that quartet
 qw.ew = t(qw) %*% quartet2edge[,internal.edge]
 score.edge = function(quartetID){
 # !Warning! external variables: cf, guidetree
 # gives score to each branch in guidetree: number/weight of 4-taxon sets 
 #       in 'quartetID' whose internal edge includes that branch
  ues = rep(1,length(quartetID)) %*% quartet2edge[quartetID,internal.edge] / total.ew
  wes = qw[quartetID] %*% quartet2edge[quartetID,internal.edge] / qw.ew
  return(rbind(ues,wes))
 }

 plot.species.tree = function(edge.keep){
  sp.tree = guidetree
  if (length(edge.keep)==0){
   warning("Complete panmixia: the tree looks like a single tip. Won't plot this.")
   return(NULL)
  }
  sp.tree$edge.length[-edge.keep]=0
  plot(sp.tree)
 }

 test.species.tree = function(edge.keep,plot=FALSE){
  # !Warning! external variables:  guidetree, tab0, M, cf.obs, dominant
  #----------- get expected concordance factors -----------------------#
  if (length(edge.keep)==0) { tk = numeric(M) } else { 
   tk = quartet2edge[,edge.keep] %*% cbind(guidetree$edge.length[edge.keep]) }
  cf.exp = matrix(exp(-tk)/3, M,3) # minor CFs only, so far
  # cf.exp[i,j] also = cf.exp[i+(j-1)*M] 
  cf.exp[1:M + (dominant-1)*M] = 1-exp(-tk)*2/3
  # expected CFs: q12.34, q13.24, q14.23
  # dat.names[1336,]; cf.exp[1336,]; tk[1336] # A_Lyr Bik_1 Stw_0 Ting_1, 0.484 0.269 0.247
  #--------------------- estimate alpha -------------------------#
  if (plot) { # fixit: find better visualization with color gradients
   plot(cf.obs,cf.exp); abline(a=0,b=1,col="red")
  }
  logCFtilde = sum(log(cf.obs)*cf.exp)/M 
  tmp = rle(sort(tk))
  tk.unique = tmp$values
  nk = tmp$lengths # number of quartets with internal edge of length tk.
  wk = nk/M  # weight for quartets along edge k. They sum to 1.
  pk.minor = exp(-tk.unique)/3
  pk.domin = 1-2*pk.minor
  # f1 = -logPL/M with constant term logCFbar ommitted
  # f2 = d(f1)/d(alpha)
  f1 = function(alpha){-lgamma(alpha) - alpha * logCFtilde +
                         sum(wk*lgamma(alpha*pk.domin)) +
                       2*sum(wk*lgamma(alpha*pk.minor))   }
  f2 = function(alpha){-digamma(alpha) - logCFtilde +
                         sum(wk*pk.domin*digamma(alpha*pk.domin)) +
                       2*sum(wk*pk.minor*digamma(alpha*pk.minor)) }
  # now need to solve f2=0, or minimize f1
  # sa = (2/9)*(3*M)/sum((cf.obs-cf.exp)^2) # starting alpha value
  # optim(243,f1,gr=f2, method="L-BFGS-B",lower=0,upper=+Inf) #248.2907
  res=optimize(f1,interval=c(0,10^5))
  alpha=res$minimum
  minus.pll = res$objective * M +logCFsum # - pseudo-log-likelihood
  # cat("estimated alpha:",alpha,", -pseudo log-lik:",minus.pll,"\n")
  # xx = seq(1,4*alpha,by=1); yy=sapply(xx,f1)
  # plot(xx,yy,type="l"); points(alpha,res$objective,col="red",pch=16)
  #--------------------- get a p-value for each quartet -------------------#
  p.exp = pmax(cf.exp[,1], cf.exp[,2],cf.exp[,3])
  ind.tree = which(tk>0)  
  ind.panmixia = which(tk==0)
  p.obs = pmax(cf.obs[,1], cf.obs[,2],cf.obs[,3]) # max for panmixia test
  # cf.obs[i+(j-1)*M] is same as cf.obs[i,j]
  p.obs[ind.tree] = cf.obs[ind.tree+(dominant[ind.tree]-1)*M] 

  pval = rep(NA,M) # p-values
  pval[ind.panmixia] = 3*pbeta(p.obs[ind.panmixia],lower.tail=F,
                                  shape1=alpha/3,shape2=alpha*2/3)
  dev = abs(p.exp[ind.tree]-p.obs[ind.tree]) # deviations
  pe = p.exp[ind.tree] # expected proportions
  pval[ind.tree] =
    pbeta(pe+dev,shape1=alpha*pe,shape2=alpha*(1-pe),lower.tail=F) +
    pbeta(pe-dev,shape1=alpha*pe,shape2=alpha*(1-pe),lower.tail=T)
  pval[pval>1] = 1
  if (plot) { hist(pval, breaks=100, xlim=c(0,.5),col="tan") }

  pcat = cut(pval, breaks=c(0,.01,.05,.10,1),
                   labels=c(".01",".05",".10","large"),include.lowest=T)
  tab = table(pcat)
  res = chisq.test(tab, p=tab0) # p-value from X-squared and df=3
  #return(cbind(tk,cf.exp,pval,pcat))
  return(list(alpha=alpha,minus.pll=minus.pll,X2=res$statistic,
              p.01=tab[".01"],p.05=tab[".05"],p.10=tab[".10"],pval=pval))
 }

 # main: heuristic search through partial trees
 if (search=="heuristic"){
  mystat = data.frame(number.edges=0:length(internal.edge))
  mystat$newedge=NA; mystat$alpha=NA; mystat$negPseudoLoglik=NA; mystat$X2=NA;
  mystat$p.01=NA; mystat$p.05=NA; mystat$p.10=NA;
  for (Nedge in mystat$number.edges){ 
   if (Nedge==0){ edge2keep = c(); # starting from panmixia 
   } else { 
    if (Nedge==1){ nextedge = as.numeric(names(which.max(res2[2,]))) }
    if (Nedge> 1){
     edge2keep.2 = which(colnames(res2) %in% as.character(edge2keep))
     nextedge = as.numeric(names(which.max(res2[2,-edge2keep.2])))
    }
    mystat$newedge[Nedge+1] = nextedge
    edge2keep = c(edge2keep, nextedge)
    cat("new edge to keep:",nextedge,"\n")
   }
   res1 = test.species.tree(edge2keep,plot=F)
   mystat$negPseudoLoglik[Nedge+1] = res1$minus.pll
   mystat$alpha[Nedge+1] = res1$alpha
   mystat$X2[Nedge+1]    = res1$X2
   mystat$p.01[Nedge+1]  = res1$p.01
   mystat$p.05[Nedge+1]  = res1$p.05
   mystat$p.10[Nedge+1]  = res1$p.10
   res2 = score.edge(which(res1$pval<.01));
   # print(res2[,order(res2[2,],decreasing=T)]) # to see edge scores
  }
  return(mystat)
 }
 # main: stepwise, forward+backward search through partial trees
 if (search=="both"){
  if (all(startT=="panmixia")){ edge2keep.current = c() }         else {
    if (all(startT=="fulltree")){edge2keep.current=internal.edge} else {
      if (is.numeric(startT)){ edge2keep.current = startT }       else {
        stop('problem with bad startT. Options: "panmixia", "fulltree", or numeric vector of edges')
      }
    }
  }
  Nedge = length(edge2keep.current)
  bestRes1 = test.species.tree(edge2keep.current,plot=F)
  res2 = score.edge(which(bestRes1$pval<.01))
  if (method=="PLL"){ bestCriterion = bestRes1$minus.pll } # criterion
  iter = 0
  {cat("iter=",iter,"\n\tNedge=",Nedge,"\n\tedges=",edge2keep.current)
  cat("\n\talpha=",bestRes1$alpha,"\n\t-PLL =",bestRes1$minus.pll)
  cat("\n\tX2   =",bestRes1$X2,"\n\tcrit.=",bestCriterion,"\n")}
  while (iter<maxiter){
   iter = iter+1
   bestAction = "stop"
   bestEdge = NULL
   # forward search
   if (Nedge<length(internal.edge)){
    candidates = setdiff(internal.edge,edge2keep.current)
    myk = min(kbest, length(candidates))
    candidates = res2[2,as.character(candidates)]
    edges2add = sort.int(candidates,index.return=T,decreasing=T)$ix[1:myk]
    edges2add = as.numeric(names(candidates[edges2add]))
    #cat("\tForward: trying edges",edges2add,"\n")
    for (ed in edges2add){
     edge2keep = c(edge2keep.current, ed)
     res1 = test.species.tree(edge2keep,plot=F)
     #cat("\t edge:",ed,"-PLL=",res1$minus.pll,"\n")
     if (method=="PLL"){
      if (res1$minus.pll < bestCriterion) {
       bestCriterion = res1$minus.pll
       bestAction="add"
       bestEdge = ed
       bestRes1 = res1
     }}
    }
   }
   # backward search
   if (Nedge >0){
    #myk = min(kbest,Nedge)
    #candidates = res2[2,as.character(edge2keep.current)]
    #edges2remove = sort.int(candidates,index.return=T,decreasing=T)$ix[1:myk]
    #edges2remove = as.numeric(names(candidates[edges2remove]))
    #edges2removeI = which(edge2keep.current %in% edges2remove)
    #cat("\tBackward: trying edges",edges2remove,"\n")
    #cat("\tBackward:\n")
    for (i in 1:length(edge2keep.current)){ # or in edges2removeI
     edge2keep = edge2keep.current[-i]
     res1 = test.species.tree(edge2keep,plot=F)
     #cat("\t edge:",edge2keep.current[i],"-PLL=",res1$minus.pll,"\n")
     if (method=="PLL"){
      if (res1$minus.pll < bestCriterion) {
       bestCriterion = res1$minus.pll
       bestAction="remove"
       bestEdge = i # index. the edge itself is edge2keep.current[i]
       bestRes1 = res1
     }}
    }
   }
   cat("iter=",iter,"\n\taction=",bestAction,", edge=",bestEdge,"\n")
   if (bestAction=="stop"){ break }
   if (bestAction=="add"){ edge2keep.current=c(edge2keep.current,bestEdge)}
   if (bestAction=="remove"){ edge2keep.current=edge2keep.current[-bestEdge]}
   res2 = score.edge(which(bestRes1$pval<.01))
   Nedge = length(edge2keep.current)
   #{cat("\tNedge=",Nedge,"\n\tedges=",edge2keep.current)
   #cat("\n\talpha=",bestRes1$alpha,"\n\t-PLL =",bestRes1$minus.pll)
   #cat("\n\tX2   =",bestRes1$X2,"\n\tcrit.=",bestCriterion,"\n")}
  }
  return(list(Nedge=Nedge,edges=edge2keep.current,
              notincluded = setdiff(internal.edge,edge2keep.current),
              alpha=bestRes1$alpha,
              negPseudoLoglik=bestRes1$minus.pll, X2=bestRes1$X2,
              p.01=bestRes1$p.01, p.05=bestRes1$p.05, p.10=bestRes1$p.10))
 }
}


#--------------------------------------------------------#
#      Same but standalone functions                     #
#      to analyze just one tree                          #
#--------------------------------------------------------#

test.tree.preparation <- function(cf, guidetree){
 # change all external edges of the guide tree to 0
 external.edge <- which(guidetree$edge[,2]<=length(guidetree$tip.label))
 internal.edge <- which(guidetree$edge[,2]> length(guidetree$tip.label))
 guidetree$edge.length[external.edge]=0
 # this is modifying guidetree in the current environment only, not in gloval env.
 Ntip <- length(guidetree$tip.label)
 tax2id <- 1:Ntip; names(tax2id)=guidetree$tip.label
 M = nrow(cf) # number of quartets with data
 dat.names <- as.matrix(cf[,1:4]) # -> taxon names as characters, not factors
 dat.leaves <- matrix(tax2id[dat.names],nrow=M,ncol=4)
 #cat("dat.names: \n"); print(head(dat.names))
 #cat("dat.leaves:\n"); print(head(dat.leaves))

 tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
 names(tab0)=c(".01",".05",".10","large") # tab0*M:  .01    .05     .10   large
                                          #       274.05 1096.2 1370.25 24664.5 
 cat("determining tree traversal node post-order... ")
 nodes.postorder = rep(NA,Ntip+guidetree$Nnode)
 guidetree$node.order = rep(NA,Ntip+guidetree$Nnode)
 nextnode <- 1
 set.node.postordertraversal = function(nodeid){
  # external: nextnode, nodes.postorder, guidetree$node.order
  if (nodeid > Ntip) { # not a leaf
   ces = which(guidetree$edge[,1]==nodeid) # child edges
   for (ce in ces){
    set.node.postordertraversal(guidetree$edge[ce,2])
   }
  }
  guidetree$node.order[nodeid] <<- nextnode
  nodes.postorder[nextnode] <<- nodeid
  nextnode <<- nextnode+1
 }
 set.node.postordertraversal(Ntip+1) # ID for the root
 cat("done.\n")
 #plot(guidetree);tiplabels(adj=0);nodelabels(adj=1);edgelabels()
 #cat("nodes.postorder:\n"); print(nodes.postorder)
 #cat("guidetree$node.order:\n"); print(guidetree$node.order)

 cat("calculating matrix of descendant relationships... ")
 node2descendant = matrix(FALSE, Ntip+guidetree$Nnode, Ntip+guidetree$Nnode)
 # entry [i,j] is TRUE if j is a descendant to (or equal to) i.
 for (n in nodes.postorder){
  node2descendant[n,n] <- TRUE
  if (n > Ntip) { # not a leaf
   ces <- which(guidetree$edge[,1]==n) # child edges     
   for (cn in guidetree$edge[ces,2]){ # child nodes
    node2descendant[n,] <- node2descendant[cn,] | node2descendant[n,]
   }
  }
 }
 cat("done.\n")

 get.mrca <- function(nodes){
   i1 <- max(guidetree$node.order[nodes]) # first candidate node (postorder index)
   for (n in nodes.postorder[i1:(Ntip+guidetree$Nnode)])
     if (all(node2descendant[n,nodes])){ return(n) }
 }
 cat("calculating matrix of edges spanned by each quartet... ")
 dominant.table  = c(1:3,1:3) # ie tree displays 12|34 if dominant = 1
                              #                  13|24               2
                              #                  14|23               3
 dominant = rep(NA,M)
 quartet2edge = matrix(FALSE,M,length(guidetree$edge.length)) # big matrix!

 for (i in 1:M){
  x <- dat.leaves[i,]
  ancestor <- sapply(list(x[1:2],x[c(1,3)],x[c(1,4)],x[3:4],x[c(2,4)],x[2:3]),get.mrca)
  ind = which.min(guidetree$node.order[ancestor])
  a1 = ancestor[ind]
  dominant[i] = dominant.table[ind]
  mrca = ancestor[which.max(guidetree$node.order[ancestor])]
  ancestor[ind] = Ntip+1 # replacing by the root to find 2nd best, next
  a2 = ancestor[which.min(guidetree$node.order[ancestor])]
  mynode = a1
  while (mynode!=a2 & mynode!=mrca){
   parentE = which(guidetree$edge[,2]==mynode)
   quartet2edge[i,parentE]=TRUE
   mynode = guidetree$edge[parentE,1]
  }
  if (mynode == mrca){ # the mrca of all 4 taxa is on path between both ancestors
   mynode = a2
   while (mynode!=mrca){
    parentE = which(guidetree$edge[,2]==mynode)
    quartet2edge[i,parentE]=TRUE
    mynode = guidetree$edge[parentE,1]
   }
  }
 } # i=4288; dat.names[i,]; quartet2edge[i,]
 colnames(quartet2edge)=as.character(1:length(guidetree$edge.length))
 cat("done.\n")

 total.ew = rep(1,M) %*% quartet2edge[,internal.edge] # number quartet per edge
 qw =  1/quartet2edge %*% rep(1,dim(quartet2edge)[2])
 # quartet weight: 1/Ne for each edge, where Ne= # edges for that quartet
 qw.ew = t(qw) %*% quartet2edge[,internal.edge]
 score.edge = function(quartetID){
  ues = rep(1,length(quartetID)) %*% quartet2edge[quartetID,internal.edge] / total.ew
  wes = qw[quartetID] %*% quartet2edge[quartetID,internal.edge] / qw.ew
  return(rbind(ues,wes))
 }
 return(list(quartet2edge = quartet2edge,
             dominant     = dominant))
}

plot.species.tree <- function(guidetree,edge.keep){
  sp.tree = guidetree
  if (length(edge.keep)==0){
   warning("Complete panmixia: the tree looks like a single tip. Won't plot this.")
   return(NULL)
  }
  sp.tree$edge.length[-edge.keep]=0
  plot(sp.tree)
}

test.one.species.tree <- function(cf,guidetree,prep,edge.keep,plot=FALSE){

  # prep should be the result of this, which takes a little while
  # so it's best to do it once and re-use it multiple times later:
  # prep = test.tree.preparation(cf,guidetree)
  
  M = nrow(cf) # number of quartets with data
  cf.obs = as.matrix(cf[,5:7])
  logCFsum = sum(log(cf.obs))
  tab0 = c(.01,.04,.05,.90)                # expected proportions of p-values
  names(tab0)=c(".01",".05",".10","large") # tab0*M:  .01    .05     .10   large
                                           #       274.05 1096.2 1370.25 24664.5
  if (plot) plot.species.tree(guidetree,edge.keep)
  
  #----------- get expected concordance factors -----------------------#
  if (length(edge.keep)==0) { tk = numeric(M) } else { 
    tk = prep$quartet2edge[,edge.keep] %*% cbind(guidetree$edge.length[edge.keep]) }
  cf.exp = matrix(exp(-tk)/3, M,3) # minor CFs only, so far
  cf.exp[1:M + (prep$dominant-1)*M] = 1-exp(-tk)*2/3
  
  #--------------------- estimate alpha -------------------------#
  if (plot){ plot(cf.exp,cf.obs,xlab="Observed CFs",ylab="Expected CFs"); abline(a=0,b=1,col="red")}

  logCFtilde = sum(log(cf.obs)*cf.exp)/M 
  tmp = rle(sort(tk))
  tk.unique = tmp$values
  nk = tmp$lengths # number of quartets with internal edge of length tk.
  wk = nk/M  # weight for quartets along edge k. They sum to 1.
  pk.minor = exp(-tk.unique)/3
  pk.domin = 1-2*pk.minor
  # f1 = -logPL/M with constant term logCFbar ommitted
  # f2 = d(f1)/d(alpha)
  f1 = function(alpha){-lgamma(alpha) - alpha * logCFtilde +
                         sum(wk*lgamma(alpha*pk.domin)) +
                       2*sum(wk*lgamma(alpha*pk.minor))   }
  f2 = function(alpha){-digamma(alpha) - logCFtilde +
                         sum(wk*pk.domin*digamma(alpha*pk.domin)) +
                       2*sum(wk*pk.minor*digamma(alpha*pk.minor)) }
  # now need to solve f2=0, or minimize f1
  # sa = (2/9)*(3*M)/sum((cf.obs-cf.exp)^2) # starting alpha value
  # optim(sa,f1,gr=f2, method="L-BFGS-B",lower=0,upper=+Inf)
  # alpha=res$par;  val=res$value
  res=optimize(f1,interval=c(0,10^5))
  alpha=res$minimum; val=res$objective
  minus.pll = val * M +logCFsum # - pseudo-log-likelihood
  # cat("estimated alpha:",alpha,", -pseudo log-lik:",minus.pll,"\n")
  # xx = seq(1,4*alpha,by=1); yy=sapply(xx,f1)
  # plot(xx,yy,type="l"); points(alpha,res$objective,col="red",pch=16)
  
  #--------------------- get a p-value for each quartet -------------------#
  p.exp = pmax(cf.exp[,1], cf.exp[,2],cf.exp[,3])
  ind.tree = which(tk>0)  
  ind.panmixia = which(tk==0)
  p.obs = pmax(cf.obs[,1], cf.obs[,2],cf.obs[,3]) # max for panmixia test
  # cf.obs[i+(j-1)*M] is same as cf.obs[i,j]
  p.obs[ind.tree] = cf.obs[ind.tree+(prep$dominant[ind.tree]-1)*M] 

  pval = rep(NA,M) # p-values
  pval[ind.panmixia] = 3*pbeta(p.obs[ind.panmixia],lower.tail=F,
                                  shape1=alpha/3,shape2=alpha*2/3)
  dev = abs(p.exp[ind.tree]-p.obs[ind.tree]) # deviations
  pe = p.exp[ind.tree] # expected proportions
  pval[ind.tree] =
    pbeta(pe+dev,shape1=alpha*pe,shape2=alpha*(1-pe),lower.tail=F) +
    pbeta(pe-dev,shape1=alpha*pe,shape2=alpha*(1-pe),lower.tail=T)
  pval[pval>1] = 1
  if (plot) { hist(pval, breaks=100, xlim=c(0,.5),col="tan") }

  pcat = cut(pval, breaks=c(0,.01,.05,.10,1),
                   labels=c(".01",".05",".10","large"),include.lowest=T)
  tab = table(pcat)
  res = chisq.test(tab, p=tab0) # p-value from X-squared and df=3
  #return(cbind(tk,cf.exp,pval,pcat))
  return(list(alpha=alpha,minus.pll=minus.pll,X2=res$statistic,
              p.01=tab[".01"],p.05=tab[".05"],p.10=tab[".10"],pval=pval,
              cf.exp=cf.exp))
}
