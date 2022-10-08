

# routines for slide-seq simulations

require(Matrix)
require(Rcpp)
require(gtools)
sourceCpp("/home/hsarkar/Projects/SpatialContext/R_library/fast_sampling.cpp")


##' Simluate a slide-seq count matrix
##'
##' @param counts unnormalized count matrix (colums - cells, rows - genes)
##' @param cell.groups a cell-named factor specifying subpopulations of cells in the counts matrix
##' @param nbeads number of beads to simulate
##' @param cell.group.proportions proportions of different groups to simulate (defaults to the distribution of cell.groups)
##' @param bead.sizes an optional vector of bead sizes (defaults to sizes normally distributed around 1000 molecules with sd=100)
##' @param mixing.alpha a vector of dirichlet alpha parameters for cell mixing
##' @param ncores number of cores to use
##' @return a list including the following slots:
##'  $bead.counts - a simulated bead molecule count matrix (gene x bead)
##'  $bead.mixtures - a cell by bead matrix specifying mixing proportions for each bead
##'  $bead.entropy - entorpy vector across beads
slide.sim <- function(
  counts, 
  cell.groups=as.factor(setNames(colnames(counts),colnames(counts))), 
  nbeads=1e4, 
  cell.group.proportions=table(cell.groups)/length(na.omit(cell.groups)), 
  bead.sizes=rnorm(nbeads,1e3,1e2), 
  mixing.alpha=c(1,2), 
  ncores=10
  ) {

  # omit NA cells
  cell.groups <- as.factor(na.omit(cell.groups[names(cell.groups) %in% colnames(counts)]));
  counts <- counts[,names(cell.groups)]
  # the profile (count) matrix needs to be normalized so that the total sum of values for each cell is equal to 1
  pm <- t(t(counts)/colSums(counts));
  # check other arugments

  if(length(bead.sizes)!=nbeads) stop("specified bead.sizes vector must be equal to the nbeads")
  
  print(cell.group.proportions)
  if(!all(names(cell.group.proportions) == levels(cell.groups))) {
    if(any(! names(cell.group.proportions) %in% levels(cell.groups))) stop("cell.group.proportions specifies cell group names that are not part of cell.groups levels")
    if(length(intersect(names(cell.group.proportions),levels(cell.groups)))<1)  stop("cell.group.proportions need to include at least one level of cell groups")
    # fix cell.group.proportions to match the levels
    x <- setNames(rep(0,length(levels(cell.groups))), levels(cell.groups));
    x[names(cell.group.proportions)] <- cell.group.proportions;
    cell.group.proportions <- x;
  }

  # construct mixtures
  # two-level sampling: first we sample cell types that will participate in a mixture; then representative profiles for each cell type
  # a mapping from type to a set of profiles of that type
  pgl <- tapply(1:ncol(pm),cell.groups,function(ii) {
    if(length(ii)>1) { ii } else { c(ii,ii) } # work around for the unfortunate behavior of sample() function when x is of size 1
  },simplify=FALSE);
  # sampled profiles
  sampled.profiles <- unlist(lapply(1:nbeads,function(x) {
    sg <- sample.int(length(levels(cell.groups)), size=length(mixing.alpha), prob=cell.group.proportions, replace=FALSE)
    unlist(lapply(pgl[sg],sample,1))
  }))
  
  bead.mixtures <- sparseMatrix(
    i=sampled.profiles, # sample cells
    p=seq(0,nbeads*length(mixing.alpha),by=length(mixing.alpha)), # column offsets - easy, as every bead has a fixed number of cells mixing in - equal to length(mixing.alpha)
    x=as.vector(t( gtools::rdirichlet(nbeads, alpha=mixing.alpha))), # dirichlet-sampled mixing values
    dims=c(ncol(pm),nbeads)
  )
  rownames(bead.mixtures) <- colnames(counts); colnames(bead.mixtures) <- paste0('b',rep(1:ncol(bead.mixtures)));

  bead.entropy <- bead.mixtures <- drop0(bead.mixtures)
  bead.entropy@x <- bead.entropy@x * log(bead.entropy@x)
  bead.entropy <- -1*colSums(bead.entropy)

  # sample molecules
  bead.counts <- sample_mixtures(pm,bead.mixtures,bead.sizes,ncores=ncores)
  rownames(bead.counts) <- rownames(counts); colnames(bead.counts) <- colnames(bead.mixtures);

  return(list(bead.counts=bead.counts, bead.mixtures=bead.mixtures, bead.entropy=bead.entropy))
  
}

##' get pseudo-bulk profiles of clusters
##'
##' @param counts molecular count matrix (rows-genes, columns -cells)
##' @param cell.groups a named factor across cells
##' @return a sparse gene by level matrix
collapse.clusters <- function(counts,cell.groups) {
  cn <- intersect(names(cell.groups),colnames(counts))
  cell.groups <- as.factor(cell.groups[cn])
  if(length(cn)<2) stop("the matrix and the factor have too few common cells")
  x <- t(sccore::colSumByFactor(t(counts[,cn]),cell.groups)[-1,,drop=FALSE]) # remove the first NA column
  colnames(x) <- levels(cell.groups); rownames(x) <- rownames(counts);
  return(as(x,'dgCMatrix'))
}

##' Perturb an expression profile
##'
##' The procedure introduces an expression change into a profile by taking n.genes from upper quantile (defined by upper.q or q by default), and downregulating them to match genes sampled from the lower quantile (lower.q). And the other way around.
##' @param profile a count vector of genes
##' @param n.genes number of genes to change (divied evenly between upregulate and downregulate as long as they're true)
##' @param q quantile (to set upper.q and lower.q defaults) (default:0.2)
##' @param upper.q upper quantile, defaults to q
##' @param lower.q lower quantile, defualts to q
##' @param downregulate whether to decrease expression of select genes in the upper quantile
##' @param upregulate whether to increase expression of select genes in the lower quantile
##' @param swaps optional data frame with two columns ($target and $source gene profile integer ids), specifying which swaps should be made
##' @return a list with
##'   $profile - the modified count profile
##'   $altered.genes - a table specifying which genes were altered and by how much
##'   $swaps - gene swaps that were performed
simulate.expression.shift <- function(
  profile, n.genes=100, q=0.2, 
  upper.q=q, lower.q=q, 
  downregulate=TRUE, upregulate=TRUE,swaps=NULL) {
  if(is.null(swaps)) { # determine which genes should be swapped
    if(downregulate && upregulate) {
      n.genes.up <- n.genes.down <- round(n.genes/2)
    } else if(upregulate) {
      n.genes.up <- n.genes; n.genes.down <- 0;
    }  else { n.genes.down <- n.genes; n.genes.up <- 0; }
    
    eq <- quantile(profile,p=c(1-upper.q,lower.q))
    gns <- list(high=which(profile>=eq[1]), low=which(profile<=eq[2]))
    n.genes.up <- min(n.genes.up, length(gns$high), length(gns$low))
    n.genes.down <- min(n.genes.down, length(gns$high), length(gns$low))
    
    mp <- profile;
    swaps <- data.frame(source=c(),target=c());
    # sample genes for upregulation
    if(n.genes.up>0) {
      swaps <- rbind(swaps, 
      data.frame(target=sample(gns$low,n.genes.up), 
      source=sample(gns$high,n.genes.up)))
    }
    if(n.genes.down>0) {
      swaps <- rbind(swaps, 
      data.frame(target=sample(gns$high,n.genes.down), 
      source=sample(gns$low,n.genes.down)))
    }
  }
  
  # do the swaps
  mp[swaps$target] <- profile[swaps$source]

  # anontate altered genes
  cgi <- which(mp!=profile);
  df <- data.frame(cgi,delta=(mp-profile)[cgi]/sum(profile)*1e3,M=log2(mp[cgi]+1)-log2(profile[cgi]+1),row.names=names(profile)[cgi])
  
  return(list(profile=mp,altered.genes=df,swaps=swaps))
}
