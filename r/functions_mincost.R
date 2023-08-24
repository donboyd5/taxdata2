
# packages needed:
#   tidyverse: dplyr, tibble, tidyr
#   FNN
#   rlemon


get_distances <- function(afile, bfile, xvars, k=NULL){

  # determine how many near-neighbors to get, if not specified
  if(is.null(k)){
    maxrecs <- max(nrow(afile), nrow(bfile))
    minrecs <- min(nrow(afile), nrow(bfile))
    k = round(maxrecs) * .05
    k <- min(k, 1000) # k can never be more than 1000
    k <- max(k, 10) # k always must be at least 10
    k <- min(k, minrecs) # but k cannot be less than the number of rows in the shorter file 
  }
  print(paste0("k: ", k))
  
  # distance computations:
    # scale input data to mean=0, sd=1 before computing distances
    # compute nearest neighbors two ways to ensure that we have:
    #   arcs from each B record to k A records, 
    #   arcs from each A record to k B records
    # this does not guarantee feasibility but should help

  
  # k nearest distances for donating from file B to file A (dbtoa)
  # result matrices have same # rows as afile
  dbtoa <- FNN::get.knnx(bfile |> dplyr::select(!!xvars) |> scale(), 
                         afile |> dplyr::select(!!xvars) |> scale(),
                         k=k, algorithm="brute")    
  
  # k nearest distances for donating from file A to file B (datob)
  # result matrices have same # rows as bfile
  datob <- FNN::get.knnx(afile |> dplyr::select(all_of(xvars)) |> scale(), # scale to mean=0, sd=1 and compute distances
                    bfile |> dplyr::select(all_of(xvars)) |> scale(),
                    k=k, algorithm="brute") # brute is fastest algorithm based on testing

  
  return(list(dbtoa=dbtoa, datob=datob))
}


get_nodes <- function(afile, bfile){
  
  nodes <- dplyr::bind_rows(
    afile |> 
      dplyr::select(id, node, abrow=arow, weight, weightadj, iweight) |> 
      dplyr::mutate(file="A", supply=-iweight), # note the minus sign because the A file demands weights
    bfile |> 
      dplyr::select(id, node, abrow=brow, weight, weightadj, iweight) |> 
      dplyr::mutate(file="B", supply=iweight)) |> # note NO minus sign because the B file supplies weights
    dplyr::select(id, node, file, abrow, everything())
  
  return(nodes)
}


get_arcs <- function(dbtoa, datob, nodes){
  
  # create tibbles for btoa:
  #   idx has crossing of arow=1:nrow(a), neighbor=1:k indexes, brow=relevant row of b
  #   dist has crossing of arow=1:nrow(a), neighbor=1:k indexes, dist=relevant distance
  # create tibble for btoa with arow=1:nrow(a), neighbor=1:k, brow=relevant row or b
  
  dbtoa_idx <- tibble(arow=rep(1:nrow(dbtoa$nn.index), each=ncol(dbtoa$nn.index)),
                      neighbor=rep(1:ncol(dbtoa$nn.index), nrow(dbtoa$nn.index)),
                      brow=c(t(dbtoa$nn.index)))  
  
  dbtoa_dist <- tibble(arow=rep(1:nrow(dbtoa$nn.dist), each=ncol(dbtoa$nn.dist)),
                       neighbor=rep(1:ncol(dbtoa$nn.dist), nrow(dbtoa$nn.dist)),
                       dist=c(t(dbtoa$nn.dist)))

  dbtoa_arcs <- dplyr::left_join(dbtoa_idx, dbtoa_dist, by = join_by(arow, neighbor))
  # dbtoa_arcs
  
  datob_idx <- tibble(brow=rep(1:nrow(datob$nn.index), each=ncol(datob$nn.index)),
                      neighbor=rep(1:ncol(datob$nn.index), nrow(datob$nn.index)),
                      arow=c(t(datob$nn.index)))  
  
  datob_dist <- tibble(brow=rep(1:nrow(datob$nn.dist), each=ncol(datob$nn.dist)),
                       neighbor=rep(1:ncol(datob$nn.dist), nrow(datob$nn.dist)),
                       dist=c(t(datob$nn.dist)))
  
  datob_arcs <- dplyr::left_join(datob_idx, datob_dist, by = join_by(brow, neighbor))  
  
  # arcs: combine and keep unique arcs, keeping track of their neighbor status
  arcs1 <- dplyr::bind_rows(
    dbtoa_arcs |> 
      dplyr::mutate(src="btoa"), 
    datob_arcs |> 
      dplyr::mutate(src="atob")) |> 
    dplyr::select(arow, brow, neighbor, dist, src)
  
  # keeping the neighbor number can help in figuring out the quality of a match
  # the fastest way I could find to do this is 
  #   (1) get unique arow, brow arcs and their distances, which by definition are unique
  #   (2) merge back to get the atob and btoa neighbor numbers; 
  # note that I waste some memory by creating interim files
  
  arcs_distinct <- arcs1 |>
    dplyr::select(arow, brow, dist) |> 
    dplyr::distinct()
  
  arcs_neighbors <- arcs_distinct |>
    dplyr::left_join(arcs1 |> dplyr::filter(src=="btoa") |> 
                       dplyr::select(arow, brow, btoa_neighbor=neighbor),
              by = join_by(arow, brow)) |> 
    dplyr::left_join(arcs1 |> dplyr::filter(src=="atob") |> 
                       dplyr::select(arow, brow, atob_neighbor=neighbor),
              by = join_by(arow, brow)) |> 
    # prefer btoa_neighbor for later analysis of how far we had to go to find matches
    dplyr::mutate(neighbor=ifelse(is.na(btoa_neighbor), atob_neighbor, btoa_neighbor)) 
    
  # create the final arcs file: bring in node numbers
  arcs <- arcs_neighbors |> 
    dplyr::left_join(nodes |> dplyr::filter(file=="B") |> dplyr::select(bnode=node, brow=abrow),
                     by = join_by(brow)) |> 
    dplyr::left_join(nodes |> dplyr::filter(file=="A") |> dplyr::select(anode=node, arow=abrow), 
                     by = join_by(arow)) |> 
    dplyr::select(anode, bnode, arow, brow, dist, neighbor, btoa_neighbor, atob_neighbor) |> 
    
    # Convert distances, which are in standard deviation units because of scaling, 
    # to integers because the minimum cost flow algorithms require integer inputs.
    # Multiply by 100 to spread them out (otherwise we might have 0, 1, 2, 3 standard deviations)
    # I use 100 rather than a larger number, to keep the costs relatively small because
    # small costs to be important for minimum cost flow solvers.
    dplyr::mutate(dist=as.integer(dist*100.)) |> 
    dplyr::arrange(anode, bnode)
  
  return(arcs)
}


prepab <- function(afile, bfile, idvar, wtvar, xvars, k=NULL){
  
  a <- proc.time()
  
  # flows are from B to A
  # create a node file
  awtsum <- sum(afile[[wtvar]])
  bwtsum <- sum(bfile[[wtvar]])
  abratio <- awtsum / bwtsum
  print(paste0("initial ratio of sum of afile weights to bfile weights is: ", round(abratio, digits=3)))
  print("bfile weights will be adjusted as needed so that bfile weight sum equals afile weight sum")
  if(abratio < 0.75 || abratio > 1.25){
    print("however, large difference in sums suggests caution needed")
  }
  
  afile1 <- afile |> 
    dplyr::select(all_of(c(idvar, wtvar, xvars))) |> 
    dplyr::rename(id = !!as.symbol(idvar), # investigate a consistent way of converting strings to symbols
           weight = !!as.symbol(wtvar)) |>
    dplyr::mutate(file="A",
           arow=row_number(),
           node=row_number(),
           weightadj=weight,
           iweight=round(weightadj) |> as.integer())
  
  bfile1 <- bfile |> 
    dplyr::select(all_of(c(idvar, wtvar, xvars))) |> 
    dplyr::rename(id = !!as.symbol(idvar),
           weight = !!as.symbol(wtvar)) |>
    dplyr::mutate(file="B",
           brow=row_number(),
           node=row_number() + nrow(afile),
           weightadj=weight * sum(afile[[wtvar]]) / sum(weight),
           iweight=round(weightadj) |> as.integer())
  
  # print("balancing integer weights by adjusting bfile...")
  # this is rough - come up with a better way later
  awtsum <- sum(afile1$iweight)
  bwtsum <- sum(bfile1$iweight)
  diffba <- bwtsum - awtsum
  
  addval <- dplyr::case_when(diffba < 0 ~  1,
                      diffba > 0 ~ -1,
                      TRUE ~ 0)
  
  bfile1 <- bfile1 |> 
    dplyr::mutate(iweight=ifelse(row_number() <= abs(diffba),
                          iweight + addval,
                          iweight))
  
  # get distances
  dists <- get_distances(afile1, bfile1, xvars, k)
  
  nodes <- get_nodes(afile1, bfile1)
  arcs <- get_arcs(dbtoa=dists$dbtoa, datob=dists$datob, nodes=nodes)
  
  b <- proc.time()
  preptime <- (b - a)[3]
  
  return(list(nodes=nodes, arcs=arcs, preptime=preptime))
}


get_abfile <- function(arcs, nodes, flows, afile, bfile, idvar, wtvar, xvars, yvars, zvars){
  
  print("preparing base abfile...")
  abfile <- arcs |> 
    dplyr::mutate(weight=flows) |> 
    dplyr::filter(weight > 0) |>  # drop potential matches that weren't used
    
    # get the id and weight variables for the a and b files
    dplyr::left_join(nodes |> 
                       dplyr::filter(file=="A") |> 
                       dplyr::select(aid=id, arow=abrow, a_weight=iweight),
                     by = join_by(arow)) |> 
    dplyr::left_join(nodes |> 
                       dplyr::filter(file=="B") |> 
                       dplyr::select(bid=id, brow=abrow, b_weight=iweight),
                     by = join_by(brow)) |> 
    
    # convert the a and b id variable names to user-recognizable names
    dplyr::select(anode, bnode, aid, bid, neighbor, a_weight, b_weight, dist, weight) |> 
    dplyr::rename(!!paste0("a_", idvar):=aid,
                  !!paste0("b_", idvar):=bid) |> 
    
    # bring in each file's xvars, plus the yvars from a and zvars from b
    left_join(afile |> 
                dplyr::select(-all_of(wtvar)) |> 
                dplyr::rename(!!paste0("a_", idvar):=sym(idvar)) |> 
                dplyr::rename(!!!setNames(xvars, paste0("a_", xvars))), # give afile xvars an a prefix - do I need this many !!! ?
              by=join_by(!!paste0("a_", idvar))) |> 
    left_join(bfile |> 
                dplyr::select(-all_of(wtvar)) |> 
                dplyr::rename(!!paste0("b_", idvar):=sym(idvar)) |> 
                dplyr::rename(!!!setNames(xvars, paste0("b_", xvars))), # give bfile xvars a b prefix
              by=join_by(!!paste0("b_", idvar))) |> 
    
    # move variables around so that it is easier visually to compare the afile xvars to the bfile xvars
    dplyr::relocate(dist, .after=sym(paste0("b_", idvar))) |>  
    dplyr::relocate(all_of(yvars), .after = last_col()) |> 
    dplyr::relocate(all_of(zvars), .after = last_col()) |> 
    dplyr::arrange(anode, dist)
  
  return(abfile)
}


matchab <- function(afile, bfile, idvar, wtvar, xvars, yvars, zvars, k=NULL){
  
  print("preparing nodes and arcs...")
  prep_list <- prepab(afile,
                bfile,
                idvar=idvar, 
                wtvar=wtvar, 
                xvars=xvars,
                k=k)
  print(paste0("# seconds to prepare nodes and arcs: ", round(prep_list$preptime, 3)))
  
  a <- proc.time()
  # allowable_algorithms <- c("NetworkSimplex", "CostScaling", "CapacityScaling", "CycleCancelling")
  mcfresult <- rlemon::MinCostFlow(
    # flows are from B to A -- B has supply nodes, A has demand nodes
    arcSources=prep_list$arcs$bnode,
    arcTargets=prep_list$arcs$anode,
    arcCapacities=rep(max(abs(prep_list$nodes$supply)), nrow(prep_list$arcs)),
    arcCosts=prep_list$arcs$dist,
    nodeSupplies=prep_list$nodes$supply,
    numNodes=nrow(prep_list$nodes),
    algorithm = "NetworkSimplex" # NetworkSimplex seems fastest for these problems
  )
  b <- proc.time()
  mcfresult$mcftime <- (b - a)[3]
  print(paste0("# seconds to solve minimum cost flow problem: ", round(mcfresult$mcftime, 3)))
  print(paste0("Solution status: ", mcfresult$feasibility))
  
  abfile <- get_abfile(arcs=prep_list$arcs, 
                       nodes=prep_list$nodes, 
                       flows=mcfresult$flows, 
                       afile=afile, bfile=bfile, idvar=idvar, wtvar=wtvar,
                       xvars=xvars, yvars=yvars, zvars=zvars)
  
  return(list(prep_list=prep_list, mcfresult=mcfresult, abfile=abfile))
}



