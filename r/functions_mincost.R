
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
  dbtoa <- FNN::get.knnx(afile |> dplyr::select(all_of(xvars)) |> scale(), # scale to mean=0, sd=1 and compute distances
                    bfile |> dplyr::select(all_of(xvars)) |> scale(),
                    k=k, algorithm="brute") # brute is fastest algorithm based on testing
  
  # k nearest distances for donating from file A to file B (datob)
  datob <- FNN::get.knnx(bfile |> dplyr::select(!!xvars) |> scale(), 
                    afile |> dplyr::select(!!xvars) |> scale(),
                    k=k, algorithm="brute")
  
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
  # arcs: all demands, nearest neighbor supply indexes (supply to demand, i.e., b to a)
  # suppress messages because as_tibble is verbose with the neighbor names created by get_knnx
  suppressMessages({
  dbtoa_idx <- tibble::as_tibble(dbtoa$nn.index, .name_repair="unique") |> 
    dplyr::mutate(brow=row_number()) |> 
    tidyr::pivot_longer(-brow, names_to = "neighbor", values_to = "arow")
  })
  
  suppressMessages({
  dbtoa_dist <- tibble::as_tibble(dbtoa$nn.dist, .name_repair="unique") |> 
    dplyr::mutate(brow=row_number()) |> 
    tidyr::pivot_longer(-brow, names_to = "neighbor", values_to = "dist")
  })
  
  dbtoa_arcs <- dplyr::left_join(dbtoa_idx, dbtoa_dist, by = join_by(brow, neighbor)) |> 
    dplyr::mutate(neighbor=str_remove_all(neighbor, coll(".")) |> as.integer())
  # dbtoa_arcs
  
  # arcs: all supplies, nearest neighbor demand indexes (a to b)
  suppressMessages({
  datob_idx <- as_tibble(datob$nn.index, .name_repair="unique") |> 
    dplyr::mutate(arow=row_number()) |> 
    tidyr::pivot_longer(-arow, names_to = "neighbor", values_to = "brow")
  })
  
  suppressMessages({
  datob_dist <- as_tibble(datob$nn.dist, .name_repair="unique") |> 
    dplyr::mutate(arow=row_number()) |> 
    tidyr::pivot_longer(-arow, names_to = "neighbor", values_to = "dist")
  })
  
  datob_arcs <- left_join(datob_idx, datob_dist, by = join_by(arow, neighbor)) |> 
    dplyr::mutate(neighbor=str_remove_all(neighbor, coll(".")) |> as.integer())
  # datob_arcs
  
  # arcs: combine and keep unique arcs
  arcs1 <- bind_rows(dbtoa_arcs, datob_arcs)
  
  arcs1 <- arcs1 |> 
    dplyr::select(-neighbor) |> 
    dplyr::distinct()
  
  arcs <- arcs1 |> 
    dplyr::left_join(nodes |> dplyr::filter(file=="B") |> dplyr::select(bnode=node, brow=abrow),
                     by = join_by(brow)) |> 
    dplyr::left_join(nodes |> dplyr::filter(file=="A") |> dplyr::select(anode=node, arow=abrow), 
                     by = join_by(arow)) |> 
    dplyr::select(anode, bnode, arow, brow, dist) |> 
    
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
  print(paste0("ratio of sum of afile weights to bfile weights is: ", round(abratio, digits=3)))
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
  arcs <- get_arcs(dists$dbtoa, dists$datob, nodes)
  
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
    dplyr::select(anode, bnode, aid, bid, a_weight, b_weight, dist, weight) |> 
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



