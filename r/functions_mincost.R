
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
      dplyr::mutate(file="A", supply=-iweight), # note the minus sign --- the A file receives weights
    bfile |> 
      dplyr::select(id, node, abrow=brow, weight, weightadj, iweight) |> 
      dplyr::mutate(file="B", supply=iweight)) |> # note the minus sign --- the A file receives weights
    dplyr::select(id, node, file, abrow, everything())
  
  return(nodes)
}


get_arcs <- function(dbtoa, datob, nodes){
  # arcs: all demands, nearest neighbor supply indexes (supply to demand -- stod)
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
  
  # arcs: all supplies, nearest neighbor demand indexes
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
    dplyr::mutate(dist=as.integer(dist*100.)) |> 
    dplyr::arrange(anode, bnode)
  
  return(arcs)
}


prepab <- function(afile, bfile, idvar, wtvar, xvars, k=NULL){
  a <- proc.time()
  # flows are from B to A
  # create a node file
  afile1 <- afile |> 
    dplyr::select(all_of(c(idvar, wtvar, xvars))) |> 
    dplyr::rename(id = !!as.symbol(idvar),
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


get_abfile <- function(arcs, nodes, flows, afile, bfile, idvar, wtvar, xvars){
  
  print("preparing base abfile with xvars from the afile")
  print("merge with afile for yvars and with bfile for zvars")
  print("merge on a_ and b_ followed by your id variable name and ...")
  abfile <- arcs |> 
    dplyr::mutate(weight=flows) |> 
    dplyr::filter(weight > 0) |>  # important
    dplyr::left_join(nodes |> 
                       dplyr::filter(file=="A") |> 
                       dplyr::select(aid=id, arow=abrow, a_weight=iweight),
                     by = join_by(arow)) |> 
    dplyr::left_join(nodes |> 
                       dplyr::filter(file=="B") |> 
                       dplyr::select(bid=id, brow=abrow, b_weight=iweight),
                     by = join_by(brow)) |> 
    # convert the a and b id variavbles to user-recognizable names
    dplyr::select(anode, bnode, aid, bid, a_weight, b_weight, dist, weight) |> 
    dplyr::rename(!!paste0("a_", idvar):=aid,
                  !!paste0("b_", idvar):=bid) |> 
    left_join(afile |> 
                select(-all_of(wtvar)) |> 
                rename(!!paste0("a_", idvar):=!!sym(idvar)) |> 
                rename(!!!setNames(xvars, paste0("a_", xvars))),
              by=join_by(!!paste0("a_", idvar))) |> 
    left_join(bfile |> 
                select(-all_of(wtvar)) |> 
                rename(!!paste0("b_", idvar):=!!sym(idvar)) |> 
                rename(!!!setNames(xvars, paste0("b_", xvars))),
              by=join_by(!!paste0("b_", idvar))) |> 
    dplyr::arrange(anode, dist)
  
  return(abfile)
}


matchab <- function(afile, bfile, idvar, wtvar, xvars, k=NULL){
  
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
  
  abfile <- get_abfile(res$prep_list$arcs, 
                       res$prep_list$nodes, 
                       mcfresult$flows, 
                       afile=afile, bfile=bfile, idvar=idvar, wtvar=wtvar, xvars=xvars)
  
  return(list(prep_list=prep_list, mcfresult=mcfresult, abfile=abfile))
}



