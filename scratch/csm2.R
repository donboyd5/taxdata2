
source(here::here("r", "libraries.r"))
source(here::here("r", "functions_mincost.r"))
# source(here::here("r", "functions.r"))
# source(here::here("r", "constants.r"))
# library(Rcpp)
# library(optmatch)
# library(RANN)
# library(FNN) # get.knnx(data, query, k, algorithm)
# library(rlemon)

# sourceCpp(here::here("cpp", "functions_matching.cpp"))

# make pums testdata ----

library(tidycensus)
pvars2021 <- pums_variables |> 
  filter(year == 2021, survey == "acs1")

vars <- pvars2021 |> 
  filter(level=="person") |> 
  distinct(var_code, var_label, data_type)

testdata <- get_pums(
  variables = c("SEX", "MAR", "AGEP", "WKHP", "PINCP", "WAGP", "INTP", "SSP", "SSIP", "RETP", "SEMP"),
  state = "MA",
  survey = "acs1",
  year = 2021
)

count(testdata, SEX) # 1 male, 2 female
count(testdata, MAR) # 1=married

td2 <- testdata |> 
  lcnames() |> 
  filter(sex==1, mar==1, agep %in% 18:80) |> 
  mutate(pid=paste0(serialno, sporder)) |> 
  select(pid, weight=pwgtp, age=agep, hoursworked=wkhp,
         income=pincp, wages=wagp, interest=intp, socsec=ssp, pension=retp, selfemploy=semp)
quantile(td2$age)
anyDuplicated(td2$pid)

idvars <- c("pid", "weight")
xvars <- c("age", "hoursworked", "income")
yvars <- c("socsec", "selfemploy")
zvars <- c("interest", "pension", "wages")

set.seed(1234)
test_afile <- td2 |> 
  sample_n(10000) |> 
  select(all_of(c(idvars, xvars, yvars)))

test_bfile <- td2 |> 
  filter(!pid %in% test_afile$pid) |> 
  select(all_of(c(idvars, xvars, zvars))) |> 
  mutate(weight=weight * sum(test_afile$weight) / sum(weight))

sum(test_afile$weight)
sum(test_bfile$weight)

glimpse(test_afile)
glimpse(test_bfile)

saveRDS(test_afile, here::here("data", "test_afile.rds"))
saveRDS(test_bfile, here::here("data", "test_bfile.rds"))

# start file matching ----

na <- 15000
nb <- 5000

dfx <- dfa |> 
  filter(mars==2, e00100>=25e3, e00100<150e3) |> 
  sample_n(na + nb) |> 
  select(recid, weight, agi=e00100, wages=e00200, interest=e00300, dividends=e00600, businc=e00900, capgains=e01000, socsec=e02400) 
glimpse(dfx)

afile <- dfx |> 
  filter(row_number() <= na) 

bfile <- dfx |> 
  filter(row_number() %in% (na + 1):(na + nb))

a <- proc.time()
ldf <- prepab(afile,
               bfile,
               idvar="recid", 
               wtvar="weight", 
               xvars=c("agi", "capgains", "socsec"),
               k=100)
b <- proc.time()
b - a

res <- matchab(afile, bfile,
                   idvar="recid",
                   wtvar="weight", 
                   xvars=c("agi", "capgains", "socsec"), k=100)

str(res)
tmp <- res$abfile


rename(!!!setNames(paste0("a_", xvars), xvars))


result <- get_abfile(afile, bfile,
                   idvar="recid",
                   wtvar="weight", 
                   xvars=c("agi", "capgains", "socsec"), k=100)

str(result)

ab <- mac(res$prep_list$arcs, res$prep_list$nodes, res$mcfresult$flows, 
                 afile=afile, bfile=bfile, idvar="recid", wtvar="weight", xvars=c("agi", "capgains", "socsec"))


# match acs test files ----
library(tidyverse)
library(FNN) # get.knnx(data, query, k, algorithm)
library(rlemon)

source(here::here("r", "functions_mincost.r"))

# two random selections of records from American Community Survey in 2021
# Massachusetts, male, married, ages 18-80
# 10k records in afile, 4.3k records in bfile
afile <- readRDS(here::here("data", "test_afile.rds"))
bfile <- readRDS(here::here("data", "test_bfile.rds"))

idvars <- c("pid", "weight")
xvars <- c("age", "hoursworked", "income")
yvars <- c("socsec", "selfemploy") # different kinds of income
zvars <- c("interest", "pension", "wages") # different kinds of income

glimpse(afile) # idvars, xvars, yvars
glimpse(bfile) # idvars, xvars, zvars

# make sure neither file has duplicated ids
anyDuplicated(afile$pid)
anyDuplicated(bfile$pid)

# are sums of file weights reasonably close?
sum(afile$weight) / sum(bfile$weight)

# to help us choose k
nrow(afile); nrow(bfile) 
# k=100 seems to be reasonably good much of the time
# this allows each afile record to be matched to up to 100 of its nearest bfile neighbors
# and the same for each bfile record


res <- matchab(afile=afile, 
               bfile=bfile,
               idvar="pid",
               wtvar="weight", 
               xvars=xvars,
               yvars=yvars,
               zvars=zvars,
               k=100)

str(res)
res$mcfresult$cost
res$prep_list$arcs # this shows all of the allowable a-b matches; only some will have been selected

# check that the weights for all split people sum to the full person weight (as adjusted)
ab <- res$abfile 
ab |> 
  summarise(n=n(), a_weight=first(a_weight), weight=sum(weight), .by=a_pid) |> 
  mutate(diff=a_weight - weight) |> 
  arrange(desc(abs(diff)), desc(n))

ab |> 
  summarise(n=n(), b_weight=first(b_weight), weight=sum(weight), .by=b_pid) |> 
  mutate(diff=b_weight - weight) |> 
  arrange(desc(abs(diff)), desc(n))

# how far did we have to go to get matches, and what were their distances?
quantile(ab$neighbor, probs=c(0, .25, .5, .75, .9, .95, .99, 1))

ab |> 
  summarise(n=n(),
            dist=mean(dist),
            .by=neighbor) |> 
  arrange(neighbor) |> 
  mutate(pct=n / sum(n),
         cumpct=cumsum(pct)) |> 
  filter(cumpct <= .9)

# end file matching ----


# get data ----
# draw two similar subsets of the puf
dir <- r"(E:\data\puf_files\puf2015)"
dfa <- arrow::read_parquet(fs::path(dir, "puf_2015.parquet"))
glimpse(dfa)
count(dfa, mars)

# for later
idvars <- quote(c(node, suppdem, sdrow, recid, weight))
xvars <- quote(c(agi, capgains, socsec))
yvars <- quote(c(wages, businc))
zvars <- quote(c(interest, dividends))


nsupply <- 10
ndemand <- 20

df1 <- dfa |> 
  filter(mars==2, e00100>=25e3, e00100<50e3) |> 
  sample_n(nsupply + ndemand) |> 
  select(recid, weight, agi=e00100, wages=e00200, interest=e00300, dividends=e00600, businc=e00900, capgains=e01000, socsec=e02400) |> 
  mutate(node=row_number(),
         suppdem=ifelse(row_number() <= nsupply, "supply", "demand")) |>  # first recs are supplies, then demands
  mutate(sdrow=row_number(), .by=suppdem) |> 
  # adjust real demand weights to equal real supply weights
  mutate(adjweight=ifelse(suppdem=="demand",
                          weight*sum(weight[suppdem=="supply"]) / sum(weight[suppdem=="demand"]),
                          weight)) |> 
  # calculate unadjusted integer weights
  mutate(iweight1=round(adjweight) |> as.integer())

# is adjustment needed?
df1 |> 
  summarise(across(c(weight, adjweight, iweight1), sum), .by=suppdem)

# balance integer weights if needed
nodes <- df1 |> 
  mutate(supply=ifelse(suppdem=="supply", iweight1, -iweight1))
sum(nodes$supply)

# get rann nearest-neighbor costs -----

# get nearest neighbors in both directions

supplies <- nodes |> 
  filter(suppdem=="supply") |> 
  select(agi, capgains, socsec)

demands <- nodes |> 
  filter(suppdem=="demand") |> 
  select(agi, capgains, socsec)

# all demands, indexes for nearest supply neighbors
stod_nn <- nn2(supplies, demands, k=nsupply, eps=0.0) # rows are demands, cols are supplies; eps 0.0 is exact
stod_nn

# all supplies, indexes for nearest demand neighbors
dtos_nn <- nn2(demands, supplies, k=ndemand, eps=0.0) # rows are supplies, cols are demands; eps 0.0 is exact
dtos_nn

# get arcs both ways, then keep unique arcs ----

# arcs: all demands, nearest neighbor supply indexes (supply to demand -- stod)
stod_idx <- as_tibble(stod_nn$nn.idx, .name_repair="unique") |> 
  mutate(drow=row_number()) |> 
  pivot_longer(-drow, names_to = "neighbor", values_to = "srow")

stod_dist <- as_tibble(stod_nn$nn.dist, .name_repair="unique") |> 
  mutate(drow=row_number()) |> 
  pivot_longer(-drow, names_to = "neighbor", values_to = "dist")

stod_arcs <- left_join(stod_idx, stod_dist, by = join_by(drow, neighbor)) |> 
  mutate(neighbor=str_remove_all(neighbor, coll(".")) |> as.integer())
stod_arcs

# arcs: all supplies, nearest neighbor demand indexes
dtos_idx <- as_tibble(dtos_nn$nn.idx, .name_repair="unique") |> 
  mutate(srow=row_number()) |> 
  pivot_longer(-srow, names_to = "neighbor", values_to = "drow")

dtos_dist <- as_tibble(dtos_nn$nn.dist, .name_repair="unique") |> 
  mutate(srow=row_number()) |> 
  pivot_longer(-srow, names_to = "neighbor", values_to = "dist")

dtos_arcs <- left_join(dtos_idx, dtos_dist, by = join_by(srow, neighbor)) |> 
  mutate(neighbor=str_remove_all(neighbor, coll(".")) |> as.integer())
dtos_arcs

# arcs: combine and keep unique arcs
arcs1 <- bind_rows(stod_arcs, dtos_arcs)

arcs2 <- arcs1 |> 
  select(-neighbor) |> 
  distinct()

arcs <- arcs2 |> 
  left_join(nodes |> filter(suppdem=="demand") |> select(dnode=node, drow=sdrow), by = join_by(drow)) |> 
  left_join(nodes |> filter(suppdem=="supply") |> select(snode=node, srow=sdrow), by = join_by(srow)) |> 
  select(snode, dnode, srow, drow, dist) |> 
  mutate(dist=as.integer(dist)) |> 
  arrange(snode, dnode)
ht(arcs)


# lemon min cost flow ----

# https://errickson.net/rlemon/
# https://rdrr.io/cran/rlemon/man/
# https://rdrr.io/cran/rlemon/man/MinCostFlow.html
# https://github.com/somjit101/Min-Cost-Network-Flow-Lemon
# https://gist.github.com/Zhouxing-Su/6bb471f885dc29e5d25837abf78e7e5a
# http://lemon.cs.elte.hu/pub/doc/latest-svn/quicktour.html
# https://rdrr.io/cran/rlemon/src/R/mincostflow.R # "NetworkSimplex", "CostScaling", "CapacityScaling", "CycleCancelling"

# MinCostFlow(
#   arcSources,
#   arcTargets,
#   arcCapacities,
#   arcCosts,
#   nodeSupplies,
#   numNodes,
#   algorithm = "NetworkSimplex" # 
# )
# arcSources	Vector corresponding to the source nodes of a graph's edges
# arcTargets	Vector corresponding to the destination nodes of a graph's edges
# arcCapacities	Vector corresponding to the capacities of nodes of a graph's edges
# arcCosts	Vector corresponding to the capacities of nodes of a graph's edges
# nodeSupplies	Vector corresponding to the supplies of each node
# numNodes	The number of nodes in the graph
# algorithm	Choices of algorithm include "NetworkSimplex", "CostScaling", "CapacityScaling", and "CycleCancelling". NetworkSimplex is the default.


# be sure to convert all to integers
max(abs(supplies$isupply))
max(abs(nodes$supply))

res <- MinCostFlow(
  # flows are from B to A
  arcSources=arcs$snode,
  arcTargets=arcs$dnode,
  arcCapacities=rep(max(abs(nodes$supply)), nrow(arcs)),
  arcCosts=arcs$dist,
  nodeSupplies=nodes$supply,
  numNodes=nrow(nodes),
  algorithm = "NetworkSimplex"
)

res$feasibility
res$flows
res$potentials
res$cost

# put the results together and check
flows <- arcs |>
  mutate(flow=res$flows)
  
flows |> 
  summarise(flow=sum(flow), .by=snode)

nodes |> 
  filter(suppdem=="supply") |> 
  select(node, iweight1)

flows |> 
  summarise(flow=sum(flow), .by=dnode)

nodes |> 
  filter(suppdem=="demand") |> 
  select(node, iweight1)




small_graph_example # note that all of the following are integers
s1 <- small_graph_example$startnodes # supply labels
t1 <- small_graph_example$endnodes # termination labels
cap1 <- small_graph_example$arccapacity # capacities
costs1 <- small_graph_example$arccosts
n1 <- small_graph_example$nodedemand  # nodeSupplies, sums to zero
nnodes1 <- 34
# pick a problem
s <- s1; t <- t1; cap <- cap1; costs <- costs1; n <- n1; nnodes <- nnodes1 # lemon small example


tibble(s=s, t=t)

length(unique(s)) # 32
length(unique(t)) # 24
length(s); length(t); length(cap); length(costs); length(n)
# 274, 274, 274, 274, 34
sort(unique(c(s, t))) # 1:34
n
sum(n) # sums to zero
quantile(costs) # ranges from zero to 3,837
cbind(s, t) # first is from 2 to 14, cap appears to be 1, costs are 1194, s is 22, t is -22, flow is zero

CountBiEdgeConnectedComponents(s, t, nnodes)
CountConnectedComponents(s, t, nnodes)
CountStronglyConnectedComponents(s, t, nnodes) # 34 or 40
FindBiEdgeConnectedComponents(s, t, nnodes)
FindBiEdgeConnectedCutEdges(s, t, nnodes)
FindConnectedComponents(s, t, nnodes)
GetAndCheckTopologicalSort(s, t, nnodes)
GetTopologicalSort(s, t, nnodes)

IsAcyclic(s, t, nnodes)
IsBiEdgeConnected(s, t, nnodes)
IsBipartite(s, t, nnodes) # mine is bipartite, lemon small example is not
IsConnected(s, t, nnodes)
IsEulerian(s, t, nnodes) # lemon example false, my example true
IsParallelFree(s, t, nnodes)
IsSimpleGraph(s, t, nnodes)

s <- s1; t <- t1; cap <- cap1; costs <- costs1; n <- n1; nnodes <- nnodes1 # lemon small example
s <- s2; t <- t2; cap <- cap2; costs <- costs2; n <- n2; nnodes <- nnodes2 # my problem

IsStronglyConnected(s, t, nnodes)
MaxCardinalityMatching(s, t, nnodes) # lemon simple has 12, mine has 20
IsSimpleGraph(s, t, nnodes)
IsSimpleGraph(s, t, nnodes)
IsSimpleGraph(s, t, nnodes)
IsSimpleGraph(s, t, nnodes)
 
# https://lemon.cs.elte.hu/pub/doc/1.3.1/a00612.html
# In general, NetworkSimplex and CostScaling are the most efficient
# implementations. NetworkSimplex is usually the fastest on relatively small
# graphs (up to several thousands of nodes) and on dense graphs, while
# CostScaling is typically more efficient on large graphs (e.g. hundreds of
# thousands of nodes or above), especially if they are sparse. However, other
# algorithms could be faster in special cases. For example, if the total supply
# and/or capacities are rather small, CapacityScaling is usually the fastest
# algorithm (without effective scaling).

# These classes are intended to be used with integer-valued input data
# (capacities, supply values, and costs), except for CapacityScaling, which is
# capable of handling real-valued arc costs (other numerical data are required
# to be integer).


# let's look at arc 230 -- flow is 20
cbind(s, t)[230, ] # s 10, t 33 -- these are the nodes
costs[230] # zero cost
cap[230] # capacity 21
n[c(10, 33)] # supplies are 22, -198
out$potentials[c(10, 33)] # -1416 -1416
# so, arc 230, from node 10 to 33:
#   node 10 can supply 22 in total
#   node 33 wants 198 in total
#   arc has capacity 21
#   flow is 20 from node 10 to node 33


# NetworkSimplex CostScaling CapacityScaling CycleCancelling
out <- MinCostFlow(s, t, cap, costs, n, nnodes, algorithm = "CycleCancelling")
out$feasibility
out$flows # length is number of arcs
out$potentials # length is number of nodes
out$cost



# --- djb start here next ----  

library(dplyr)

# Sample data
set.seed(123)
nrecs <- 50
data <- tibble(supply = 1000+rnorm(nrecs)*100) |> 
  mutate(supply=ifelse(supply < median(supply), -supply, supply),
         supply=ifelse(supply < 0,
                        supply * sum(supply*(supply >= 0)) / sum(-supply*(supply < 0)),
                        supply))
data
sum(data$supply)

data <- data %>%
  mutate(
    isupply = floor(supply), # initial integer approximation
    frac = supply - isupply  # fractional part
  )

# If sum of isupply is too negative, adjust upwards
n_adj = -sum(data$isupply)
data = data %>%
  arrange(desc(frac)) %>%
  mutate(isupply = ifelse(row_number() <= n_adj, isupply + 1, isupply)) %>%
  select(-frac) 

sum(data$supply)
sum(data$isupply)
data
data |> 
  filter(isupply != round(supply))


# cur_data()` was deprecated in dplyr 1.1.0. ℹ Please use `pick()` instead
data <- data %>%
  mutate(isupply = floor(supply)) %>%
  arrange(supply - isupply) %>%
  mutate(
    isupply = ifelse(
      row_number() <= (-sum(cur_data()$isupply)), 
      isupply + 1, 
      isupply
    )
  )

data
sum(data$supply)
sum(data$isupply)


# test on the same data I used in julia ----
adf <- read_delim("
ida weighta ranka
4 79.9297 0.045742
1 93.4751 0.410937
8 74.8385 0.518504
7 49.0379 0.624424
10 5.8328 0.630887
5 66.3594 0.729227
9 56.304 0.747482
2 11.9394 0.858867
6 67.1042 0.979532
3 21.3759 0.986625",
delim=" ", col_names = TRUE, trim_ws = TRUE)

bdf <- read_delim("
idb weightb rankb
5 32.3814 0.0658773
3 110.878 0.171634
2 117.388 0.219503
4 170.165 0.580475
1 95.3846 0.759547",
delim=" ", col_names = TRUE, trim_ws = TRUE)

adf <- adf |> 
  arrange(desc(ranka))

bdf <- bdf |> 
  arrange(desc(rankb))

ab <- matchrecs(adf, bdf)
ab |> 
  arrange(ida, desc(rankb))

#    ida idb  weight    ranka     rankb
# 1    1   2 30.1455 0.410937 0.2195030
# 2    1   3 63.3296 0.410937 0.1716340
# 3    2   1  6.9045 0.858867 0.7595470
# 4    2   4  5.0349 0.858867 0.5804750
# 5    3   1 21.3759 0.986625 0.7595470
# 6    4   3 47.5484 0.045742 0.1716340
# 7    4   5 32.3813 0.045742 0.0658773
# 8    5   4 66.3594 0.729227 0.5804750
# 9    6   1 67.1042 0.979532 0.7595470
# 10   7   4 36.6339 0.624424 0.5804750
# 11   7   2 12.4040 0.624424 0.2195030
# 12   8   2 74.8385 0.518504 0.2195030
# 13   9   4 56.3040 0.747482 0.5804750
# 14  10   4  5.8328 0.630887 0.5804750

# julia result:
# Row │ ida    idb    weight    ranka     rankb
# 1 │     1      2  30.1459   0.410937  0.219503
# 2 │     1      3  63.3293   0.410937  0.171634
# 3 │     2      1   6.90455  0.858867  0.759547
# 4 │     2      4   5.03481  0.858867  0.580475
# 5 │     3      1  21.3759   0.986625  0.759547
# 6 │     4      3  47.5482   0.045742  0.171634
# 7 │     4      5  32.3814   0.045742  0.0658773
# 8 │     5      4  66.3594   0.729227  0.580475
# 9 │     6      1  67.1042   0.979532  0.759547
# 10 │     7      4  36.6342   0.624424  0.580475
# 11 │     7      2  12.4037   0.624424  0.219503
# 12 │     8      2  74.8385   0.518504  0.219503
# 13 │     9      4  56.304    0.747482  0.580475
# 14 │    10      4   5.8328   0.630887  0.580475


# test on real data ----

dir <- r"(E:\data\puf_files\puf2015)"

dfa <- arrow::read_parquet(fs::path(dir, "puf_2015.parquet"))

set.seed(12)
dfb <- dfa |> 
  sample_n(100000) |> 
  select(recid, weight, agi=e00100, wages=e00200, interest=e00300, dividends=e00600, businc=e00900, capgains=e01000, socsec=e02400)
dfb

idvars <- quote(c(recid, weight))
xvars <- quote(c(agi, capgains, socsec))
yvars <- quote(c(wages, businc))
zvars <- quote(c(interest, dividends))

afile <- dfb |> 
  select(!!idvars, !!xvars, !!yvars)
glimpse(afile)

bfile <- dfb |> 
  select(!!idvars, !!xvars, !!zvars)

set.seed(45)
afile2 <- afile |> 
  mutate(across(!!xvars, ~ .x * (1 + rnorm(n(), mean=0, sd=0.03))))
afile2

set.seed(78)
bfile2 <- bfile |> 
  mutate(across(!!xvars, ~ .x * (1 + rnorm(n(), mean=0, sd=0.05))))
bfile2

adf <- afile2 |> 
  select(ida=recid, weighta=weight, ranka=agi) |> 
  arrange(ranka)

bdf <- bfile2 |> 
  select(idb=recid, weightb=weight, rankb=agi) |> 
  mutate(weightb=weightb * sum(adf$weighta) / sum(weightb)) |> 
  arrange(rankb)

sum(adf$weighta)
sum(bdf$weightb)

# do the match ----
df <- matchrecs(adf, bdf)

# examine the match ----
sum(df$weight)
glimpse(df)

# check - are all weight sums correct?
checka <- df |> 
  summarise(weight=sum(weight),
            n=n(),
            .by=ida)
checka |> 
  count(n, order=TRUE)

adf |> 
  left_join(checka |> select(ida, acalc=weight, n), by = join_by(ida)) |> 
  mutate(diff=acalc - weighta) |> 
  arrange(desc(abs(diff))) |> 
  head()

checkb <- df |> 
  summarise(weight=sum(weight),
            n=n(),
            .by=idb)
bdf |> 
  left_join(checkb |> select(idb, bcalc=weight, n), by = join_by(idb)) |> 
  mutate(diff=bcalc - weightb) |> 
  arrange(desc(abs(diff))) |> 
  head()

abfile <- df |> 
  left_join(afile2 |> 
              select(ida=recid, !!xvars, !!yvars),
            by = join_by(ida)) |> 
  left_join(bfile2 |> 
              select(idb=recid, !!zvars),
            by = join_by(idb))

abfile |> 
  # this has the afile versions of the xvars, which won't be same as the bfile versions
  # but sums of yvars and zvars hit the targets!
  summarise(across(c(agi:dividends),
                   ~ sum(.x * weighta)))
afile2 |> 
  summarise(across(c(agi:businc),
                   ~ sum(.x * weight)))

bfile2 |> 
  summarise(across(c(agi:dividends),
                   ~ sum(.x * weight)))

df <- df |> 
  arrange(ida, idb)

tmp <- count(df, ida, sort=TRUE)
count(tmp, n)
recs <- tmp |> filter(n==4)
recs
id <- 152119 # 102800 152119 145909
check <- df |> 
  filter(ida==id)
check
sum(check$weighta)
sum(check$weightb)
adf |> 
  filter(ida==id)























# djb ----
afile <- dfb |> 
  select(!!idvars, !!xvars, !!yvars)
glimpse(afile)

bfile <- dfb |> 
  select(!!idvars, !!xvars, !!zvars)

set.seed(456)
afile2 <- afile |> 
  mutate(across(!!xvars, ~ .x * (1 + rnorm(n(), mean=0, sd=0.01))))
afile2

set.seed(789)
bfile2 <- bfile |> 
  mutate(across(!!xvars, ~ .x * (1 + rnorm(n(), mean=0, sd=0.01))))
bfile2

res <- mahalanobis.dist(data.x=afile2 |> select(!!xvars), data.y=bfile2 |> select(!!xvars))
res
dim(res)

res2 <- as_tibble(res) |> 
  mutate(obsa=row_number()) |> 
  pivot_longer(-obsa, names_to = "obsb") |> 
  arrange(obsa, value)

res3 <- res2 |> 
  group_by(obsa) |> 
  filter(row_number() <= 4) |> 
  ungroup()

stack <- bind_rows(afile2 |> mutate(ab=0),
                   bfile2 |> mutate(ab=1))


# optmatch ----------------------------------------------------------------
# https://cran.r-project.org/web/packages/optmatch/vignettes/fullmatch-vignette.html

library(optmatch)
data(nuclearplants)
head(nuclearplants)
help("nuclearplants")

nuclearplants$pt
table(nuclearplants$pt)
with(nuclearplants, table(pt))
nuke.nopt <- subset(nuclearplants, pt == 0)
ht(nuke.nopt)

# treated, control
table(nuke.nopt$pr)
# To get the pair match minimizing the mean paired distance on cap, among all collections of 7 non-overlapping pairs, do
tmp <- pairmatch(pr ~ cap, data = nuke.nopt) # a factor
print(pairmatch(pr ~ cap, data = nuke.nopt), grouped = TRUE)

abmatch <- pairmatch(ab ~ agi, data = stack)
abmatch <- pairmatch(ab ~ agi + capgains + socsec, data = stack)
print(abmatch, grouped=TRUE)
str(abmatch)

abmatch <- pairmatch(ab ~ agi + capgains + socsec, controls=4, data = stack)


#> 
#> Attaching package: 'Rcpp'
#> The following object is masked from 'package:inline':
#> 
#>     registerPlugin
cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')
# add works like a regular R function
add
#> function (x, y, z) 
#> .Call(<pointer: 0x7f96ecb3ef20>, x, y, z)
add(1, 2, 3)



library(Rcpp)
set.seed(123)

# Example data for adf
adf <- tibble(
  ida = 1:10,
  weighta = runif(10, 1, 5),
  ranka = rnorm(10)
)

# Example data for bdf
bdf <- tibble(
  idb = 1:8,
  weightb = runif(8, 2, 6),
  rankb = rnorm(8)
)
bdf <- bdf |> 
  mutate(weightb=weightb * sum(adf$weighta) / sum(weightb))
sum(adf$weighta)
sum(bdf$weightb)

adf <- adf |> 
  arrange(ranka)

bdf <- bdf |> 
  arrange(rankb)

sum(adf$weighta)
sum(bdf$weightb)


# df <- matchrecs(adf, bdf)
# sum(df$weighta)
# sum(df$weightb)

# now do it with real data -----


# test knn speed ----

# CONCLUSIONS:
#   d1d <- get.knnx(data, query, k=k, algorithm="brute") from FNN is fastest and best
#   n2 (fileb) small has a big impact
#   m # xvars small also big impact
#   k small also has a big impact
#   20k n, 5k n2, 10 m, k 10% solves in < 2 seconds -- reasonable real-world parameters
#   20k, 5k is 2x faster than 5k, 20k and MUCH faster if k is % of 20k the RHS file

# data<- query<- cbind(1:n, 1:n)
# get.knnx(data, query, k=10, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))

n <- 20000 # 1a 107 seconds, 1b 214, 81 1c , 74 1d, 105 2a , 108 2b
n2 <- 5000
m <- 10 # number of x vars used in distance calc - no material impact on speed
k <- round(n2 * .1) # for 20k, 20: time drops from 107 (10%) to 33 (5%) to  (1%)

set.seed(123)
data <- matrix(rnorm(n*m), ncol=m)
query <- matrix(rnorm(n2*m), ncol=m)

system.time(d1a <- get.knnx(data, query, k=k, algorithm="kd_tree")) # much faster than cover_tree when good; seems to hang sometimes
system.time(d1b <- get.knnx(data, query, k=k, algorithm="cover_tree"))
system.time(d1c <- get.knnx(data, query, k=k, algorithm="CR")) # DON'T USE - GIVES DIFFERENT INDEXES
system.time(d1d <- get.knnx(data, query, k=k, algorithm="brute")) # fast, maybe BEST
k2 <- round(n * .1) 
k2 <- 500
system.time(d1e <- get.knnx(query, data, k=k2, algorithm="brute")) # query, data about 2x as long as data, query, MUCH longer if we change k!!!
system.time(d2a <- nn2(data, query, k=k, eps=0.0)) # fast
system.time(d2b <- nn2(data, query, k=k, eps=1e-5))

d1a$nn.index[1:5, 1:5]
d1b$nn.index[1:5, 1:5]
d1c$nn.index[1:5, 1:5]
d1d$nn.index[1:5, 1:5]

d2a$nn.idx[1:5, 1:5]
d2b$nn.idx[1:5, 1:5]

d1a$nn.dist[1:5, 1:5]
d1b$nn.dist[1:5, 1:5]
d1c$nn.dist[1:5, 1:5]
d1d$nn.dist[1:5, 1:5]

tmp <- d1e$nn.dist
tmp[1, ]
plot(tmp[100, ])

d2a$nn.dist[1:5, 1:5]
d2b$nn.dist[1:5, 1:5]

# don't bother with fastknn - not maintained, based on ANN like other packages
library(devtools)
install_github("davpinto/fastknn")
library(mlbench)
library(caTools)
library(fastknn)

data(Ionosphere)

x <- data.matrix(subset(Ionosphere, select = -Class))
y <- Ionosphere$Class

set.seed(2048)
tr.idx <- which(sample.split(Y = y, SplitRatio = 0.7))
x.tr <- x[tr.idx,]
x.te <- x[-tr.idx,]
y.tr <- y[tr.idx]
y.te <- y[-tr.idx]

knn.out <- fastknn(xtr = x.tr, ytr = y.tr, xte = x.te, k = 10)

knn.out$class
knn.out$prob


# cppFunction('
# #include <Rcpp.h>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# DataFrame matchrecs(DataFrame adf, DataFrame bdf) {
#   
#   // Extracting columns
#   IntegerVector ida = adf["ida"];
#   NumericVector weighta = adf["weighta"];
#   NumericVector ranka = adf["ranka"];
#   
#   IntegerVector idb = bdf["idb"];
#   NumericVector weightb = bdf["weightb"];
#   NumericVector rankb = bdf["rankb"];
#   
#   int ia = 0; // 0-based in C++
#   int ib = 0;
# 
#   double dweighta = weighta[0]; // seed with info from the first record
#   double dweightb = weightb[0];
#   
#   std::vector<int> records_ida, records_idb;
#   std::vector<double> records_weighta, records_weightb, records_ranka, records_rankb;
#   
#   while (ia < ida.size() && ib < idb.size()) {
#     if (dweightb > dweighta) { // done with the a record so write result
#       records_ida.push_back(ida[ia]);
#       records_idb.push_back(idb[ib]);
#       records_weighta.push_back(dweighta);
#       records_weightb.push_back(dweighta);
#       records_ranka.push_back(ranka[ia]);
#       records_rankb.push_back(rankb[ib]);
#       
#       dweightb = dweightb - dweighta;
#       ia++;
#       dweighta = weighta[ia]; // get the next weighta
#       
#     } else {
#       records_ida.push_back(ida[ia]);
#       records_idb.push_back(idb[ib]);
#       records_weighta.push_back(dweightb);
#       records_weightb.push_back(dweightb);
#       records_ranka.push_back(ranka[ia]);
#       records_rankb.push_back(rankb[ib]);
#       
#       dweighta = dweighta - dweightb;
#       ib++;
#       dweightb = weightb[ib]; // get the next weightb
#     }
#   }
#   
#   return DataFrame::create(_["ida"]=records_ida, _["idb"]=records_idb, _["weighta"]=records_weighta, 
#                           _["weightb"]=records_weightb, _["ranka"]=records_ranka, _["rankb"]=records_rankb);
# }
# ')


# TPC approach -- my numbering
# https://www.taxpolicycenter.org/sites/default/files/alfresco/publication-pdfs/411136-The-Urban-Brookings-Tax-Policy-Center-Microsimulation-Model.PDF
# We implement predictive mean matching by using taxable income as the dependent variable and thus:

# 1) run the following regression separately for each cell using the PUF data:
#   Taxable Income = β0 + β1*(Dummy for the Aged Status) + β2*(Wages and Salaries) +
#   β3*(Taxable Interest) + β4*(Dividend Income) + β5*(Business income or loss) + β6*(Farm
#                                                                                     income or loss) + β7*(Schedule E Income) + β8*(Pensions) + β9*(Social Security Income) +
#   β10*(Unemployment Compensation) + β12*(Alimony) + β13*(Wage Share of Total Income) +
#   β14*(Capital Income Share of Total Income) + β15*(Dummy for Presence of Wage or Salary
#                                                     Income) 

# 2) calculate the fitted values of the Y and/or Z variables for
# both the host and donor files. That is, the coefficients on the X variables that were obtained in the
# regression using the host file, are then used to calculate fitted values for the Y and/or Z variables
# in both the host and donor files. Specifically, in the case of the tax model, we use the coefficients
# from the regression described above, along with the actual values of the explanatory variables in
# each file, to construct fitted values for taxable income for each record in both the CPS and the
# PUF. 

# 3) , the weights on each CPS record are multiplied by a factor such that the total of the CPS
# weights in each partitioned cell adds up to the total SOI weight for that partition

# 4) the records in each cell would then be sorted in
# descending order by the predicted values of one of the Y and/or Z variables. In the case of our tax
# model match, the cells are sorted by the predicted values of taxable income. Corresponding
# records from the PUF and the CPS are then matched within each partitioned cell. Of the two
# records, the one with the higher weight must be split or duplicated and matched with the next
# record or next several records in the other file until all of its weight has been “used up.” Thus,
# each record in the host PUF file is matched to that record in the donor CPS file that is “closest”
# in terms of having the most similar predicted value of taxable income among all records within
# the partition. Since the weights on the CPS file have been adjusted to equal the total PUF
# weights, all records are used in the match. One possible disadvantage of using all the records to
# perform the match is that some records might be matched despite having a large difference
# between the predicted values of taxable income in each of the files. 

# https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html

library(MatchIt)

data(lalonde)

vignette("matching-methods")
vignette("assessing-balance")
vignette("estimating-effects")
vignette("sampling-weights")

# Generalized Full Matching (method = "quick")
#  ?method_quick

# Generalize full PS matching
m.out1 <- matchit(treat ~ age + educ + race + nodegree +
                    married + re74 + re75, data = lalonde,
                  method = "quick")
m.out1
summary(m.out1)
str(m.out1)
str(m.out1$subclass)
levels(m.out1$subclass)
labels(m.out1$subclass)

tibble(labs=labels(m.out1$subclass)) |> 
  count(labs, sort=TRUE)


library(quickmatch)
# https://cran.r-project.org/web/packages/quickmatch/index.html



lalonde_matchit_nn <-
  matchit(
    treat ~ age + educ + black + hispan + nodegree + married + re74 + re75,
    baseline.group = 1,
    data = lalonde,
    method = "nearest",
    distance = "mahalanobis",
    subclass = T
  )


# Load package
library("quickmatch")

# Construct example data
my_data <- data.frame(y = rnorm(100),
                      x1 = runif(100),
                      x2 = runif(100),
                      treatment = factor(sample(rep(c("T", "C"), c(25, 75)))))

# Make distances
my_distances <- distances(my_data, dist_variables = c("x1", "x2"))

### Average treatment effect (ATE)

# Make matching
my_matching_ate <- quickmatch(my_distances, my_data$treatment)

# Covariate balance
covariate_balance(my_data$treatment, my_data[c("x1", "x2")], my_matching_ate)

# Estimate effect
lm_match(my_data$y, my_data$treatment, my_matching_ate)


### Average treatment effect of the treated (ATT)

# Make matching
my_matching_att <- quickmatch(my_distances, my_data$treatment, target = "T")

# Covariate balance
covariate_balance(my_data$treatment, my_data[c("x1", "x2")], my_matching_att, target = "T")

# Estimate effect
lm_match(my_data$y, my_data$treatment, my_matching_att, target = "T")


# https://stackoverflow.com/questions/63242715/nearest-neighbor-matching-with-the-mahalanobis-distance-in-r

# I would like to use the MatchIt package in R to perform nearest neighbor
# matching using the Mahalanobis distance withing some caliper. Which of the
# following two parameters of the matchit function that are related to the
# Mahalanobis distance should I use:
   
#   the distance="mahalanobis" param, or
#   the mahvars param (e.g., mahvars = c("X1", "X2")?
# What's the difference between the two?

# The documentation is terse about this (see pages 16 and 19): https://imai.fas.harvard.edu/research/files/matchit.pdf.

# You should use the latter. You need the distance argument to identify the
# propensity score that will be used to form the caliper. Setting mahvars will
# perform Mahalanobis distance matching on the mahvars variables, and the
# propensity score will be estimated based on the variables in the main formula.
# The caliper argument can then be specified, which defines the width of the
# caliper in units of the standard deviation of the propensity score.

# From https://cran.r-project.org/web/packages/MatchIt/vignettes/matching-methods.html:
# https://kosukeimai.github.io/MatchIt/reference/match.data.html

# Setting the distance="mahalanobis" and method="nearest" make MatchIt to run
# Nearest-neighbour matching with the Mahalanobis distance, without the
# consideration of the propensity score. And the covariates supplied in the main
# formula are used.

library(MatchIt)
data(lalonde)
glimpse(lalonde)

m.out1 <- matchit(treat ~ age + educ + married +
                    race + nodegree + re74 + re75,
                  data = lalonde, replace = TRUE,
                  caliper = .05, ratio = 4)

m.data1 <- match.data(m.out1, data = lalonde,
                      distance = "prop.score")

glimpse(m.data1)

g.matches1 <- get_matches(m.out1, data = lalonde,
                          distance = "prop.score")

dim(g.matches1) #multiple rows per matched unit

glimpse(g.matches1)

# https://github.com/kosukeimai/MatchIt/blob/master/src/nn_matchC.cpp

