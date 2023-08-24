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
