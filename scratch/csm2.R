
source(here::here("r", "libraries.r"))
# source(here::here("r", "functions.r"))
# source(here::here("r", "constants.r"))
library(optmatch)
# library(rlemon)


dir <- r"(E:\data\puf_files\puf2015)"

dfa <- arrow::read_parquet(fs::path(dir, "puf_2015.parquet"))

set.seed(123)
dfb <- dfa |> 
  sample_n(100) |> 
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

library(Rcpp)
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

cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame matchrecs(DataFrame adf, DataFrame bdf) {
  
  // Extracting columns
  IntegerVector ida = adf["ida"];
  NumericVector weighta = adf["weighta"];
  NumericVector ranka = adf["ranka"];
  
  IntegerVector idb = bdf["idb"];
  NumericVector weightb = bdf["weightb"];
  NumericVector rankb = bdf["rankb"];
  
  int ia = 0; // 0-based in C++
  int ib = 0;

  double dweighta = weighta[0]; // seed with info from the first record
  double dweightb = weightb[0];
  
  std::vector<int> records_ida, records_idb;
  std::vector<double> records_weighta, records_weightb, records_ranka, records_rankb;
  
  while (ia < ida.size() && ib < idb.size()) {
    if (dweightb > dweighta) { // done with the a record so write result
      records_ida.push_back(ida[ia]);
      records_idb.push_back(idb[ib]);
      records_weighta.push_back(dweighta);
      records_weightb.push_back(dweighta);
      records_ranka.push_back(ranka[ia]);
      records_rankb.push_back(rankb[ib]);
      
      dweightb = dweightb - dweighta;
      ia++;
      dweighta = weighta[ia]; // get the next weighta
      
    } else {
      records_ida.push_back(ida[ia]);
      records_idb.push_back(idb[ib]);
      records_weighta.push_back(dweightb);
      records_weightb.push_back(dweightb);
      records_ranka.push_back(ranka[ia]);
      records_rankb.push_back(rankb[ib]);
      
      dweighta = dweighta - dweightb;
      ib++;
      dweightb = weightb[ib]; // get the next weightb
    }
  }
  
  return DataFrame::create(_["ida"]=records_ida, _["idb"]=records_idb, _["weighta"]=records_weighta, 
                          _["weightb"]=records_weightb, _["ranka"]=records_ranka, _["rankb"]=records_rankb);
}
')

# df <- matchrecs(adf, bdf)
# sum(df$weighta)
# sum(df$weightb)

# now do it with real data -----

set.seed(123)
dfb <- dfa |> 
  sample_n(10000) |> 
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

set.seed(456)
afile2 <- afile |> 
  mutate(across(!!xvars, ~ .x * (1 + rnorm(n(), mean=0, sd=0.03))))
afile2

set.seed(789)
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

df <- matchrecs(adf, bdf)
sum(df$weighta)
sum(df$weightb)
glimpse(df)

# check - are all weight sums correct?
checka <- df |> 
  summarise(weighta=sum(weighta),
            weightb=sum(weightb),
            n=n(),
            .by=ida)
adf |> 
  left_join(checka |> select(ida, acalc=weighta, n), by = join_by(ida)) |> 
  mutate(diff=acalc - weighta) |> 
  arrange(desc(abs(diff)))

checkb <- df |> 
  summarise(weighta=sum(weighta),
            weightb=sum(weightb),
            n=n(),
            .by=idb)
bdf |> 
  left_join(checkb |> select(idb, bcalc=weightb, n), by = join_by(idb)) |> 
  mutate(diff=bcalc - weightb) |> 
  arrange(desc(abs(diff)))

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
