
source(here::here("r", "libraries.r"))
# source(here::here("r", "functions.r"))
# source(here::here("r", "constants.r"))
library(Rcpp)
library(optmatch)
# library(rlemon)
sourceCpp(here::here("cpp", "functions_matching.cpp"))

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
  left_join(checka |> select(ida, acalc=weighta, n), by = join_by(ida)) |> 
  mutate(diff=acalc - weighta) |> 
  arrange(desc(abs(diff))) |> 
  head()

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

