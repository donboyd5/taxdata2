library(StatMatch)
library(clue)

# reproduce the statistical matching framework
# starting from the iris data.frame
suppressWarnings(RNGversion("3.5.0"))
set.seed(98765)
pos <- sample(1:150, 50, replace=FALSE)
ir.A <- iris[pos,c(1,3:5)]
ir.B <- iris[-pos, 2:5]
xx <- intersect(colnames(ir.A), colnames(ir.B))
xx # common variables "Petal.Length" "Petal.Width"  "Species"  

intersect(colnames(ir.A), colnames(ir.B)) # "Petal.Length" "Petal.Width"  "Species"  
union(colnames(ir.A), colnames(ir.B)) # "Sepal.Length" "Petal.Length" "Petal.Width"  "Species"      "Sepal.Width" 
dplyr::symdiff(colnames(ir.A), colnames(ir.B)) # "Sepal.Length" "Sepal.Width" 


# mixed.mtc
# ML estimation method under CIA ((rho_YZ|X=0));
# only parameter estimates (micro=FALSE)
# only continuous matching variables
xx.mtc <- c("Petal.Length", "Petal.Width") # don't use species
mtc.1 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width")
# estimated correlation matrix
mtc.1$cor
# ML estimation method under CIA ((rho_YZ|X=0));
# only parameter estimates (micro=FALSE)
# categorical variable 'Species' used as matching variable
xx.mtc <- xx
mtc.2 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width")
# estimated correlation matrix
mtc.2$cor
# ML estimation method with partial correlation coefficient
# set equal to 0.5 (rho_YZ|X=0.5)
# only parameter estimates (micro=FALSE)
mtc.3 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width",
                   rho.yz=0.5)
# estimated correlation matrix
mtc.3$cor
# ML estimation method with partial correlation coefficient
# set equal to 0.5 (rho_YZ|X=0.5)
# with imputation step (micro=TRUE)
mtc.4 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width",
                   rho.yz=0.5, micro=TRUE, constr.alg="Hungarian")
# first rows of data.rec filled in with z
head(mtc.4$filled.rec)
#
# Moriarity and Scheuren estimation method under CIA;
# only with parameter estimates (micro=FALSE)
mtc.5 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width",
                   method="MS")
# the starting value of rho.yz and the value used
# in computations
NND.hotdeck
mtc.5$rho.yz
# estimated correlation matrix
mtc.5$cor
# Moriarity and Scheuren estimation method
# with correlation coefficient set equal to -0.15 (rho_YZ=-0.15)
# with imputation step (micro=TRUE)
mtc.6 <- mixed.mtc(data.rec=ir.A, data.don=ir.B, match.vars=xx.mtc,
                   y.rec="Sepal.Length", z.don="Sepal.Width",
                   method="MS", rho.yz=-0.15,
                   micro=TRUE, constr.alg="lpSolve")
# the starting value of rho.yz and the value used
# in computations
mtc.6$rho.yz
# estimated correlation matrix
mtc.6$cor
# first rows of data.rec filled in with z imputed values
head(mtc.6$filled.rec)


str(mtc.6)
mtc.6$rho.yz
mtc.6$mu
mtc.6$vc
mtc.6$cor
mtc.6$phi
mtc.6$res.var
mtc.6$filled.rec
mtc.6$mtc.ids
mtc.6$dist.rd





