# library("qgraph")
library("IsingFit")
library("IsingSampler")

n.models <- 100
sampsize <- "100000"
estimates <- array(NA, dim=c(n.models,10))

for (i in 1:n.models) {
  
  df <- read.csv(paste("data_equilibria/basecase/equilibrium_nTypes5_k3_h1_pIntra0_25_c3_nObjects",sampsize,"_set",i-1,".csv",sep=""))[,2:6]
  
  if (sum(apply(df, 2, sum) != 0) == 6 & sum(apply(df, 2, sum) != n.models) == 6) {
    # regularized Ising network estimation
    m1 <- IsingFit(as.matrix(df), AND=TRUE, progressbar=FALSE, plot=FALSE)
    m1$weiadj <- m1$weiadj[1:5,1:5]
    estimates[i,] <- m1$weiadj[lower.tri(m1$weiadj)]
  }
  
}
est.value <- estimates[,2:11]
x11(); hist(est.value)
summary(as.vector(est.value))
select <- !is.na(apply(est.value, 1, sum))

truth <- read.csv("data_parameter_sets/basecase_subpopulation_success_parameter_sets_log_nTypes5_k3_h1_pIntra0_25_c3.csv")
true.value <- as.matrix(truth[,2:11]-truth[,12:21])
hist(true.value)
summary(as.vector(true.value))

x11(); plot(est.value ~ true.value,
            xlab="True parameter value", ylab="Estimated parameter value")
abline(a=0, b=1)


fspec <- function(d, i){
  d2 <- d[i,]
  spec <- sum(d2[,11:20][d2[,1:10] == 0] == 0)/sum(d2[,1:10] == 0)
  return(spec*100)
}

fsens <- function(d, i){
  d2 <- d[i,]
  sens <- sum(d2[,11:20][d2[,1:10] != 0] != 0)/sum(d2[,1:10] != 0)
  return(sens*100)
}

fppv <- function(d, i){
  d2 <- d[i,]
  ppv <- sum(d2[,11:20][d2[,1:10] != 0] != 0)/sum(d2[,11:20] != 0)
  return(ppv*100)
}

fscore <- function(d, i){
  d2 <- d[i,]
  sens <- sum(d2[,11:20][d2[,1:10] != 0] != 0)/sum(d2[,1:10] != 0)
  ppv <- sum(d2[,11:20][d2[,1:10] != 0] != 0)/sum(d2[,11:20] != 0)
  f.score <- 2 * (ppv * sens)/(ppv + sens)
  return(f.score)
}

fcorr <- function(d, i){
  d2 <- d[i,]
  pearsonr <- cor(d2[,1:10][d2[,11:20] != 0],d2[,11:20][d2[,11:20] != 0])
  return(pearsonr)
}

frmse <- function(d, i){
  d2 <- d[i,]
  residual <- d2[,11:20][d2[,11:20] != 0] - d2[,1:10][d2[,11:20] != 0]
  rmse <- sqrt(sum(residual^2)/length(residual))
  return(rmse)
}

testset <- cbind(true.value[select,], est.value[select,])

library(boot)

#turn off set.seed() if you want the results to vary
set.seed(626)

bootspec <- boot(testset, fspec, R=500)
bootspec
boot.ci(boot.out = bootspec, type = c("perc"))

bootsens <- boot(testset, fsens, R=500)
bootsens
boot.ci(boot.out = bootsens, type = c("perc"))

bootppv <- boot(testset, fppv, R=500)
bootppv
boot.ci(boot.out = bootppv, type = c("perc"))

bootfscore <- boot(testset, fscore, R=500)
bootfscore
boot.ci(boot.out = bootfscore, type = c("perc"))

bootcorr <- boot(testset, fcorr, R=500)
bootcorr
boot.ci(boot.out = bootcorr, type = c("perc"))

bootrmse <- boot(testset, frmse, R=500)
bootrmse
boot.ci(boot.out = bootrmse, type = c("perc"))
