#.............
# COI Data Adding and fitting
#.............
#

coi <- readRDS(paste0(gdrive,"/data/raw_data/coi_data.rds"))
data <- round(coi$coi)

dnbinom_fit <- function(par, data = data){

  sum(dnbinom(x = data, size = par[1], mu = par[2], log = TRUE))

}

fit <- optim(par = c("size"=1,"mu"=1),
             fn = dnbinom_fit,
             data = data,
             lower=c(0.001,0.001),upper=c(100000,25),
             method="L-BFGS-B",
             control = list(
               trace = TRUE,
               fnscale = -1,
               maxit = 10000
             ))
