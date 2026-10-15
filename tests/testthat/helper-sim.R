# Simulated data sets shared by the tests. All simulations use the drcTMB parameterisation,
# f(x) = b1 + (b2 - b1) F(x), via the mean function of the package.

models <- c("logistic", "loglogistic", "weibull1", "weibull2", "lognormal",
            "michaelismenten", "asymptoticregression", "exponentialdecay", "gompertz")

# mean function for model name m and parameters b = c(b1, ..., b5)
drcmean <- function(m, x, b){
  drcTMB:::drcmeanR(match(m, models), x, b[1], b[2], b[3], b[4], b[5])
}

# true parameters and doses for a well-identified curve of each model
truepars <- list(
  logistic             = list(b=c(1, 8, 1.2, 5, 0),      dose=0:10),
  loglogistic          = list(b=c(1, 8, -1.5, 3, 0),     dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  weibull1             = list(b=c(1, 8, -1.2, 4, 0),     dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  weibull2             = list(b=c(1, 8, -1, 2, 0),       dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  lognormal            = list(b=c(1, 8, -1, 3, 0),       dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  michaelismenten      = list(b=c(1, 8, 0, 3, 0),        dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  asymptoticregression = list(b=c(1, 8, 0, 4, 0),        dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  exponentialdecay     = list(b=c(8, 0.5, 0, 4, 0),      dose=c(0, 0.5, 1, 2, 4, 8, 16, 32)),
  gompertz             = list(b=c(1, 8, 0.9, 4.5, 0),    dose=0:10)
)

# gaussian data with nrep replicates per dose
simgauss <- function(m, nrep=4, sd=0.4, seed=1){
  set.seed(seed)
  tp <- truepars[[m]]
  d <- data.frame(dose=rep(tp$dose, each=nrep))
  d$y <- drcmean(m, d$dose, tp$b) + rnorm(nrow(d), 0, sd)
  d
}

# grouped data: ng groups, every dose in every group, nrep replicates
simgroups <- function(dose=c(0, 0.5, 1, 2, 4, 8, 16), ng=10, nrep=2, seed=1){
  set.seed(seed)
  d <- expand.grid(rep=seq_len(nrep), dose=dose, grp=factor(seq_len(ng)))
  d$rep <- NULL
  d
}

# fit wrapper for the four-parameter versions of logistic and loglogistic (b5 fixed at 0)
fit4 <- function(..., start=c(b5=0), fix=c(b5=TRUE)){
  drmTMB(..., start=start, fix=fix)
}

# drmTMB coefficients without the variance parameters
bcoef <- function(fit){
  cf <- coef(fit)
  cf[grepl("^b[1-5]:", names(cf))]
}
