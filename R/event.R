# Time-to-event (germination) data: family="event" of drmTMB, with the formula count ~ start + end. Each
# row is a time interval of a series (e.g. a Petri dish) with the number of events in it; a row with
# end = Inf holds the units without an event by the end of the experiment. The curve f(t) = b2 F(t),
# with b1 fixed at 0, is the cumulative proportion of events by time t, b2 the proportion of units with
# an event. The log-likelihood is the multinomial sum(count log(f(end) - f(start))), with f(Inf) = 1,
# without the multinomial coefficients, as drc's type="event".

# models whose F is a distribution function of the time, from 0 at time 0 to 1
eventmodels <- c("loglogistic", "weibull1", "weibull2", "lognormal", "gamma", "asymptoticregression", "michaelismenten")

# series of the rows of event data: the formula series (e.g. ~ dish), or by default the rows with the
# same covariates of the curve (the design matrices) and the same levels of the random effect terms
eventseries <- function(series, data, Xs, re){
  if (!is.null(series)){
    if (!inherits(series, "formula") || length(series) != 2) stop("series must be a one-sided formula, e.g. ~ dish.")
    return(droplevels(groupfactor(series[[2]], data)))
  }
  key <- do.call(paste, c(lapply(Xs, function(X) apply(X, 1, paste, collapse=",")),
                          lapply(re$terms, function(tm) as.character(groupfactor(str2lang(tm$group), data))), sep="|"))
  factor(key, levels=unique(key))
}

# template data of event data: end times (Inf as 0, with the indicator xinf), rows by series and start
# (indices from 0) and the series lengths, for the multinomial simulation
eventdata <- function(xs, xe, series){
  if (is.null(series)) return(list(xe=numeric(0), xinf=integer(0), eord=integer(0), elen=integer(0)))
  o <- order(series, xs)
  list(xe=ifelse(is.finite(xe), xe, 0), xinf=as.integer(!is.finite(xe)), eord=as.integer(o - 1),
       elen=as.integer(table(series)))
}

# whether every series is complete: intervals from 0 to Inf without gaps or overlaps, so that the counts
# of a series are multinomial with the probabilities of its intervals (simulate, residuals)
eventcomplete <- function(xs, xe, series){
  all(vapply(split(seq_along(xs), series), function(i){
    i <- i[order(xs[i])]
    xs[i[1]] == 0 && !is.finite(xe[i[length(i)]]) && all(abs(xe[i[-length(i)]] - xs[i[-1]]) < 1e-8)
  }, logical(1)))
}

# data for the starting values of event data: the cumulative counts of each series by the end of each
# interval, as binomial counts of the series totals, at the end times; rows with end = Inf get weight 0
eventstartdd <- function(dd, y=dd$y){
  s <- dd$series
  o <- order(s, dd$x)
  cum <- numeric(length(y))
  cum[o] <- ave(y[o], s[o], FUN=cumsum)
  tot <- ave(y, s, FUN=sum)
  fin <- is.finite(dd$xe)
  dd$x <- ifelse(fin, dd$xe, dd$x)
  dd$y <- cum
  dd$bn <- tot
  dd$w <- ifelse(fin & tot > 0, dd$w, 0)
  dd$family <- "binomial"
  dd$fam <- 3L
  dd
}

# the rows of dd with positive weight, for lmstart
ddrows <- function(dd, keep){
  dd$x <- dd$x[keep]
  dd$y <- dd$y[keep]
  dd$w <- dd$w[keep]
  dd$bn <- dd$bn[keep]
  dd$Xs <- lapply(dd$Xs, function(X) X[keep, , drop=FALSE])
  dd
}

# probabilities of the intervals of event data, f(end) - f(start) with f(Inf) = 1, at the curves of the
# observations (random=TRUE with the fitted random effects)
eventprob <- function(object, random=TRUE){
  dat <- object$obj$env$data
  d <- object$data
  fs <- predict(object, random=random)
  tv <- all.vars(object$formula[[3]])
  de <- d
  de[[tv[1]]] <- ifelse(dat$xinf == 1, 0, dat$xe)
  fe <- predict(object, newdata=de, random=random)
  fe[dat$xinf == 1] <- 1
  fe - fs
}

# whether the series of an event fit are complete (see eventcomplete)
eventcompletefit <- function(object){
  dat <- object$obj$env$data
  eventcomplete(dat$x, ifelse(dat$xinf == 1, Inf, dat$xe), object$series)
}

# residuals of event data: response residuals, observed minus expected proportions of the series
# totals N; pearson residuals (count - N p)/sqrt(N p (1 - p)); randomised quantile residuals of the
# binomial distribution of each count, with N trials and the probability p of its interval
eventresid <- function(object, type, random){
  if (!eventcompletefit(object)) stop("residuals of event data need complete series: intervals of each series from 0 to Inf without gaps; see the argument series of drmTMB.", call.=FALSE)
  y <- as.vector(object$obj$env$data$y)
  p <- eventprob(object, random)
  N <- ave(y, object$series, FUN=sum)
  if (type == "response") return(y/N - p)
  if (type == "pearson") return((y - N*p)/sqrt(N*p*(1 - p)))
  qnorm(runif(length(y), pbinom(y - 1, N, p), pbinom(y, N, p)))
}
