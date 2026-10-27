#' Diagnostic plots for drmTMB models
#'
#' Residuals against the fitted values, residuals against the dose, a normal QQ plot of the
#' residuals, and normal QQ plots of the predicted random effects, drawn as ggplot2 graphs on one
#' page.
#'
#' The default randomised quantile residuals are standard normal under the model for all families,
#' so the QQ plot is compared with the identity line; for the discrete families they are drawn
#' uniformly within the jump of the distribution function at each observation, and \code{seed}
#' makes the plots reproducible. For \code{type="response"}, the QQ plot line goes through the
#' quartiles. A lowess smoother (red) shows trends of the residuals against the fitted values and
#' the dose.
#'
#' For the models on the log-dose scale, the dose axis is logarithmic; observations at dose zero are
#' drawn left of the smallest positive dose, separated by a dashed line.
#'
#' The QQ plots of the random effects (\code{which=4}) show the predicted random effects of each
#' parameter and grouping factor. They are shrunk towards zero, so their spread is smaller than the
#' estimated standard deviation, and with few groups they say little about the normality of the
#' random effects.
#'
#' @param x a drmTMB object
#' @param which the plots to draw: 1 residuals against fitted values, 2 residuals against dose, 3
#'   normal QQ plot of the residuals, 4 normal QQ plots of the predicted random effects (models with
#'   random effects only)
#' @param type residual type, see \code{\link{residuals.drmTMB}}: "quantile" (default), "pearson" or
#'   "response"
#' @param random if TRUE, residuals and fitted values include the predicted random effects
#' @param logdose logarithmic dose axis; by default for the models on the log-dose scale
#'   (loglogistic, weibull1, weibull2, lognormal)
#' @param ncol number of columns of the page
#' @param seed seed for the randomised quantile residuals
#' @param draw if FALSE, the graphs are returned without drawing them
#' @param ... not used
#'
#' @return a list of ggplot objects, invisibly, named fitted, dose, qq and ranef, for changing or
#'   drawing single graphs
#' @export
plot.drmTMB <- function(x, which=1:3, type=c("quantile", "pearson", "response"), random=TRUE, logdose=NULL, ncol=NULL, seed=NULL, draw=TRUE, ...){
  type <- match.arg(type)
  if (!all(which %in% 1:4)) stop("which must be a subset of 1:4.")
  if (4 %in% which && length(x$reterms) == 0){
    if (length(which) == 1) stop("the model has no random effects.")
    which <- setdiff(which, 4)
  }
  if (is.null(logdose)) logdose <- modelindex(x$model) %in% 2:5
  if (!is.null(seed)) set.seed(seed)

  res <- residuals(x, type=type, random=random)
  d <- data.frame(fitted=fitted(x, random=random), dose=as.vector(x$obj$env$data$x), res=res)
  ylab <- switch(type, quantile="Quantile residuals", pearson="Pearson residuals", response="Residuals")

  plots <- list()
  if (1 %in% which){
    plots$fitted <- residplot(d, "fitted", ylab) + ggplot2::labs(x="Fitted values", title="Residuals vs fitted")
  }
  if (2 %in% which){
    plots$dose <- residplot(d, "dose", ylab, logdose=logdose) + ggplot2::labs(x="Dose", title="Residuals vs dose")
  }
  if (3 %in% which){
    plots$qq <- qqplot(data.frame(res=res), identity=type != "response") +
      ggplot2::labs(x="Theoretical quantiles", y=ylab, title="Normal QQ plot")
  }
  if (4 %in% which){
    plots$ranef <- qqplot(ranefdata(x), identity=FALSE) +
      ggplot2::facet_wrap(~ term, scales="free_y") +
      ggplot2::labs(x="Theoretical quantiles", y="Predicted random effects", title="Random effects")
  }
  if (draw) drawplots(plots, ncol)
  invisible(plots)
}

# residuals against a variable, with a lowess smoother; zero doses left of the positive doses on a
# log axis
residplot <- function(d, var, ylab, logdose=FALSE){
  d$xv <- d[[var]]
  zero <- logdose && any(d$xv == 0)
  if (logdose){
    xp <- d$xv[d$xv > 0]
    if (length(xp) == 0) stop("all doses are zero; use logdose=FALSE.")
    # dose zero at a gap of a fifth of the range of the positive doses, at least half a decade
    gap <- max(0.2*diff(log10(range(xp))), 0.5)
    x0 <- 10^(log10(min(xp)) - gap)
    d$xv[d$xv == 0] <- x0
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x=xv, y=res)) +
    ggplot2::geom_hline(yintercept=0, colour="grey50") +
    ggplot2::geom_point(alpha=0.7)
  if (length(unique(d$xv)) >= 3){
    lw <- as.data.frame(lowess(if (logdose) log10(d$xv) else d$xv, d$res))
    if (logdose) lw$x <- 10^lw$x
    if (zero) lw <- lw[lw$x > x0, ]
    p <- p + ggplot2::geom_line(data=lw, ggplot2::aes(x=x, y=y), colour="red")
  }
  if (logdose){
    ux <- sort(unique(xp))
    br <- if (length(ux) <= 8) ux else 10^(floor(log10(min(ux))):ceiling(log10(max(ux))))
    lab <- format(br, trim=TRUE, drop0trailing=TRUE)
    if (zero){
      br <- c(x0, br)
      lab <- c("0", lab)
      p <- p + ggplot2::geom_vline(xintercept=sqrt(x0*min(xp)), linetype="dashed", colour="grey50")
    }
    p <- p + ggplot2::scale_x_log10(breaks=br, labels=lab)
  }
  p + ggplot2::labs(y=ylab) + ggplot2::theme_bw()
}

# normal QQ plot of res, with the identity line or a line through the quartiles; by term if d has one
qqplot <- function(d, identity=TRUE){
  grp <- if (is.null(d$term)) rep(1, nrow(d)) else d$term
  d$theo <- NA
  for (g in unique(grp)){
    i <- which(grp == g)
    d$theo[i[order(d$res[i])]] <- qnorm(ppoints(length(i)))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x=theo, y=res))
  if (identity){
    p <- p + ggplot2::geom_abline(intercept=0, slope=1, colour="red")
  } else {
    ql <- do.call(rbind, lapply(unique(grp), function(g){
      r <- d$res[grp == g]
      slope <- diff(quantile(r, c(0.25, 0.75), names=FALSE)) / diff(qnorm(c(0.25, 0.75)))
      data.frame(term=g, slope=slope, intercept=quantile(r, 0.25, names=FALSE) - slope*qnorm(0.25))
    }))
    if (is.null(d$term)) ql$term <- NULL
    p <- p + ggplot2::geom_abline(data=ql, ggplot2::aes(intercept=intercept, slope=slope), colour="red")
  }
  p + ggplot2::geom_point(alpha=0.7) + ggplot2::theme_bw()
}

# predicted random effects, one panel for each parameter (and covariate of random slopes) and term
ranefdata <- function(object){
  u <- object$sdrl$par.random
  uoff <- 0
  do.call(rbind, lapply(object$reterms, function(tm){
    nl <- length(tm$levels)
    out <- do.call(rbind, lapply(seq_along(tm$pars), function(k){
      data.frame(term=paste(relabels(tm)[k], "|", tm$group), res=u[uoff + (k - 1)*nl + seq_len(nl)])
    }))
    uoff <<- uoff + nl*length(tm$pars)
    out
  }))
}

# draw ggplot objects on one page, in ncol columns
drawplots <- function(plots, ncol=NULL){
  n <- length(plots)
  if (is.null(ncol)) ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n/ncol)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow, ncol)))
  for (i in seq_len(n)){
    print(plots[[i]], vp=grid::viewport(layout.pos.row=(i - 1) %/% ncol + 1, layout.pos.col=(i - 1) %% ncol + 1))
  }
  grid::popViewport()
}

# columns of the plot data used in the aesthetics
utils::globalVariables(c("xv", "res", "x", "y", "theo", "intercept", "slope"))
