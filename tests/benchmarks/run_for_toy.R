library(dplyr)
library(ggplot2)
library(gridExtra)
library(streambugs)

benchmarks.basename <- file.path("output","states_vs_time")

# handler to suppress specific warnings from the streambugs library
.streambugs.suppress.warning.handler <- function(w) {
    if(any(grepl("(no (envcond|parameters) for .*|limitation by (\\w)+( \\w+)? not considered)", w))) {
        invokeRestart("muffleWarning")
    }
}

.count.n.var <- function(sys.def) length(sys.def$y.names$y.names)
# FIXME: always constant (56) when scaling reaches and habitats, and it should scale linearly
.count.n.interact <- function(sys.def) {
    par.stoich = c(unlist(sys.def$par.stoich.taxon), unlist(sys.def$par.stoich.web))
    length(par.stoich)
}

.update.list <- function(lst1, lst2) {
    keys <- unique(c(names(lst1), names(lst2)))
    as.list(stats::setNames(
        mapply(
            function(x,y) if (is.null(y)) x else y,
            lst1[keys], lst2[keys]),
        keys))
}

#' Time `n.repeats` times `run.streambugs` (C version).
#'
#' @param n.repeats number of re-runs of `run.streambugs`
#' @param streambugs.model.fun model definition function in the list format as
#'    returned e.g. by \code{\link{streambugs.example.model.toy}
#' @param ... parameter of the \code{streambugs.model.fun} function
#'
#' @return Data frame with \code{n.repeats} rows and columns:\describe{
#'      \item{\code{"n.var"}}{number of ODE variables (constant over rows)}
#'      \item{\code{"n.interact"}}{number of interactions between ODE variables
#'          (constant over rows)}
#'      \item{\code{"time"}}{ellapsed simulation time}
#'    }
.time.model <- function(n.repeats, streambugs.model.fun,
    streambugs.model.fun.args=list(), run.streambugs.args=list()) {
    model <- do.call(streambugs.model.fun, streambugs.model.fun.args)

    run.args <- .update.list(list(
            y.names  = model$y.names,
            times    = model$times,
            par      = model$par,
            inp      = model$inp,
            C        = TRUE
        ), run.streambugs.args)

    # Note: warning handler does not significantly affect timing
    withCallingHandlers({
        time <- replicate(n.repeats, system.time({
            res.C <- do.call(run.streambugs, run.args)
        })["elapsed"])
        sys.def <- streambugs.get.sys.def(model$y.names, model$par, model$inp)
    }, warning = .streambugs.suppress.warning.handler)
    return(data.frame(
        n.var=.count.n.var(sys.def),
        n.interact=.count.n.interact(sys.def),
        time=time)
    )
}

.time.toy <- function(n.repeats, ...)
    .time.model(n.repeats, streambugs.example.model.toy, ...)

benchmark.states.vs.time <- function(n.repeats.per.run=10, n.Reach.max=10, n.Hab.max=10, ...) {
    # burn-in (lib loading etc)
    .time.toy(3)

    n.Reach.vec <- seq(1,n.Reach.max,by=1)
    n.Hab.vec <- seq(1,n.Hab.max,by=1)
    grid <- expand.grid("n.Reach"=n.Reach.vec, "n.Hab"=n.Hab.vec)
    tdf.list <- apply(grid, 1,
        function(n) data.frame(t(n),
            .time.toy(n.repeats.per.run, streambugs.model.fun.args=list(
                n.Reaches=n["n.Reach"], n.Habitats = n["n.Hab"]
            ), run.streambugs.args=list(...))
        )
    )
    tdf <- Reduce(rbind, tdf.list, NULL)
    return(tdf)
}


.get.benchmarks.file.name <- function(method, ext, basename=benchmarks.basename) {
  paste0(benchmarks.basename, "_", method, ext)
}
.get.benchmarks.csv.name <- function(method, ...) {
    .get.benchmarks.file.name(method, ".csv", ...)
}
.get.benchmarks.pdf.name <- function(method, ...) {
    .get.benchmarks.file.name(method, ".pdf", ...)
}


# TODO: benchmark w/ pracma and w/ cOde "cvode" implementation
# cf. https://cran.r-project.org/web/views/DifferentialEquations.html
# pracma:
# ode23      ode23s     ode45      ode78      odregress
#ode23(f, t0, tf, y0, rtol=1e-5, atol=1e-10)
# rk4 rk4sys
#rk4sys(f, a, b, y0, n)

methods.all = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk", "euler",
    "rk4", "ode23", "ode45", "radau",  "bdf", "bdf_d", "adams", "impAdams",
    "impAdams_d")
#methods.todo = c("ode23", "lsode", "lsodes", "lsodar", "daspk", "euler", "bdf_d",
#    "adams", "impAdams", "impAdams_d")
methods.none = c()
# TODO: select fastest and compare perf. for scaling invert and algae
methods.fastest = c("lsoda", "lsodar", "rk4", "ode23", "bdf_d", "adams", "impAdams_d")
# Adjust to run benchmarks
methods = methods.all
for (method in methods) {
  tdf <- benchmark.states.vs.time(n.repeats.per.run=10, n.Reach.max=10, n.Hab.max=10, method=method)
  benchmarks.csv.name <- .get.benchmarks.csv.name(method)
  write.csv(tdf, benchmarks.csv.name, row.names = FALSE)
}




.plot.states.vs.time <- function(tdf, plot.title = NULL, color.factor=NULL) {
  ttbl <- tibble::as_tibble(tdf)
  #ttbl.gReachHab <- ttbl %>% group_by(n.Reach,n.Hab,n.var,n.interact) %>% summarise_at("time", funs(mean, sd))

  p0 <- ggplot(ttbl, aes(x=n.var, y=time))
  p <- p0 + theme_bw() + xlab("#ODE variables") + ylab("time (s)")
  if (is.null(color.factor)) {
      p <- p +
          geom_point() + # Points
          geom_smooth(method=lm, formula = y ~ x + I(x^2), se=FALSE, fullrange=FALSE) # Add linear regression line
  } else if (color.factor == "n.Reach") {
      p <- p +
          geom_point(aes(color=factor(n.Reach))) + # Points
          scale_color_discrete(name="#Reaches") +  # Legend name (opt: labels)
          geom_smooth(aes(color=factor(n.Reach)), method=lm, formula = y ~ x + I(x^2), se=FALSE, fullrange=FALSE) # Add linear regression line
  } else if (color.factor == "n.Hab") {
      p <- p +
          geom_point(aes(color=factor(n.Hab))) + # Points
          scale_color_discrete(name="#Habitats") +  # Legend name (opt: labels)
          geom_smooth(aes(color=factor(n.Hab)), method=lm, formula = y ~ x + I(x^2), se=FALSE, fullrange=FALSE) # Add quadratic regression line
  } else {
    stop("Invalid value of the `color.factor` argument")
  }
  if (!is.null(plot.title)) p <- p + ggtitle(plot.title)
  return(p)
}

write.benchmarks.plot <- function(tdf, pdf.name, ...) {
    p0 <- .plot.states.vs.time(tdf, ...)
    p1 <- .plot.states.vs.time(tdf, color.factor = "n.Reach", ...)
    p2 <- .plot.states.vs.time(tdf, color.factor = "n.Hab", ...)
    pdf(pdf.name, width=8, height=3*6, useDingbats=FALSE)
    grid.arrange(p0, p1, p2, nrow=3, ncol=1)
    dev.off()
}


for (method in methods) {
  benchmarks.csv.name <- .get.benchmarks.csv.name(method)
  tdf <- read.csv(benchmarks.csv.name, header=TRUE)
  # TODO pdf.basename re-factor
  benchmarks.pdf.name <- .get.benchmarks.pdf.name(method)
  write.benchmarks.plot(tdf, benchmarks.pdf.name,
    plot.title = "vary #Reaches and #Habitats")
}
