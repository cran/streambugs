################################################################
#
# streambugs 1.2
# ==============
#
# -------------------
# auxiliary functions
# -------------------
#
# creation:      06.04.2020
# modifications: 08.04.2020
#
################################################################

# get stoichiometric coefficients for "taxon-based" processes:
# ------------------------------------------------------------

# Returns a list of stoichiometric vectors of processes of the type given
# by the character string "proc".
# The list elements are labelled according to the taxon labelling the process,
# the stoichiometric vectors contain stoichiometric coefficients of the taxa/
# substances indicated by their labels.

get.par.stoich.taxon <- function(proc,par)
{
  par.names <- names(par)
  par.splitted <- strsplit(par.names,split="_")
  stoich <- list()
  for ( i in 1:length(par) )
  {
    if ( par.splitted[[i]][1] == proc )
    {
      if ( length(par.splitted[[i]]) < 3 )
      {
        warning("Parameter name \"",par.names[i],"\" starts with \"",
                proc,"\""," but does not have three components")
      }
      else
      {
        if ( length(stoich[[par.splitted[[i]][2]]]) == 0 )
        {
          stoich[[par.splitted[[i]][2]]] <- numeric(0)
        }
        stoich[[par.splitted[[i]][2]]][par.splitted[[i]][3]] <- par[i]
      }
    }
  }
  return(stoich)
}


# get stoichiometric coefficients for food web processes:
# -------------------------------------------------------

# Returns a list of lists of stoichiometric vectors of food web processes
# of the type given by the character string "proc".
# The list elements are labelled according to the lead taxon labelling the
# process (typically the consumer), the list elements of the lower-level
# lists are labelled according to the second taxon characterizing the process
# (typically the food). The stoichiometric vectors contain stoichiometric
# coefficients of the taxa/substances indicated by their labels.

get.par.stoich.web <- function(proc,par)
{
  par.names <- names(par)
  par.splitted <- strsplit(par.names,split="_")
  stoich <- list()
  for ( i in 1:length(par) )
  {
    if ( par.splitted[[i]][1] == proc )
    {
      if ( length(par.splitted[[i]]) < 4 )
      {
        warning("Parameter name \"",par.names[i],"\" starts with \"",
                proc,"\""," but does not have four components")
      }
      else
      {
        if ( length(stoich[[par.splitted[[i]][2]]]) == 0 )
        {
          stoich[[par.splitted[[i]][2]]] <- list()
        }
        if ( length(stoich[[par.splitted[[i]][2]]][[par.splitted[[i]][3]]]) == 0 )
        {
          stoich[[par.splitted[[i]][2]]][[par.splitted[[i]][3]]] <- numeric(0)
        }
        stoich[[par.splitted[[i]][2]]][[par.splitted[[i]][3]]][par.splitted[[i]][4]] <- par[i]
      }
    }
  }
  return(stoich)
}


# get process stoichiometries list:
# ---------------------------------

# get process stoichiometries list from a list of kinetic parameters using kinetic
# parameters names vector and corresponding single param query function:
# get.par.stoich.taxon or get.par.stoich.web

get.par.stoich.list <- function(par, kinparnames.vec, get.par.stoich.func) {
  # sapply(..., simplify = FALSE) instead of lapply(...) to pass `USE.NAMES = TRUE`
  sapply(kinparnames.vec, get.par.stoich.func, par, simplify = FALSE, USE.NAMES = TRUE)
}


# add state variable index row to stoichiometric matrix:
# ------------------------------------------------------

# for a given reach and habitat complement a given stoichiometric matrix row with
# indices of corresponding state variables as indicated by row

stoich.add.ind.statevar.row <- function(y.names, reach, habitat, stoichvec) {
  n <- length(stoichvec)
  stoichmat <- matrix(NA,nrow=2,ncol=n)
  colnames(stoichmat) <- names(stoichvec)
  rownames(stoichmat) <- c("statevar","coeff")
  stoichmat[2,] <- stoichvec
  for ( j in 1:n )
  {
    taxon <- names(stoichvec)[j]
    ind.statevar <-
      y.names$y.taxa == taxon &
      y.names$y.reaches == reach &
      y.names$y.habitats == habitat
    if ( sum(ind.statevar) > 1 )
    {
      warning("state vector contains multiple entries for",
              "(Reach,Habitat,Taxon) = (",reach,",",habitat,",",taxon,")", sep="")
    }
    else
    {
      if ( sum(ind.statevar) == 0 )
      {
        warning("state vector has no entry for",
                "(Reach,Habitat,Taxon) = (",reach,",",habitat,",",taxon,")", sep="")
      }
      else
      {
        stoichmat[1,j] <- which(ind.statevar)
      }
    }
  }
  return(stoichmat)
}


# get kinetic parameters of taxon-based processes:
# ------------------------------------------------

get.proc.inpinds.mask.and.offset <- function(kinparnames, inpinds.j) {
  kinparnames.length <- sapply(kinparnames, length)
  kinparnames.cumsum <- cumsum(kinparnames.length)
  kinparnames.offset <- kinparnames.cumsum - kinparnames.length

  # match inpinds for a given proc working under assumption that inpinds matrix was
  # obtained using `unlist(kinparamnames)`, i.e. using cumulative sizes of queried
  # parameters sub-vectors
  proc.inpinds.mask <- outer(1:length(kinparnames), inpinds.j,
                             function(i, j.inpinds) {
                               kinparnames.offset[i] < j.inpinds & j.inpinds <= kinparnames.cumsum[i]
                             })
  rownames(proc.inpinds.mask) <- names(kinparnames)
  return(list(mask=proc.inpinds.mask, offset=kinparnames.offset))
}

get.i.proc.inpinds <- function(inpinds, i.proc.inpinds.mask, i.proc.inpinds.offset = 0) {
  i.proc.inpinds <- inpinds[i.proc.inpinds.mask, c("inpind", "j"), drop=FALSE]
  colnames(i.proc.inpinds)[2] <- "i"
  i.proc.inpinds[, "i"] <- i.proc.inpinds[, "i"] - i.proc.inpinds.offset
  i.proc.inpinds
}

get.i.proc.inpinds.parvals.list <- function(
  inpind.parval,
  i.y.name,  # single row index or row (var) name
  proc.par.names,
  i.proc.inpinds.mask,
  i.proc.inpinds.offset)
{
  i.proc.par.inpinds.parvals.list <- list(inpinds=NULL, parvals=NULL)

  # set $parvals
  i.proc.parvals <- if (is.na(i.y.name)) {
    rep(NA, length(proc.par.names))
  } else {
    inpind.parval$parvals[i.y.name, proc.par.names]
  }
  names(i.proc.parvals) <- proc.par.names
  i.proc.par.inpinds.parvals.list$parvals <- i.proc.parvals

  # set $inpinds
  i.proc.par.inpinds.parvals.list$inpinds <- get.i.proc.inpinds(
    inpinds = inpind.parval$inpinds,
    i.proc.inpinds.mask = i.proc.inpinds.mask,
    i.proc.inpinds.offset = i.proc.inpinds.offset)

  i.proc.par.inpinds.parvals.list
}

get.par.proc.taxon <- function(y.names, par, inp, kinparnames.taxon, par.stoich.taxon) {

  # 1. Query parameters and inputs
  inpind.parval <- get.inpind.parval.taxon.group.reach.habitat(
    par.names = unlist(kinparnames.taxon),
    y.names = y.names,
    par = par,
    inp = inp,
    defaults = c(fgrotax=1)
  )

  # 2. Prepare result list (re-shaping)
  par.proc.taxon <- list()

  inpinds.j.mask.and.offset <- get.proc.inpinds.mask.and.offset(
    kinparnames.taxon, inpind.parval$inpinds[, "j"])

  for ( i in 1:length(y.names$y.names) )
  {
    i.par.proc.taxon <- list()

    i.reach <- y.names$y.reaches[i]
    i.habitat <- y.names$y.habitats[i]

    i.inpinds.mask.i <- (inpind.parval$inpinds[,"i"] == i)

    for ( proc.name in names(kinparnames.taxon) )
    {

      proc.par.names <- kinparnames.taxon[[proc.name]]
      proc.par.stoich <- par.stoich.taxon[[proc.name]]

      proc.inpinds.j.mask <- inpinds.j.mask.and.offset$mask[proc.name,]
      proc.inpinds.j.offset <- inpinds.j.mask.and.offset$offset[proc.name]

      # taxon-based processes:
      stoich <- proc.par.stoich[[y.names$y.taxa[i]]]
      if ( length(stoich) > 0 ) {

        i.proc.inpinds.parvals.list <- get.i.proc.inpinds.parvals.list(
          inpind.parval,
          i,
          proc.par.names,
          i.inpinds.mask.i & proc.inpinds.j.mask,
          proc.inpinds.j.offset)

        check.inpind.parval.required(i.proc.inpinds.parvals.list,
          required = proc.par.names,
          common.no.val.suffix = paste0(
            " for state variable #", i, ' and process "', proc.name, '"'))

        # set $stoich
        i.proc.inpinds.parvals.list$stoich <- stoich.add.ind.statevar.row(
          y.names, i.reach, i.habitat, stoich)

        i.par.proc.taxon[[proc.name]] <- i.proc.inpinds.parvals.list

      }

    }

    par.proc.taxon[[i]] <- i.par.proc.taxon
  }

  names(par.proc.taxon) <- y.names$y.names
  return(par.proc.taxon)
}


# get kinetic parameters of food web processes:
# ---------------------------------------------

get.ind.1sthabinreach <- function(ind.fA) {
  first.ind.fA.ind.1sthab <- ind.fA[[1]]$ind.1sthab
  ind.1sthabinreach <- sapply(
    ind.fA, function(i.ind.fA) i.ind.fA$ind.reach[first.ind.fA.ind.1sthab])
  dim(ind.1sthabinreach) <- NULL
  names(ind.1sthabinreach) <- NULL  # Note: sapply(..., USE.NAMES = FALSE) will not work
  ind.1sthabinreach
}

.group.int.values <- function(x) {
  dup.x <- duplicated(x)
  l <- list()
  for (k in seq_along(x)) {
    k.x <- x[k]
    k.l <- if(dup.x[k]) l[[k.x]] else NULL  # Note: ifelse won't work here with NULL
    l[[k.x]] <- append(k.l, k)
  }
  return(l)
}

group.inpinds.by.i <- function(inpind.parval) {

  inpinds.list <- .group.int.values(inpind.parval$inpinds[,"i"])
  n <- length(inpinds.list)
  if (n > 0) {
    names(inpinds.list) <- rownames(inpind.parval$parvals)[1:n]
    inpinds.list <- inpinds.list[!sapply(inpinds.list, is.null)]
  }

  # equivalent dplyr variant (w/o own helper function and omitting non-shows)
  #
  # inpinds <- inpind.parval$inpinds
  # inpinds.group_data <-
  #   dplyr::group_by(data.frame(inpinds), i) %>%
  #   dplyr::group_data()
  # inpinds.list <- as.list(t(inpinds.group_data$.rows))
  # names(inpinds.list) <- rownames(inpind.parval$parvals)[inpinds.group_data$i]

  return(inpinds.list)
}

.names.to.indices <- function(names.vec) {
  list2env(
    stats::setNames(
      object = as.list(1:length(names.vec)),
      nm = names.vec),
    hash = TRUE
  )
}

get.par.proc.web <- function(y.names, par, inp, kinparnames.web, par.stoich.web) {

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  n0 <- length(y.names$y.names)

  # 1st in habitat and reach get second taxa: "Fish"
  ind.1sthabinreach <- get.ind.1sthabinreach(y.names$ind.fA)

  # TODO: re-factor into a function
  # extend taxa for "Fish" incl. group, reaches and habitats for querying as well as
  # y.names for further processing
  y.names$taxa <- c(y.names$taxa, "Fish")
  # y.names$y.taxa <- factor(c(y.names$y.taxa, rep("Fish", length(ind.1sthabinreach))), levels=y.names$taxa)
  y.names$y.taxa <- c(y.names$y.taxa, rep("Fish", length(ind.1sthabinreach)))
  ind.y.names.ext <- c(1:n0, ind.1sthabinreach)
  # y.names$y.groups <- factor(y.names$y.groups[ind.y.names.ext], levels=y.names$groups)
  y.names$y.groups <- y.names$y.groups[ind.y.names.ext]
  # y.names$y.reaches <- factor(y.names$y.reaches[ind.y.names.ext], levels=y.names$reaches)
  y.names$y.reaches <- y.names$y.reaches[ind.y.names.ext]
  # y.names$y.habitats <- factor(y.names$y.habitats[ind.y.names.ext], levels=y.names$habitats)
  y.names$y.habitats <- y.names$y.habitats[ind.y.names.ext]
  # beware: prefixing "Fish" taxa y.names to perserve unique var names
  y.names$y.names <- c(y.names$y.names, paste0("Fish_", y.names$y.names[ind.1sthabinreach]))

  n1 <- length(y.names$y.names)

  # 1. Query parameters and inputs
  kinparnames.web.taxon1 <- sapply(kinparnames.web,
                                   function(proc.par.names) proc.par.names$taxon1)
  par.names.taxon1 <- unlist(kinparnames.web.taxon1)
  inpind.parval.taxon1 <- get.inpind.parval.taxon.group.reach.habitat(
    par.names = par.names.taxon1,
    y.names = y.names,
    par = par,
    inp = inp,
    defaults = c(fgrotax=1, q=1))
  inpinds.j.taxon1.mask.and.offset <- get.proc.inpinds.mask.and.offset(
    kinparnames.web.taxon1, inpind.parval.taxon1$inpinds[, "j"])

  kinparnames.web.taxa2 <- sapply(kinparnames.web,
                                  function(proc.par.names) proc.par.names$taxa2)
  par.names.taxa2 <- unlist(kinparnames.web.taxa2)
  defaults.taxa2 <- c(Pref=1)
  inpind.parval.taxa2 <- get.inpind.parval.taxon.group.pairs.reach.habitat(
    par.names = par.names.taxa2,
    y.names = y.names,
    par = par,
    inp = inp,
    defaults = defaults.taxa2)
  # performance tweak: use environment as a fast hash map instead of rownames for
  # a string-based row lookup
  parval.taxa2.rowname.to.ind <- .names.to.indices(
    rownames(inpind.parval.taxa2$parvals))
  inpinds.j.taxa2.mask.and.offset <- get.proc.inpinds.mask.and.offset(
    kinparnames.web.taxa2, inpind.parval.taxa2$inpinds[, "j"])
  taxa2.inpinds.ind.list <- group.inpinds.by.i(inpind.parval.taxa2)

  # 2. Prepare result list (re-shaping)
  par.proc.web <- rep(list(list()), n0)

  for ( i in 1:n1 )
  {
    i.par.proc.web <- list()
    i0 <- ind.y.names.ext[i]

    i.reach   <- y.names$y.reaches[i]
    i.habitat <- y.names$y.habitats[i]
    taxon1  <- y.names$y.taxa[i]

    i.inpinds.taxon1.mask.i <- (inpind.parval.taxon1$inpinds[,"i"] == i)

    # logical mask with positions matching same reaches and habitats as the current ones
    i.mask.reach.habitat <- y.names$y.reaches == i.reach & y.names$y.habitats == i.habitat
    i.y.names.reach.habitat <- y.names$y.taxa[i.mask.reach.habitat]

    for ( proc.name in names(kinparnames.web) )
    {
      proc.stoich <- par.stoich.web[[proc.name]]
      proc.par.names <- kinparnames.web[[proc.name]]

      proc.inpinds.j.taxon1.mask <- inpinds.j.taxon1.mask.and.offset$mask[proc.name,]
      proc.inpinds.j.taxon1.offset <- inpinds.j.taxon1.mask.and.offset$offset[proc.name]

      proc.inpinds.j.taxa2.mask <- inpinds.j.taxa2.mask.and.offset$mask[proc.name,]
      proc.inpinds.j.taxa2.offset <- inpinds.j.taxa2.mask.and.offset$offset[proc.name]

      # food web processes:
      stoich <- proc.stoich[[taxon1]]
      # e.g. stoich == list(POM1=c(Invert1=1, POM1=-1), Alga1=c(Invert1=1, Alga1=-1))
      if ( length(stoich) > 0 )
      {

        i.proc.inpinds.parvals.list <- get.i.proc.inpinds.parvals.list(
          inpind.parval.taxon1,
          i,
          proc.par.names$taxon1,
          i.inpinds.taxon1.mask.i & proc.inpinds.j.taxon1.mask,
          proc.inpinds.j.taxon1.offset)

        check.inpind.parval.required(i.proc.inpinds.parvals.list,
          required = proc.par.names$taxon1,
          common.no.val.suffix = paste0(
            " for state variable #", i, ' and process "', proc.name, '"'))

        # TODO: re-factor into a function
        # set $taxa2
        # stoichiometrically-related taxa e.g. c("POM1", "Alga1")
        taxa2 <- names(stoich)

        # indices in parvals pairs rows of stoichiometrically-related state variables
        # matching process, taxon1, reach, habitat, source state var and target taxa2
        i.taxa2.reach.habitat.i <- match(taxa2, i.y.names.reach.habitat)
        # self-loop hack as a `c -> c + N x`  intepretation for input constants
        # representing external food sources and defined only implicitly via parameters
        # specification, hence missing from explicit taxa list (giving `NA` in the above
        # `match`), such as `"SusPOM"` taxon used in
        # `streambugs.example.model.extended()` parameters, e.g. in
        # `"Cons_Chironomidae_SusPOM_Chironomidae"` parameter.
        i.taxa2.reach.habitat.i <- ifelse(
          is.na(i.taxa2.reach.habitat.i), i, i.taxa2.reach.habitat.i)
        i.taxa2.parvals.rownames <-
          paste0(y.names$y.names[i], "->", y.names$y.names[i.taxa2.reach.habitat.i])
        i.taxa2.inpinds.ind <- taxa2.inpinds.ind.list[i.taxa2.parvals.rownames]

        taxa2.inpinds.parvals.list <- list()
        for ( k in 1:length(taxa2) )
        {
          k.i.taxa2.parvals.rowname <- i.taxa2.parvals.rownames[k]
          k.i <- parval.taxa2.rowname.to.ind[[k.i.taxa2.parvals.rowname]]
          k.i.inpinds.taxa2.ind <- i.taxa2.inpinds.ind[[k.i.taxa2.parvals.rowname]]
          # k.i.inpinds.taxa2.mask.i <- (inpind.parval.taxa2$inpinds[,"i"] == k.i)

          k.i.proc.par.proc.taxa2 <- get.i.proc.inpinds.parvals.list(
            inpind.parval = inpind.parval.taxa2,
            # i.y.name = k.i.taxa2.parvals.rowname,
            i.y.name = k.i,
            proc.par.names = proc.par.names$taxa2,
            i.proc.inpinds.mask =
              k.i.inpinds.taxa2.ind[proc.inpinds.j.taxa2.mask[k.i.inpinds.taxa2.ind]],
            # i.proc.inpinds.mask = k.i.inpinds.taxa2.mask.i & proc.inpinds.j.taxa2.mask,
            i.proc.inpinds.offset = proc.inpinds.j.taxa2.offset)

          check.inpind.parval.required(k.i.proc.par.proc.taxa2,
            required = proc.par.names$taxa2,
            common.no.val.suffix = paste0(
              ' for state variables "', k.i.taxa2.parvals.rowname, '"',
              ' and taxon1 process "', proc.name, '"'))

          # set $stoich
          k.i.proc.par.proc.taxa2$stoich <- stoich.add.ind.statevar.row(
            y.names, i.reach, i.habitat, stoich[[k]])

          taxa2.inpinds.parvals.list[[taxa2[k]]] <- k.i.proc.par.proc.taxa2
        }
        i.proc.inpinds.parvals.list$taxa2 <- taxa2.inpinds.parvals.list

        i.par.proc.web[[proc.name]] <- i.proc.inpinds.parvals.list
      }

    }

    par.proc.web[[i0]] <- c(par.proc.web[[i0]], i.par.proc.web)
  }
  names(par.proc.web) <- y.names$y.names[1:n0]

  return(web=par.proc.web)
}
