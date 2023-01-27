################################################################
#
# streambugs 1.2
# ==============
#
# -------------------
# auxiliary functions
# -------------------
#
# creation:      08.09.2012
# modifications: 08.04.2020
#
################################################################

#' Construct the streambugs ODE state variable names
#'
#' Construct encoded labels of streambugs ODE state variable names from reach
#' names, habitat names, and "taxa" names for at least one of the POM, Algae, or
#' Invertebrates groups.
#'
#' @param Reaches reach names character vector; duplicates are dropped
#' @param Habitats habitats names character vector; duplicates are dropped
#' @param POM optional ("taxa") character vector of POM (particulate organic matter) group;
#'   duplicates are dropped
#'   (note: at least one "taxon" name out of all groups has to be given)
#' @param Algae optional taxa character vector of Algae group; duplicates are
#'  dropped (note: at least one taxon name out of all groups has to be given)
#' @param Invertebrates optional taxa character vector of Invertebrates group
#'  duplicates are dropped (note: at least one taxon name out of all groups has
#;  to be given)
#'
#' @return vector with state variable names in the form of
#'    \code{"Reach_Habitat_Taxon_Group"}
#'
#' @examples
#' Reaches           <- paste0("Reach",1:2)
#' Habitats          <- paste0("Hab",1:1)
#' y.names <- construct.statevariables(Reaches,Habitats,Invertebrates=c("Baetis","Ecdyonurus"))
#'
#' @export
construct.statevariables <- function(Reaches,Habitats,POM=NULL,Algae=NULL,Invertebrates=NULL)
{
  # validate that at least one taxon or POM was provided
  if ( is.null(c(POM, Algae, Invertebrates)) )
    stop("Missing at least one \"taxon\" from either group of POM, Algae, or Invertebrates")

  # encoding paste function helper
  paste.enc <- function(...) paste(..., sep="_")

  # rm duplicates
  if ( length(unique(Reaches)) != length(Reaches) ) {
    warning("Non unique reaches names - dropping duplicates")
    Reaches = unique(Reaches)
  }
  if ( length(unique(Habitats)) != length(Habitats) ) {
    warning("Non unique reaches names - dropping duplicates")
    Habitats = unique(Habitats)
  }

  # merge taxa groups
  POM.enc <- if (is.null(POM)) NULL else paste.enc(POM, "POM")
  Algae.enc <- if (is.null(Algae)) NULL else paste.enc(Algae, "Algae")
  Invertebrates.enc <- if (is.null(Invertebrates)) NULL else paste.enc(Invertebrates, "Invertebrates")
  taxa.enc <- c(POM.enc, Algae.enc, Invertebrates.enc)

  # rm duplicates
  if ( length(unique(taxa.enc)) != length(taxa.enc) ) {
    warning("Non unique taxa names within a group - dropping duplicates")
    taxa.enc = unique(taxa.enc)
  }

  # encode state variables: (reach_i, habitat_i) \times taxon_j
  n.Reaches = length(Reaches)
  n.Habitats = length(Habitats)
  n.Taxa = length(taxa.enc)
  y.names <- character(n.Reaches*n.Habitats*n.Taxa)
  for ( i in 1:n.Reaches )
    for ( j in 1:n.Habitats )
      y.names[(i-1)*n.Habitats*n.Taxa+ (j-1)*n.Taxa + 1:n.Taxa] <-
        paste.enc(Reaches[i], Habitats[j], taxa.enc)

  # sanity check: all values set
  if ( any(y.names == "") )
    stop("problem constructing state variables")

  return(y.names)
}

# extract reach names, habitat names, taxa names and optional
# group names from labels of a state vector:
# -----------------------------------------------------------

#' Decode the streambugs ODE state variable names
#'
#' Extract reach names, habitat names, taxa names and optional group names from
#' encoded labels of streambugs ODE state variable names.
#'
#' @param y.names vector with state variable names in the form of
#'    \code{"Reach_Habitat_Taxon"} or \code{"Reach_Habitat_Taxon_Group"}
#'
#' @return List with:\describe{
#'    \item{\code{$y.names}}{names of state variables (input argument)}
#'    \item{\code{$y.reaches}, \code{$y.reaches}, \code{$y.habitats},
#'      \code{$y.taxa}, and \code{$y.groups}:}{Names of, respectively, reaches,
#'      habitats, taxa, and of groups of each state variable}
#'    \item{\code{$reaches}, \code{$habitats}, \code{$taxa}, \code{$groups}:}{
#'      Unique names of, respectively, reaches, habitats, taxa, and of groups of
#'      state variables}
#'    \item{\code{$ind.fA}:}{Indices used for the areal fractions of each reach
#'      and habitat}
#'    }
#'
#' @examples
#' y.names <- c("Reach1_Hab1_Baetis_Invertebrates","Reach1_Hab1_Ecdyonurus_Invertebrates",
#'              "Reach2_Hab1_Baetis_Invertebrates", "Reach2_Hab1_Ecdyonurus_Invertebrates")
#' decode.statevarnames(y.names)
#'
#' @export
decode.statevarnames <- function(y.names)
{
  # split names of state variables:

  y.split    <- strsplit(y.names,split="_")
  y.reaches  <- rep(NA,length(y.names))
  y.habitats <- rep(NA,length(y.names))
  y.taxa     <- rep(NA,length(y.names))
  y.groups   <- rep(NA,length(y.names))
  for ( i in 1:length(y.names) )
  {
    y.reaches[i]  <- y.split[[i]][1]
    y.habitats[i] <- y.split[[i]][2]
    y.taxa[i]     <- y.split[[i]][3]
    y.groups[i]   <- y.split[[i]][4]
  }

  # determine reach, habitat, taxon and group names:

  reaches  <- unique(y.reaches[!is.na(y.reaches)])
  habitats <- unique(y.habitats[!is.na(y.habitats)])
  taxa     <- unique(y.taxa[!is.na(y.taxa)])
  groups   <- unique(y.groups[!is.na(y.groups)])

  # issue warnings:

  reaches.na <- which(is.na(y.reaches))
  if ( length(reaches.na)  > 0 )
  {
    warning("Undefined reach name(s) in component(s) of state vector: ",
            paste(reaches.na,collapse=","))
  }
  habitats.na <- which(is.na(y.habitats))
  if ( length(habitats.na)  > 0 )
  {
    warning("Undefined habitat name(s) in component(s) of state vector: ",
            paste(habitats.na,collapse=","))
  }
  taxa.na <- which(is.na(y.taxa))
  if ( length(taxa.na)  > 0 )
  {
    warning("Undefined taxon name(s) in component(s) of state vector: ",
            paste(taxa.na,collapse=","))
  }
  y.groups[is.na(y.groups)] <- ""

  # calculate indices for normalization of habitat area fractions:

  ind.fA <- list()
  for ( reach in reaches )
  {
    ind.reach  <- which(y.reaches == reach)
    habs.reach <- unique(y.habitats[ind.reach])
    ind.1sthab <- match(habs.reach,y.habitats[ind.reach])

    ind.fA[[reach]] <- list()
    ind.fA[[reach]][["ind.reach"]]  <- ind.reach
    ind.fA[[reach]][["ind.1sthab"]] <- ind.1sthab
  }

  # return splitted state variable names and reach, habitat, taxa and group names:

  return(list(y.names     = y.names,
              y.reaches   = y.reaches,
              y.habitats  = y.habitats,
              y.taxa      = y.taxa,
              y.groups    = y.groups,
              reaches     = reaches,
              habitats    = habitats,
              taxa        = taxa,
              groups      = groups,
              ind.fA      = ind.fA))
}


# get values of global parameters:
# --------------------------------

# The function searches for time-dependent inputs or constant values
# of global parameters with names specified in the character vector
# "par.names".
# It returns a list with a 2-column matrix of input and value indices
# in its columns for the parameter names that occur as input elements
# (values of those components will need updating as a function of time)
# and a vector of numerical values of the occrrence in the parameter
# vector or in the defaults vector.
# The function issues warnings for parameters that do neither occur
# in the input, nor in the parameter vector and do not have default
# values.

get.inpind.parval.global <- function(par.names,par,inp=NA,
                                     required=NA,defaults=NA)
{
  # get input indices:

  inds <- matrix(0,nrow=0,ncol=2)       # matrix with indices of input and
  # of values to be updated by inputs
  if ( is.list(inp) )
  {
    inpind <- match(par.names,names(inp))
    inpind <- inpind[!is.na(inpind)]
    if(length(inpind)>0)
    {
      i <- match(names(inp)[inpind],par.names)
      inds <- matrix(nrow=length(inpind),ncol=2,c(inpind,i),byrow=F)
    }
  }
  colnames(inds) <- c("inpind","i")

  # get values:

  vals <- par[par.names]
  names(vals) <- par.names

  # if not available set values to defaults (if available):

  ind.na <- which(is.na(vals))
  if ( length(ind.na) > 0 )
  {
    vals[ind.na] <- defaults[par.names[ind.na]]
  }

  inds.vals <- list(inpinds=inds,parvals=vals)

  # check existence of required parameters:

  check.inpind.parval.required(inds.vals, required = required)

  # return list of matrix of input indices and vector of parameter values:

  return(inds.vals)
}

get.inpind.parval.global.envcondtraits <- function(trait.names,par,inp=NA,
                                                   required=NA,defaults=NA)

  # for global parameters that are related to a trait, where we have an unknown number of classes
  # calls get.inpind.parval.global,
  # returns a named list with one entry for each trait and a vector with corresponding parameters

{

  # search par.names for trait.names

  res.trait <- list()

  for (i in 1:length(trait.names))
  {
    par.names <- character(0)

    ind1 <- grep(paste(trait.names[i],"_",sep=""),names(par))

    if(length(ind1)>0)
    {
      ind2 <- which(trait.names[i]==unlist(strsplit(names(par[ind1]),split="_")))

      par.names <- c(par.names,
                     paste(trait.names[i],unique(unlist(strsplit(names(par[ind1]),
                                                                 split="_"))[ind2+1]),
                           sep="_") )
    } else
    {
      warning("no envcond for traitclass ", trait.names[i], " found."  )
    }


    if (length(par.names)>0)
    {
      res <- get.inpind.parval.global(par.names=par.names,par=par,inp=inp,
                                      required=required,defaults=defaults)

      # sort vals and inds

      ind.order <- order(res$parvals)

      res$parvals <- res$parvals[ind.order]

      if(nrow(res$inpinds)>0)
      {
        for(j in 1:nrow(res$inpinds))
        {
          i.new <-  which(res$inpinds[j,2]==ind.order)
          res$inpinds[j,2] <- i.new
        }
      }

      # return list of matrices with input indices and values
      # of all parameters for all state variables:
    } else
    {
      inds <- matrix(0,nrow=0,ncol=2)
      colnames(inds) <- c("inpind","i")
      vals <- matrix(0,nrow=0,ncol=0)
      res <- list(inpinds=inds,parvals=vals)
    }

    res.trait[[trait.names[i]]] <- res

  }

  return(res.trait)
}


# get values of reach-depedent environmental conditions:
# ------------------------------------------------------

# The function searches for time-dependent inputs or constant values
# of parameters with names given without reach identifiers in the
# character vector "par.names".
# Reach identifiers are added to the search string when searching for
# input and parameter values for state variables.
# The function returns a list with a 3-column matrix of input and matrix
# indices in its columns for the parameter names that occur as input
# elements (values of those components will need updating as a function
# of time)and a matrix of numerical values of the occrrence in the
# parameter vector or in the defaults vector. This matrix contains the
# values of all parameters (columns) for all state variables (rows).
# The function issues warnings for parameters that have been declared
# as required and do neither occur in the input, nor in the parameter
# vector and do not have default values.

get.inpind.parval.envcond.reach <- function(par.names,y.names,par,inp=NA,
                                            required=NA,defaults=NA)
{
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  inds.vals <- get.inpind.parval.reach(
    par.names, y.names, par, inp = inp, defaults = defaults)

  check.inpind.parval.required(inds.vals, required = required)

  return(inds.vals)
}


# get values of habitat-(and reach-)depedent environmental conditions:
# --------------------------------------------------------------------

# The function searches for time-dependent inputs or constant values
# of parameters with names given without reach and habitat identifiers
# in the character vector "par.names".
# Reach and habitat identifiers are added to the search string when
# searching for input and parameter values for state variables.
# The function returns a list with a 3-column matrix of input and matrix
# indices in its columns for the parameter names that occur as input
# elements (values of those components will need updating as a function
# of time) and a matrix of numerical values of the occrrence in the
# parameter vector or in the defaults vector. This matrix contains the
# values of all parameters (columns) for all state variables (rows).
# The function issues warnings for parameters that have been declared
# as required and do neither occur in the input, nor in the parameter
# vector and do not have default values.

get.inpind.parval.envcond.habitat <- function(par.names,y.names,par,inp=NA,
                                              required=NA,defaults=NA)
{
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  inds.vals <- get.inpind.parval.reach.habitat(
    par.names, y.names, par, inp = inp, defaults = defaults)

  check.inpind.parval.required(inds.vals, required = required)

  # normalize parameter fA:
  if ( !is.na(match("fA",par.names)) )
  {
    inds.vals$parvals[,"fA"] <- calc.fA.norm(inds.vals$parvals[,"fA"], y.names$ind.fA)
  }

  return(inds.vals)
}

# get values of habitat-(and reach-)depedent environmental conditions
# for a group of environmental conditions:
# ----------------------------------------

# Group.names can be specified to search for habitat specific environmental conditions
# that belong to the same group, e.g. areal fraction of different habitat types
# It automatically derives the parameter names that belong to the group and calls the
# function "get.inpind.parval.envcond.habitat"

get.inpind.parval.envcond.habitat.group <- function(group.names,y.names,par,inp=NA,
                                                    required=NA,defaults=NA)
{

  # decode state variable names if the function was not already called
  # with the decoded names:

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  # search par.names for trait.names

  res.group <- list()

  for (i in 1:length(group.names))
  {
    par.names <- character(0)

    ind1 <- grep(group.names[i],names(par))

    if(length(ind1)>0)
    {
      ind2 <- which(group.names[i]==unlist(strsplit(names(par[ind1]),split="_")))

      par.names <- c(par.names,
                     paste(group.names[i],unique(unlist(strsplit(names(par[ind1]),
                                                                 split="_"))[ind2+1]),
                           sep="_") )
    } else
    {
      warning("no parameters for group ", group.names[i], " found."  )
    }


    if (length(par.names)>0)
    {
      res <- get.inpind.parval.envcond.habitat(par.names=par.names,y.names,par=par,inp=inp,
                                               required=required,defaults=defaults)

      # return list of matrices with input indices and values
      # of all parameters for all state variables:

    } else
    {
      inds <- matrix(0,nrow=0,ncol=3)
      colnames(inds) <- c("inpind","i","j")
      vals <- matrix(0,nrow=0,ncol=0)
      res <- list(inpinds=inds,parvals=vals)
    }

    res.group[[group.names[i]]] <- res

  }

  return(res.group)
}


# get initial conditions and inputs:
# ----------------------------------

# The function searches for time-dependent inputs or constant values
# of parameters with names given without taxon, reach and habitat
# identifiers in the character vector "par.names".
# Taxon, reach and habitat identifiers are added to the search string
# when searching for input and parameter values for state variables.

# The function returns a list with a 3-column matrix of input and matrix
# indices in its columns for the parameter names that occur as input
# elements (values of those components will need updating as a function
# of time) and a matrix of numerical values of the occurrence in the
# parameter vector or in the defaults vector. This matrix contains the
# values of all parameters (columns) for all state variables (rows).
# The function issues warnings for parameters that have been declared
# as required and do neither occur in the input, nor in the parameter
# vector and do not have default values.

get.inpind.parval.initcondinput <- function(par.names,y.names,par,inp=NA,
                                            required=NA,defaults=NA)
{
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  inds.vals <- get.inpind.parval.taxon.group.reach.habitat(
    par.names, y.names, par, inp = inp, defaults = defaults)

  check.inpind.parval.required(inds.vals, required = required)

  return(inds.vals)
}


# get taxa properties:
# --------------------

# The function searches for time-dependent inputs or constant values
# of parameters with names given without taxon and group identifiers
# in the character vector "par.names".
# Taxon and group identifiers are added to the search string
# when searching for input and parameter values for state variables.
# The function returns a list with a 2-column matrix of input and row
# indices in its columns for the parameter names that occur as input
# elements (values of those components will need updating as a function
# of time) and a matrix of numerical values of the occrrence in the
# parameter vector or in the defaults vector. This matrix contains the
# values of all parameters (columns) for all state variables (rows).
# The function issues warnings for parameters that have been declared
# as required and do neither occur in the input, nor in the parameter
# vector and do not have default values.

# The function accepts in the character vector "par.names" parameter names
# without taxon, reach and habitat specifiers and returns a list of two
# matrices with the indices of inputs and the values of all parameters,
# respectively (columns) for all state variables (rows).
# A warning is issued if a parameter that has been declared as required
# is not available for some state variables.

get.inpind.parval.taxaprop <- function(par.names,y.names,par,inp=NA,
                                       required=NA,defaults=NA)
{
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  inds.vals <- get.inpind.parval.taxon.group(
    par.names, y.names, par, inp = inp, defaults = defaults)

  check.inpind.parval.required(inds.vals, required = required)

  return(inds.vals)
}


# get taxa properties belonging to a specific trait:
# --------------------------------------------------

# The function accepts in the character vector "trait.names" trait names
# without taxon, reach and habitat specifiers and returns two matrices with
# the indices of inputs and the values of all parameters, respectively
# (columns) for all state variables (rows) for all traits.
# In contrast to get.parval.taxaprop we do not have to specify the name of
# the parameter but only the name of the trait the parameter belongs to.
# A warning is issued if a parameter that has been declared as required
# is not available for some state variables.

get.inpind.parval.taxaprop.traits <- function(trait.names,xval.pars=NA,y.names,par,inp=NA,
                                              required=NA,defaults=NA)
{
  # decode state variable names if the function was not already called
  # with the decoded names:

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  # function to split names and return the second half of it

  get2ndname <- function(fullnames)
  {
    shortnames <- rep(NA,length(fullnames))
    splitnames <- strsplit(fullnames,split="_")
    for(i in 1:length(fullnames))
    {
      shortnames[i] <- splitnames[[i]][2]
    }
    return(shortnames)
  }

  # search par.names for trait.names
  par.names <- NULL
  res.traits <- list()
  if(is.list(xval.pars))
  {
    if(length(trait.names)!=length(xval.pars)) warning("length of trait.names does not equal length of xval.pars")
  }

  for (i in 1:length(trait.names))
  {
    if(is.list(xval.pars))
    {
      if(substring(trait.names[i],1,5)!=substring(names(xval.pars[[i]]),1,5))
      {
        warning(paste("trait.names",trait.names[i],"is not fitting to xval.names",names(xval.pars[[i]])))
      }

      if (length(xval.pars[[i]][[1]]$parvals)>0)
      {
        if(is.matrix(xval.pars[[i]][[1]]$parvals))
        {
          x.val.names <- get2ndname(colnames(xval.pars[[i]][[1]]$parvals))
        }  else {
          x.val.names <- get2ndname(names(xval.pars[[i]][[1]]$parvals))
        }
        par.names <- paste(trait.names[i],x.val.names,sep="_")
      }
    } else {

      # this else is needed to extract parameters like feedingtypes in pre-prosessiong routines (e.g. for stoichiometry)

      ind1 <- grep(paste(trait.names[i],"_",sep=""),names(par))
      if(length(ind1)>0)
      {
        ind2 <- which(trait.names[i]==unlist(strsplit(names(par[ind1]),split="_")))

        par.names <- c(par.names,
                       paste(trait.names[i],unique(unlist(strsplit(names(par[ind1]),
                                                                   split="_"))[ind2+1]),
                             sep="_") )
      } else
      {
        warning("no parameters for trait ", trait.names[i], " found."  )
      }
    }
    if (length(par.names)>0)
    {
      res <- get.inpind.parval.taxaprop(par.names=par.names,y.names,par=par,inp=inp,
                                        required=required,defaults=defaults)

      # return list of matrices with input indices and values
      # of all parameters for all state variables:

    } else
    {
      inds <- matrix(0,nrow=0,ncol=3)
      colnames(inds) <- c("inpind","i","j")
      vals <- matrix(0,nrow=0,ncol=0)
      res <- list(inpinds=inds,parvals=vals)
    }
    res.traits[[trait.names[i]]] <- res
  }

  return(res.traits)
}


# function to get a structured representation of the system definition
# --------------------------------------------------------------------

#' Get system definition of the streambugs ODE model
#'
#' Get a structured representation of the streambugs ODE system definition.
#'
#' @inheritParams run.streambugs
#'
#' @return List with definition of the model including input state variables,
#'    parameters, and inputs, as well as the derived model structure including
#'    global parameters, environmental conditions of reaches and habitats,
#'    initial conditions, taxa properties, stoichiometric coefficients of all
#'    processes, and process definitions for each state variable.
#'
#' @examples
#' m <- streambugs.example.model.toy()
#' sys.def <- streambugs.get.sys.def(y.names=m$y.names, par=m$par, inp=m$inp)
#' # Get inital conditions for all state variables
#' sys.def$par.initcond$parvals
#' # Get stoichiometric coefficients of consumption within the food web
#' sys.def$par.stoich.web$Cons
#' # Get stoichiometric coefficients of death process for each taxon
#' # (transforming them into a dead organic matter "POM")
#' sys.def$par.stoich.taxon$Death
#'
#' @export
streambugs.get.sys.def <- function(y.names,par,inp=NA)
{
  # decode state variable names if the function was not already called
  # with the decoded names:

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  # get global parameters:
  # ----------------------

  par.global <- get.inpind.parval.global(
    par.names = c("kBoltzmann",
                  "M0",
                  "ftempmax_intercept",
                  "ftempmax_curv",
                  "fcurrent_intercept",
                  "fcurrent_curv",
                  "fsapro_intercept",
                  "fsapro_curv",
                  "forgmicropoll_intercept",
                  "forgmicropoll_curv",
                  "fmicrohab_intercept",
                  "fmicrohab_curv"
    ),
    par       = par,
    inp       = inp,
    required  = c("kBoltzmann",
                  "M0",
                  "ftempmax_intercept",
                  "ftempmax_curv",
                  "fcurrent_intercept",
                  "fcurrent_curv",
                  "fsapro_intercept",
                  "fsapro_curv",
                  "forgmicropoll_intercept",
                  "forgmicropoll_curv",
                  "fmicrohab_intercept",
                  "fmicrohab_curv"

    ),
    defaults  = c(kBoltzmann              = 8.61734e-005,
                  M0                      = 1,
                  ftempmax_intercept      = 0,
                  ftempmax_curv           = 0,
                  fcurrent_intercept      = 0,
                  fcurrent_curv           = 0,
                  fsapro_intercept        = 0,
                  fsapro_curv             = 0,
                  forgmicropoll_intercept = 0,
                  forgmicropoll_curv      = 0,
                  fmicrohab_intercept     = 0,
                  fmicrohab_curv          = 0
    ))

  par.global.envtraits <- get.inpind.parval.global.envcondtraits(

    trait.names   = c("tempmaxKval",
                      "currentmsval",
                      "saprowqclassval",
                      "orgmicropollTUval"
    ),
    par=par,
    inp=inp,
    required=NA,
    defaults=NA)


  # get evironmental conditions:
  # ----------------------------

  par.envcond.reach   <- get.inpind.parval.envcond.reach(
    par.names = c("w","L"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    required  = c("w","L"))

  par.envcond.habitat <- get.inpind.parval.envcond.habitat(
    par.names = c("T",
                  "I0",
                  "fshade",
                  "CP",
                  "CN",
                  "DSusPOM",
                  "tau",
                  "taucrit",
                  "tempmaxK",
                  "currentms",
                  "orgmicropollTU",
                  "saprowqclass",
                  "fA",
                  "DFish"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    required  = c("T"),
    defaults  = c(fA=1,DFish=0))


  par.envcond.habitat.group <- get.inpind.parval.envcond.habitat.group(
    group.names = c("microhabaf"),
    y.names=y.names,
    par         = par,
    inp         = inp)

  # get initial conditions and inputs:
  # ----------------------------------

  par.initcondinput <- get.inpind.parval.initcondinput(
    par.names = c("Dini", "Input"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    defaults  = c(Dini=0, Input=0))

  par.initcond <- get.single.par.inpind.parval(par.initcondinput, "Dini")
  par.input <- get.single.par.inpind.parval(par.initcondinput, "Input")

  # get taxa properties:
  # --------------------

  par.taxaprop.direct <- get.inpind.parval.taxaprop(
    par.names = c("M",
                  "Ea",
                  "b",
                  "i0",
                  "EC",
                  "fbasaltax"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    defaults  = c(b=0.75,
                  fbasaltax=1))


  # if trait.names are extended, the xval.pars have to be extended as well!

  trait.names = c("saprotolval",
                  "orgmicropolltolval",
                  "currenttolval",
                  "tempmaxtolval",
                  "microhabtolval")

  # objects which contain the names of the x-values for interpolation that are used to
  # derive the trait-tolerance values of invertebrates in the correct order

  xval.pars <- list( par.global.envtraits["saprowqclassval"],
                     par.global.envtraits["orgmicropollTUval"],
                     par.global.envtraits["currentmsval"],
                     par.global.envtraits["tempmaxKval"],
                     par.envcond.habitat.group["microhabaf"])

  par.taxaprop.traits <- get.inpind.parval.taxaprop.traits(
    trait.names = trait.names,
    xval.pars   = xval.pars,
    y.names     = y.names,
    par         = par,
    inp         = inp)


  # get process stoichiometries:
  # ----------------------------

  # taxon-based processes:

  kinparnames.taxon <- list(Miner    = c("kminer"),
                            Drift    = c("cdet"),
                            Resp     = c("fresp"),
                            Death    = c("fdeath"),
                            Prod     = c("fprod",
                                         "fgrotax",
                                         "hdens",
                                         "KI",
                                         "KP",
                                         "KN"))
  par.stoich.taxon <- get.par.stoich.list(
    par, names(kinparnames.taxon), get.par.stoich.taxon)

  # food web processes:
  # -------------------

  kinparnames.web <- list(Cons     = list(taxon1 = c("fcons",
                                                     "fgrotax",
                                                     "hdens",
                                                     "Kfood",
                                                     "q"),
                                          taxa2  = c("Pref")),
                          FishPred = list(taxon1 = c("cfish",
                                                     "Kfood",
                                                     "q"),
                                          taxa2  = c("Pref")))
  par.stoich.web <- get.par.stoich.list(
    par, names(kinparnames.web), get.par.stoich.web)

  # get kinetic parameters of processes:
  # ------------------------------------
  par.proc.taxon <- get.par.proc.taxon(
    y.names, par, inp, kinparnames.taxon, par.stoich.taxon)
  par.proc.web <- get.par.proc.web(
    y.names, par, inp, kinparnames.web, par.stoich.web)

  # return list of system definition:
  # ---------------------------------

  return(list(y.names                   = y.names,
              par                       = par,
              inp                       = inp,
              par.global                = par.global,
              par.global.envtraits      = par.global.envtraits,
              par.envcond.reach         = par.envcond.reach,
              par.envcond.habitat       = par.envcond.habitat,
              par.envcond.habitat.group = par.envcond.habitat.group,
              par.initcond              = par.initcond,
              par.input                 = par.input,
              par.taxaprop.direct       = par.taxaprop.direct,
              par.taxaprop.traits       = par.taxaprop.traits,
              par.stoich.taxon          = par.stoich.taxon,
              par.stoich.web            = par.stoich.web,
              par.proc.taxon            = par.proc.taxon,
              par.proc.web              = par.proc.web))
}

get.single.par.inpind.parval <- function(inpind.parval, par.name) {
  par.idx <- which(colnames(inpind.parval$parvals)==par.name)
  par.inpind.parval <- list(
    inpinds = inpind.parval$inpinds[inpind.parval$inpinds[,"j"]==par.idx,, drop=FALSE],
    parvals = inpind.parval$parvals[, par.idx, drop=FALSE])
  par.inpind.parval$inpinds[, "j"] = 1
  par.inpind.parval
}

# function to write the system definition to a file:
# --------------------------------------------------

#' Write system definition of the streambugs ODE model
#'
#' Write system definition of the streambugs ODE model into a human-readable
#' text file.
#'
#' @param sys.def system definition generated by the function
#'    \code{\link{streambugs.get.sys.def}}
#' @param file file name
#'
#' @examples
#' m <- streambugs.example.model.toy()
#' sys.def <- streambugs.get.sys.def(y.names=m$y.names, par=m$par, inp=m$inp)
#' file.name <- tempfile(m$name, fileext=".dat")
#' streambugs.write.sys.def(sys.def, file.name)
#' file.show(file.name, delete.file=TRUE)
#'
#' @export
streambugs.write.sys.def <- function(sys.def,file=NA)
{
  # local function to write typical data structures:
  # ------------------------------------------------

  write.pars <- function(pars,inp.names,main=NA)
  {
    # write header:

    if ( !is.na(main) ) cat(main,"\n",sep="")

    if ( is.vector(pars$parvals) )
    {
      # define local variables:

      vals <- pars$parvals
      inds <- pars$inpinds

      # if available. replace values by inputs:

      if ( nrow(inds) > 0 )
      {
        for ( i in 1:nrow(inds) )
        {
          vals[inds[i,2]] <- paste("input:",inp.names[inds[i,1]])
        }
      }

      # write values and inputs:

      if ( length(vals) > 0 )
      {
        for ( i in 1:length(vals) )
        {
          cat(names(vals)[i],"\t",vals[i],"\n",sep="")
        }
      }
    }
    else
    {
      if ( is.matrix(pars$parvals) )
      {
        # define local variables:

        vals <- pars$parvals
        inds <- pars$inpinds

        # if available, replace values by inputs:

        if ( nrow(inds) > 0 )
        {
          for ( i in 1:nrow(inds) )
          {
            vals[inds[i,2],inds[i,3]] <- paste("input:",inp.names[inds[i,1]])
          }
        }

        # write values and inputs:

        if ( ncol(vals) > 0 )
        {
          for ( j in 1:ncol(vals) ) cat("\t",colnames(vals)[j],sep="")
          cat("\n")
          for ( i in 1:nrow(vals) )
          {
            cat(rownames(vals)[i])
            for ( j in 1:ncol(vals) )
            {
              cat("\t",vals[i,j],sep="")
            }
            cat("\n")
          }
        }
      }
      else
      {
        cat("*** parameter format not yet implemented ***\n")
      }
    }
    cat("\n")
  }

  if ( !is.na(file) ) sink(file=file) # redirect output to file

  # header:
  # -------

  cat("=======================================================\n")
  cat("Streambugs Model Definition\n")
  cat("=======================================================\n")
  cat("\n")
  cat("-------------------------------------------------------\n")
  cat("Input to the Model\n")
  cat("-------------------------------------------------------\n")
  cat("\n")

  # write state variable names:
  # ---------------------------

  y.names <- sys.def$y.names
  cat("State Variables (",length(y.names$y.names),"):\n",sep="")
  for ( i in 1:length(y.names$y.names)) cat(y.names$y.names[i],"\n",sep="")
  cat("\n")

  # write inputs:
  # -------------

  inp <- sys.def$inp
  if ( is.list(inp) )
  {
    cat("Input Definitions (",length(inp),"):\n",sep="")
    for ( i in 1:length(inp) )
    {
      cat(names(inp)[i],"\t",nrow(inp[[i]])," data pairs\n",sep="")
    }
    cat("\n")
  }

  # write parameters:
  # -----------------

  par <- sys.def$par
  cat("Parameters (",length(par),"):\n",sep="")
  for ( i in 1:length(par) )
  {
    cat(names(par)[i],"\t",par[i],"\n",sep="")
  }
  cat("\n")

  cat("-------------------------------------------------------\n")
  cat("Derived Model Structure\n")
  cat("-------------------------------------------------------\n")
  cat("\n")

  # write global parameters:
  # ------------------------

  write.pars(sys.def$par.global,names(sys.def$inp),"Global Parameters:")

  # write global parameters related to environmental conditions and traits:
  # -----------------------------------------------------------------------

  for (i in 1:length(sys.def$par.global.envtraits))
  {
    write.pars(sys.def$par.global.envtraits[[i]],names(sys.def$inp),
               paste("Global Parameters related to environmental conditions and trait: ",names(sys.def$par.global.envtraits )[i]))
  }

  # write reach-dependent environmental conditions:
  # ----------------------------------------------

  write.pars(sys.def$par.envcond.reach,names(sys.def$inp),"Reach-Dependent Environmental Conditions:")

  # write habitat-dependent environmental conditions:
  # ------------------------------------------------

  write.pars(sys.def$par.envcond.habitat,names(sys.def$inp),"Habitat-Dependent Environmental Conditions:")

  # write habitat-dependent environmental conditions belonging to a group:
  # ----------------------------------------------------------------------

  for (i in 1:length(sys.def$par.envcond.habitat.group))
  {
    write.pars(sys.def$par.envcond.habitat.group[[i]],names(sys.def$inp),
               paste("Habitat-Dependent Environmental Conditions belonging to group:", names(sys.def$par.envcond.habitat.group)[i]))
  }

  # write initial conditons:
  # ------------------------

  write.pars(sys.def$par.initcond,names(sys.def$inp),"Initial Conditions:")

  # write input:
  # ------------

  write.pars(sys.def$par.input,names(sys.def$inp),"Input:")

  # write directly defined taxa properties:
  # ---------------------------------------

  write.pars(sys.def$par.taxaprop.direct,names(sys.def$inp),"Directly Defined Taxa Properties:")

  # write trait-derived taxa properties:
  # ------------------------------------

  for(i in 1: length(sys.def$par.taxaprop.trait))
  {
    write.pars(sys.def$par.taxaprop.trait[[i]],names(sys.def$inp),paste("Trait-Derived Taxa Properties:",names(sys.def$par.taxaprop.trait)[i]))

  }

  # write process stoichiometries of taxon-based processes:
  # -------------------------------------------------------

  par.stoich.taxon <- sys.def$par.stoich.taxon
  for( i in 1:length(par.stoich.taxon) )  # process i
  {
    if( length(par.stoich.taxon[[i]]) > 0 )
    {
      cat("Stoichiometry of Process ",names(par.stoich.taxon)[i],":\n",sep="")
      for ( j in 1:length(par.stoich.taxon[[i]]) )  # process i for organism j
      {
        for( k in 1:length(par.stoich.taxon[[i]][[j]]) )  # stoich. coeff. for subst k
        {
          cat("\t",names(par.stoich.taxon[[i]][[j]])[k],sep="")
        }
        cat("\n")
        cat("  ",names(par.stoich.taxon[[i]])[j],":",sep="")
        for( k in 1:length(par.stoich.taxon[[i]][[j]]) )
        {
          cat("\t",par.stoich.taxon[[i]][[j]][k],sep="")
        }
        cat("\n")
      }
      cat("\n")
    }
  }

  # write process stoichiometry of food web processes:
  # --------------------------------------------------

  par.stoich.web <- sys.def$par.stoich.web
  for( i in 1:length(par.stoich.web) )  # process i
  {
    if( length(par.stoich.web[[i]]) > 0 )
    {
      cat("Stoichiometry of Process ",names(par.stoich.web)[i],":\n",sep="")
      for ( j in 1:length(par.stoich.web[[i]]) )  # process i
      {
        for ( l in 1:length(par.stoich.web[[i]][[j]]) )  # process i for organisms j,l
        {
          for( k in 1:length(par.stoich.web[[i]][[j]][[l]]) )  # stoich. coeff. for subst k
          {
            cat("\t",names(par.stoich.web[[i]][[j]][[l]])[k],sep="")
          }
          cat("\n")
          cat("  ",names(par.stoich.web[[i]])[j],"-",names(par.stoich.web[[i]][[j]])[l],":",sep="")
          for( k in 1:length(par.stoich.web[[i]][[j]][[l]]) )
          {
            cat("\t",par.stoich.web[[i]][[j]][[l]][k],sep="")
          }
          cat("\n")
        }
      }
    }
    cat("\n")
  }

  # write process definitions by state variables:
  # ---------------------------------------------

  par.proc.taxon <- sys.def$par.proc.taxon
  par.proc.web   <- sys.def$par.proc.web
  cat("Process Definitions by State Variables:\n")
  for ( i in 1:length(y.names$y.names) )
  {
    var <- y.names$y.names[i]
    cat("\nState Variable ",i,":\t",var,"\n",sep="")
    if ( length(par.proc.taxon[[var]]) > 0 )
    {
      for ( j in 1:length(par.proc.taxon[[var]]) )
      {
        # define local variables:

        vals   <- par.proc.taxon[[var]][[j]]$parvals
        inds   <- par.proc.taxon[[var]][[j]]$inpinds
        stoich <- par.proc.taxon[[var]][[j]]$stoich

        # if available. replace values by inputs:

        if ( nrow(inds) > 0 )
        {
          for ( k in 1:nrow(inds) )
          {
            vals[inds[k,2]] <- paste("input:",names(inp)[inds[k,1]])
          }
        }

        # write values and inputs:

        cat("  ",names(par.proc.taxon[[var]])[j],":\n",sep="")
        cat("      Parameter:")
        for ( k in 1:length(vals) ) cat("\t",names(vals)[k],sep="")
        cat("\n")
        cat("      Value:")
        for ( k in 1:length(vals) ) cat("\t",vals[k],sep="")
        cat("\n")
        cat("      Taxon:")
        for ( k in 1:ncol(stoich) ) cat("\t",colnames(stoich)[k],sep="")
        cat("\n")
        cat("      StateVar:")
        for ( k in 1:ncol(stoich) ) cat("\t",stoich[1,k],sep="")
        cat("\n")
        cat("      StoichCoeff:")
        for ( k in 1:ncol(stoich) ) cat("\t",stoich[2,k],sep="")
        cat("\n")
      }
    }
    if ( length(par.proc.web[[var]]) > 0 )
    {
      for ( j in 1:length(par.proc.web[[var]]) )
      {
        # general process parameters:

        # define local variables:

        vals <- par.proc.web[[var]][[j]]$parvals
        inds <- par.proc.web[[var]][[j]]$inpinds

        # if available. replace values by inputs:

        if ( nrow(inds) > 0 )
        {
          for ( k in 1:nrow(inds) )
          {
            vals[inds[k,2]] <- paste("input:",names(inp)[inds[k,1]])
          }
        }

        # write values and inputs:

        cat("  ",names(par.proc.web[[var]])[j],":\n",sep="")
        cat("      Parameter:")
        for ( k in 1:length(vals) ) cat("\t",names(vals)[k],sep="")
        cat("\n")
        cat("      Value:")
        for ( k in 1:length(vals) ) cat("\t",vals[k],sep="")
        cat("\n")

        # food-specific process parameters:

        for ( k in 1:length(par.proc.web[[var]][[j]]$taxa2) )
        {
          # define local variables:

          vals   <- par.proc.web[[var]][[j]][["taxa2"]][[k]]$parvals
          inds   <- par.proc.web[[var]][[j]][["taxa2"]][[k]]$inpinds
          stoich <- par.proc.web[[var]][[j]][["taxa2"]][[k]]$stoich

          # if available. replace values by inputs:

          if ( nrow(inds) > 0 )
          {
            for ( k in 1:nrow(inds) )
            {
              vals[inds[k,2]] <- paste("input:",names(inp)[inds[k,1]])
            }
          }

          # write values and inputs:

          cat("    ",names(par.proc.web[[var]][[j]][["taxa2"]])[k],":\n",sep="")
          cat("      Parameter:")
          for ( k in 1:length(vals) ) cat("\t",names(vals)[k],sep="")
          cat("\n")
          cat("      Value:")
          for ( k in 1:length(vals) ) cat("\t",vals[k],sep="")
          cat("\n")
          cat("      Taxon:")
          for ( k in 1:ncol(stoich) ) cat("\t",colnames(stoich)[k],sep="")
          cat("\n")
          cat("      StateVar:")
          for ( k in 1:ncol(stoich) ) cat("\t",stoich[1,k],sep="")
          cat("\n")
          cat("      StoichCoeff:")
          for ( k in 1:ncol(stoich) ) cat("\t",stoich[2,k],sep="")
          cat("\n")
        }
      }
    }
  }

  if ( !is.na(file) ) sink() # terminate redirection of output
}


# function to update parameter values for a given time, t, by interpolation of inputs:
# ------------------------------------------------------------------------------------

# interpolate inputs:

interpolate.inputs <- function(inp,t)
{
  if ( !is.list(inp) ) return(NA)
  if ( length(inp) == 0 ) return(NA)
  inpvals <- rep(NA,length(inp))
  for ( i in 1:length(inp) )
  {
    inpvals[i] <- approx(x=inp[[i]][,1],y=inp[[i]][,2],xout=t,rule=2)$y
  }
  names(inpvals) <- names(inp)
  return(inpvals)
}

# update reach-dependent environmental conditions:

streambugs.update.envcond.reach <- function(par.envcond.reach,inpvals)
{
  if ( nrow(par.envcond.reach$inpinds) > 0 )
  {
    for ( i in 1:nrow(par.envcond.reach$inpinds) )
    {
      par.envcond.reach$parvals[par.envcond.reach$inpinds[i,2],
                                par.envcond.reach$inpinds[i,3]] <-
        inpvals[par.envcond.reach$inpinds[i,1]]
    }
  }
  return(par.envcond.reach$parvals)
}

# update habitat (and reach) - dependent environmental conditions:

streambugs.update.envcond.habitat <- function(par.envcond.habitat,inpvals,ind.fA)
{
  if ( nrow(par.envcond.habitat$inpinds) > 0 )
  {
    for ( i in 1:nrow(par.envcond.habitat$inpinds) )
    {
      par.envcond.habitat$parvals[par.envcond.habitat$inpinds[i,2],
                                  par.envcond.habitat$inpinds[i,3]] <-
        inpvals[par.envcond.habitat$inpinds[i,1]]
    }
  }
  # normalize parameter fA:
  if ( !is.na(match("fA",colnames(par.envcond.habitat$parvals))) )
  {
    par.envcond.habitat$parvals[,"fA"] <-
      calc.fA.norm(par.envcond.habitat$parvals[,"fA"],ind.fA)
  }
  return(par.envcond.habitat$parvals)
}

# update all:

streambugs.updatepar <- function(sys.def,t)
{
  if ( is.list(sys.def$inp) )
  {
    # interpolate inputs:
    # -------------------

    inpvals <- interpolate.inputs(sys.def$inp,t)

    # substitute interpolated inputs into parameter arrays where necessary:
    # ---------------------------------------------------------------------

    # update global parameters:

    if ( nrow(sys.def$par.global$inpinds) > 0 )
    {
      for ( i in 1:nrow(sys.def$par.global$inpinds) )
      {
        sys.def$par.global$parvals[sys.def$par.global$inpinds[i,2]] <-
          inpvals[sys.def$par.global$inpinds[i,1]]
      }
    }

    # update reach-dependent environmental conditions:

    sys.def$par.envcond.reach$parvals <-
      streambugs.update.envcond.reach(sys.def$par.envcond.reach,inpvals)

    # update habitat (and reach) - dependent environmental conditions:

    sys.def$par.envcond.habitat$parvals <-
      streambugs.update.envcond.habitat(sys.def$par.envcond.habitat,inpvals,sys.def$y.names$ind.fA)

    # update microhabitat - conditions:

    sys.def$par.envcond.habitat.group$microhabaf$parvals <-
      streambugs.update.envcond.habitat(sys.def$par.envcond.habitat.group$microhabaf,inpvals,sys.def$y.names$ind.fA)

    # update initial conditions:

    if ( nrow(sys.def$par.initcond$inpinds) > 0 )
    {
      for ( i in 1:nrow(sys.def$par.initcond$inpinds) )
      {
        sys.def$par.initcond$parvals[sys.def$par.initcond$inpinds[i,2],
                                     sys.def$par.initcond$inpinds[i,3]] <-
          inpvals[sys.def$par.initcond$inpinds[i,1]]
      }
    }

    # update inputs:

    if ( nrow(sys.def$par.input$inpinds) > 0 )
    {
      for ( i in 1:nrow(sys.def$par.input$inpinds) )
      {
        sys.def$par.input$parvals[sys.def$par.input$inpinds[i,2],
                                  sys.def$par.input$inpinds[i,3]] <-
          inpvals[sys.def$par.input$inpinds[i,1]]
      }
    }

    # update directly defined taxa properties:

    if ( nrow(sys.def$par.taxaprop.direct$inpinds) > 0 )
    {
      for ( i in 1:nrow(sys.def$par.taxaprop.direct$inpinds) )
      {
        sys.def$par.taxaprop.direct$parvals[sys.def$par.taxaprop.direct$inpinds[i,2],
                                            sys.def$par.taxaprop.direct$inpinds[i,3]] <-
          inpvals[sys.def$par.taxaprop.direct$inpinds[i,1]]
      }
    }

    # update trait-dependent taxa properties:
    for(j in 1:length(sys.def$par.taxaprop.traits))
    {
      if ( nrow(sys.def$par.taxaprop.traits[[j]]$inpinds) > 0 )
      {
        for ( i in 1:nrow(sys.def$par.taxaprop.traits[[j]]$inpinds) )
        {
          sys.def$par.taxaprop.traits[[j]]$parvals[sys.def$par.taxaprop.traits[[j]]$inpinds[i,2],  sys.def$par.taxaprop.traits[[j]]$inpinds[i,3]] <-
            inpvals[sys.def$par.taxaprop.traits[[j]]$inpinds[i,1]]
        }
      }
    }


    # par.stoich.taxon and par.stoich.web do not need updating
    # as they are are not allowed to be time-dependent

    # update kinetic parameters of taxon-based processes:

    for ( i in 1:length(sys.def$par.kin.taxon) )
    {
      if ( length(sys.def$par.kin.taxon[[i]]) > 0 )
      {
        for ( j in 1:length(sys.def$par.kin.taxon[[i]]) )
        {
          if ( nrow(sys.def$par.kin.taxon[[i]][[j]]$inpinds) > 0 )
          {
            for ( k in 1:nrow(sys.def$par.kin.taxon[[i]][[j]]$inpinds) )
            {
              sys.def$par.kin.taxon[[i]][[j]]$parvals[sys.def$par.kin.taxon$inpinds[k,2]] <-
                inpvals[sys.def$par.kin.taxon[[i]][[j]]$inpinds[k,1]]
            }
          }
        }
      }
    }

    # update kinetic parameters of food web processes:

    for ( i in 1:length(sys.def$par.kin.web) )
    {
      if ( length(sys.def$par.kin.web[[i]]) > 0 )
      {
        for ( j in 1:length(sys.def$par.kin.web[[i]]) )
        {
          if ( nrow(sys.def$par.kin.web[[i]][[j]]$inpinds) > 0 )
          {
            for ( k in 1:nrow(sys.def$par.kin.web[[i]][[j]]$inpinds) )
            {
              sys.def$par.kin.web[[i]][[j]]$parvals[sys.def$par.kin.web$inpinds[k,2]] <-
                inpvals[sys.def$par.kin.web[[i]][[j]]$inpinds[k,1]]
            }
            for ( k in 1:length(sys.def$par.kin.web[[i]][[j]]$taxa2) )
            {
              if ( nrow(sys.def$par.kin.web[[i]][[j]]$taxa2$inpinds) > 0 )
              {
                for ( l in 1:nrow(sys.def$par.kin.web[[i]][[j]]$taxa2$inpinds) )
                {
                  sys.def$par.kin.web[[i]][[j]]$taxa2[[k]]$parvals[sys.def$par.kin.web$taxa2[[k]]$inpinds[l,2],
                                                                   sys.def$par.kin.web$taxa2[[k]]$inpinds[l,3]] <-
                    inpvals[sys.def$par.kin.web[[i]][[j]]$taxa2[[k]]$inpinds[l,1]]
                }
              }
            }
          }
        }
      }
    }
  }

  return(sys.def)
}


# function to normalize habitat area fractions:
# ---------------------------------------------

calc.fA.norm <- function(fA,ind.fA)
{
  for ( i in 1:length(ind.fA) )
  {
    ind.reach  <- ind.fA[[i]][["ind.reach"]]
    ind.1sthab <- ind.fA[[i]][["ind.1sthab"]]
    fact <- 1/sum(fA[ind.reach][ind.1sthab])
    fA[ind.reach] <- fact*fA[ind.reach]
  }
  return(fA)
}

# function to calculate additional output (rates, limiting factors)
# -----------------------------------------------------------------

calculate.additional.output <- function(res,par,inp,file.add=NA,tout.add=NA)
{
  y.names <- colnames(res)[-1]
  tout <- res[,1]

  if( sum(!is.na(tout.add)) > 0 )
  {
    if( length(intersect(tout,tout.add)) < length(tout.add) )
    {
      stop("tout.add: ",tout.add," not part of tout")
    } else  tout <- tout.add
  }

  sys.def <- streambugs.get.sys.def(y.names=y.names,par=par,inp=inp)


  for (i in 1:length(tout))
  {
    sys.def.i <- streambugs.updatepar(sys.def,t=tout[i])

    res.add.i <- unlist(rhs.streambugs(t=tout[i],y=res[i,-1],parms=sys.def.i,
                                       add.out=TRUE))

    if(i==1) res.add <- res.add.i else  res.add <- rbind(res.add,res.add.i,deparse.level=0)
  }

  if(length(tout)==1) res.add <- as.matrix(t(res.add))

  res.add[,1] <- tout

  colnames(res.add)[1] <- colnames(res)[1]

  if(!is.na(file.add)) write.table(res.add,file.add,col.names=T,row.names=F,sep="\t")

  return(res.add)

}


exp.transform <- function(x,intercept=0,curv=0)
{
  #!if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left
  #!if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right

  if(curv == 0)
  {
    y = intercept-(intercept-1)*x
  } else {
    y = intercept - (intercept -1) * (1 - exp(-curv * x)) / (1-exp(-curv))
  }
  return(y)
}


#' Count feeding links between taxa in streambugs ODE.
#'
#' Count number of global "Cons" stoichiometric interactions (feeding links) between
#' streambugs ODE state variables (taxa) in all habitats.
#'
#' @param y.names same as  \code{y.names} in \code{\link{run.streambugs}}
#' @param par same as  \code{par} in \code{\link{run.streambugs}}
#'
#' @return Integer number of feeding links, which is a number of feeding links in a
#'    single habitat times number of habitats and number of reaches.
#'
#' @examples
#' m <- streambugs.example.model.toy()
#' count.feeding.links(m$y.names, m$par)
#'
#' # feeding links count does not change w/ number of habitats (nor reaches)
#' m <- streambugs.example.model.toy(n.Habitats=10)
#' count.feeding.links(m$y.names, m$par)
#'
#' @export
count.feeding.links <- function(y.names, par) {

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  stoich.Cons  <- get.par.stoich.web("Cons" ,par) #xxx

  foods <-  list()
  for ( i in 1:length(y.names$taxa) ){
    foods[[y.names$taxa[i]]] <- names(stoich.Cons[[y.names$taxa[i]]])
  }
  n.links <- length(unlist(foods))

  # simpler alternative that works only if par is consistent with states
  # n.links <- sum(sapply(stoich.Cons, length))

  n.habitats <- length(y.names$habitats)
  n.reaches <- length(y.names$reaches)

  return(n.reaches*n.habitats*n.links)
}
