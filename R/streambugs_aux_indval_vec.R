################################################################
#
# streambugs 1.2
# ==============
#
# -------------------
# auxiliary functions
# -------------------
#
# creation:      08.04.2020
# modifications: 08.04.2020
#
################################################################

# sparsify inputs logical matrix to a matrix with true indices in rows
# --------------------------------------------------------------------
sparsify.inpmatch <- function(inds) {
  inds.ij <- which(!is.na(inds), arr.ind = TRUE)
  # inpind is a position in a 1D column-wise representation of inds matrix, i.e.
  #  for c(i, j) row in inds.ij inpind is (i + n * (j - 1)), where n is #rows in inds
  inds <- cbind(inds[inds.ij[,"row"] + nrow(inds) * (inds.ij[,"col"] - 1)], inds.ij)
  colnames(inds) <- c("inpind","i","j")
  inds
}


# check parameter vector and input list for a vector of potential names:
# ----------------------------------------------------------------------

# auxiliary functions
.na.if.null <- function (x) ifelse(is.null(x), NA, x) # Beware: not a vector function
.first.not.na <- function(x) .na.if.null(Find(stats::complete.cases, x))

# Check the input list and the parameter vector for occurrence of a set
# of input or parameter names using all comibnations of these with provided parameter
# names prefixes (e.g. `"Reach1_Hab1_"``, `"Reach1_Hab2_"``, ..., `"Reach2_Hab1_"``,
# ...).
#
# If the search is not successful use the default value if provided.
#
# Prefixes may encode order of search using columns, with first row prefixes having
# priority. For instance if first column of paramter names prefixes contain
# `c("Reach1_Hab1", "Hab1", "Reach1", "")`` and parameters vector contains entries for
# both `"Hab1_w"` and `"w"`, then, when quering for parameter `"w"` value of `"Hab1_w"`
# will be returned as a first value.
#
# Return a list with the matching indices of the input and the values of the parameters.
#
get.inpind.parval.vec <- function(
  prefix.names.cols,
  par.names,
  y.names,
  par,
  inp,
  defaults,
  extend.results.intvec = NULL,
  if.inpmatch = FALSE) # instead of input indices return full logical array of input matches
{

  if (is.list(y.names))
    y.names <- y.names$y.names

  if (is.null(extend.results.intvec))
    extend.results.intvec <- 1:ncol(prefix.names.cols)

  # generate 3D array of ordered parameters lookup templates (e.g. Reach_Habitat_Par >
  # Habitat_Par > Reach_Par > Par) by lookup variables (e.g. Reach1_Hab1, Reach1_Hab2,
  # ..., ReachN_HabM) by unique parameter names (e.g. w, L)

  par.names.factor <- factor(par.names)
  par.names.unique <- levels(par.names.factor)
  idx.par.names.unique <- as.integer(par.names.factor)
  # we have: all(par.names.unique[idx.par.names.unique] == par.names)

  query.names.arr <- outer(prefix.names.cols, par.names.unique, paste0)

  # Query parameter values in perfix lookup order using plain indexing
  # Note: append default values to par for indexing - it doesn't matter if param is
  #       already defined in `par`  as first occurence is returned in case of duplicate
  #       vector element names, cf. `c(a=1, b=2, a=3)["a"] == 1`

  par.values.arr <- c(par,defaults)[query.names.arr]
  dim(par.values.arr) <- dim(query.names.arr)
  # Find first not `NA` value for each variables and parameters lookup order
  vals <- apply(par.values.arr, c(2,3), .first.not.na)

  # Lookup input parameter indices in prefix lookup order using `match`

  if ( is.list(inp) ) {
    inp.match.arr <- match(query.names.arr, names(inp))
    dim(inp.match.arr) <- dim(query.names.arr)
    # Find first not `NA` value for each variables and parameters lookup order
    inds <- apply(inp.match.arr, c(2,3), .first.not.na)
  } else {
    inds <- array(dim=dim(query.names.arr)[2:3])
  }

  # extend parameters vector

  vals <- vals[extend.results.intvec, idx.par.names.unique]
  if (is.null(dim(vals))) # 1st or 2nd dim was dropped
    dim(vals) <- c(length(y.names), length(par.names))
  colnames(vals) <- par.names
  rownames(vals) <- y.names

  # extend indices vector

  inds <- inds[extend.results.intvec, idx.par.names.unique]
  if (is.null(dim(inds))) # in case 1st or 2nd dim was dropped
    dim(inds) <- c(length(y.names), length(par.names))

  if (if.inpmatch) {
    # Note: if this would be put before if, then the sparsified rows would have names
    #       of variables given in column "i" (and could be looked-up this way)
    colnames(inds) <- par.names
    rownames(inds) <- y.names
  } else {
    inds <- sparsify.inpmatch(inds)
  }

  # Note: `inds` would be better renamed to `inp` to designate either match or indices
  #       array, not only the latter, same for inpinds in subsequent functions etc.
  return(list(inpinds=inds, parvals=vals))
}


# get values of reach-depedent parameters and inputs
# --------------------------------------------------

get.inpind.parval.reach <- function(par.names, y.names, par, inp = NA, defaults = NA)
{
  # get unique reaches to save runtime and query values only for these
  # use factor for easy mapping from vars (y.names) to unique reaches
  query.strvec <- y.names$y.reaches
  # pass unique levels to keep factor levels ordered as they apper in original vectors
  query.factor <- factor(query.strvec, levels=unique(query.strvec))
  query.intvec <- as.integer(query.factor)
  idx.query.y.names <- which(!duplicated(query.intvec))

  # define search order
  pref.y.reaches <- paste0(y.names$y.reaches[idx.query.y.names], "_")
  pref.names.cols <- rbind(
    pref.y.reaches, # 1: Reach_Par
    "")             # 2: Par

  vals.inds <- get.inpind.parval.vec(
    pref.names.cols, par.names, y.names, par, inp, defaults, query.intvec)

  return(vals.inds)
}


# get values of reach- and habitat-depedent parameters and inputs
# ---------------------------------------------------------------

get.inpind.parval.reach.habitat <- function(
  par.names, y.names, par, inp = NA, defaults = NA)
{
  # get unique (reach, habitat) pairs to save runtime and query values only for these
  # use factor for easy mapping from vars (y.names) to unique pairs
  query.strvec <- paste(y.names$y.reaches, y.names$y.habitats, sep="_")
  # pass unique levels to keep factor levels ordered as they apper in original vectors
  query.factor <- factor(query.strvec, levels=unique(query.strvec))
  query.intvec <- as.integer(query.factor)
  idx.query.y.names <- which(!duplicated(query.intvec))

  # define search order
  pref.y.habitats <- paste0(y.names$y.habitats[idx.query.y.names], "_")
  pref.y.reaches <- paste0(y.names$y.reaches[idx.query.y.names], "_")
  pref.names.cols <- rbind(
    paste0(pref.y.reaches, pref.y.habitats),  # 1: Reach_Habitat_Par
    pref.y.habitats,                          # 2: Habitat_Par
    pref.y.reaches,                           # 3: Reach_Par
    "")                                       # 4: Par

  vals.inds <- get.inpind.parval.vec(
    pref.names.cols, par.names, y.names, par, inp, defaults, query.intvec)

  return(vals.inds)
}


# get values of taxon-, and group-depedent parameters and inputs
# --------------------------------------------------------------

get.inpind.parval.taxon.group <- function(
  par.names, y.names, par, inp = NA, defaults = NA)
{
  # get unique (taxa, groups) pairs to save runtime and query values only for these
  # use factor for easy mapping from vars (y.names) to unique pairs
  query.strvec <- paste(y.names$y.taxa, y.names$y.groups, sep="_")
  # pass unique levels to keep factor levels ordered as they apper in original vectors
  query.factor <- factor(query.strvec, levels=unique(query.strvec))
  query.intvec <- as.integer(query.factor)
  idx.query.y.names <- which(!duplicated(query.intvec))

  # define search order
  pref.y.taxa <- paste0(y.names$y.taxa[idx.query.y.names], "_")
  pref.y.groups <- paste0(y.names$y.groups[idx.query.y.names], "_")
  pref.names.cols <- rbind(
    pref.y.taxa,    # 1: Taxon_Par
    pref.y.groups,  # 2: Group_Par
    "")             # 3: Par

  vals.inds <- get.inpind.parval.vec(
    pref.names.cols, par.names, y.names, par, inp, defaults, query.intvec)

  return(vals.inds)
}


# get values of taxon-, group-, reach- and habitat-depedent parameters and inputs
# -------------------------------------------------------------------------------

get.inpind.parval.taxon.group.reach.habitat <- function(
  par.names, y.names, par, inp = NA, defaults = NA)
{
  # get unique query tuples to save runtime and query values only for these
  # use factor for easy mapping from vars (y.names) to unique tuples
  query.strvec <- paste(
    y.names$y.taxa, y.names$y.groups, y.names$y.reaches, y.names$y.habitats, sep="_")
  # pass unique levels to keep factor levels ordered as they apper in original vectors
  query.factor <- factor(query.strvec, levels=unique(query.strvec))
  query.intvec <- as.integer(query.factor)
  idx.query.y.names <- which(!duplicated(query.intvec))

  # define search order
  pref.y.taxa <- paste0(y.names$y.taxa[idx.query.y.names], "_")
  pref.y.groups <- paste0(y.names$y.groups[idx.query.y.names], "_")
  pref.y.habitats <- paste0(y.names$y.habitats[idx.query.y.names], "_")
  pref.y.reaches <- paste0(y.names$y.reaches[idx.query.y.names], "_")
  pref.names.cols <- rbind(
    paste0(pref.y.taxa, pref.y.reaches, pref.y.habitats),    #  1: Taxon_Reach_Habitat_Par
    paste0(pref.y.groups, pref.y.reaches, pref.y.habitats),  #  2: Group_Reach_Habitat_Par
    paste0(pref.y.taxa, pref.y.habitats),                    #  3: Taxon_Habitat_Par
    paste0(pref.y.groups, pref.y.habitats),                  #  4: Group_Habitat_Par
    paste0(pref.y.taxa, pref.y.reaches),                     #  5: Taxon_Reach_Par
    paste0(pref.y.groups, pref.y.reaches),                   #  6: Group_Reach_Par
    pref.y.taxa,                                             #  7: Taxon_Par
    pref.y.groups,                                           #  8: Group_Par
    paste0(pref.y.reaches, pref.y.habitats),                 #  9: Reach_Habitat_Par
    pref.y.habitats,                                         # 10: Habitat_Par
    pref.y.reaches,                                          # 11: Reach_Par
    "")                                                      # 12: Par

  vals.inds <- get.inpind.parval.vec(
    pref.names.cols, par.names, y.names, par, inp, defaults, query.intvec)

  return(vals.inds)
}



# get values of pairs taxon-taxon-, group-taxon-, taxon-group-, group-group-, as well as
# reach- and habitat-depedent parameters and inputs
# --------------------------------------------------------------------------------------

get.inpind.parval.taxon.group.pairs.reach.habitat <- function(
  par.names, y.names, par, inp = NA, defaults = NA)
{
  n <- length(y.names$y.names)

  # search values are pairs of:
  # * taxon-taxon
  # * group-taxon
  # * taxon-group
  # * group-group
  # note: "->" var-var sperator as "_" may appear in var name
  pairs = list(
    y.taxa.taxa = as.vector(t(outer(y.names$y.taxa, y.names$y.taxa, paste, sep="_"))),
    y.groups.taxa = as.vector(t(outer(y.names$y.groups, y.names$y.taxa, paste, sep="_"))),
    y.taxa.groups = as.vector(t(outer(y.names$y.taxa, y.names$y.groups, paste, sep="_"))),
    y.groups.groups = as.vector(t(outer(y.names$y.groups, y.names$y.groups, paste, sep="_"))),
    y.taxa = rep(y.names$y.taxa, each=n),
    y.groups = rep(y.names$y.groups, each=n),
    y.habitats = rep(y.names$y.habitats, each=n),
    y.reaches = rep(y.names$y.reaches, each=n),
    y.names = as.vector(t(outer(y.names$y.names, y.names$y.names, paste, sep="->"))))

  # get unique query tuples to save runtime and query values only for these
  # use factor for easy mapping from var (y.names) pairs to unique tuples
  query.strvec <- paste(
    pairs$y.taxa.taxa, pairs$y.groups.taxa, pairs$y.taxa.groups, pairs$y.groups.groups,
    pairs$y.taxa, pairs$y.group, pairs$y.reaches, pairs$y.habitats, sep="_")
  # pass unique levels to keep factor levels ordered as they apper in original vectors
  query.factor <- factor(query.strvec, levels=unique(query.strvec))
  query.intvec <- as.integer(query.factor)
  idx.query.y.names <- which(!duplicated(query.intvec))

  # define search order
  pref <- lapply(pairs, function(y) paste0(y[idx.query.y.names], "_"))

  pref.names.cols <- rbind(
    paste0(pref$y.taxa.taxa, pref$y.reaches, pref$y.habitats),      #  1: Taxon1_Taxon2_Reach_Habitat_Par
    paste0(pref$y.groups.taxa, pref$y.reaches, pref$y.habitats),    #  2: Group1_Taxon2_Reach_Habitat_Par
    paste0(pref$y.taxa.groups, pref$y.reaches, pref$y.habitats),    #  3: Taxon1_Group2_Reach_Habitat_Par
    paste0(pref$y.groups.groups, pref$y.reaches, pref$y.habitats),  #  4: Group1_Group2_Reach_Habitat_Par
    paste0(pref$y.taxa.taxa, pref$y.habitats),                      #  5: Taxon1_Taxon2_Habitat_Par
    paste0(pref$y.groups.taxa, pref$y.habitats),                    #  6: Group1_Taxon2_Habitat_Par
    paste0(pref$y.taxa.groups, pref$y.habitats),                    #  7: Taxon1_Group2_Habitat_Par
    paste0(pref$y.groups.groups, pref$y.habitats),                  #  8: Group1_Group2_Habitat_Par
    paste0(pref$y.taxa.taxa, pref$y.reaches),                       #  9: Taxon1_Taxon2_Reach_Par
    paste0(pref$y.groups.taxa, pref$y.reaches),                     # 10: Group1_Taxon2_Reach_Par
    paste0(pref$y.taxa.groups, pref$y.reaches),                     # 11: Taxon1_Group2_Reach_Par
    paste0(pref$y.groups.groups, pref$y.reaches),                   # 12: Group1_Group2_Reach_Par
    paste0(pref$y.taxa.taxa),                                       # 13: Taxon1_Taxon2_Par
    pref$y.groups.taxa,                                             # 14: Group1_Taxon2_Par
    paste0(pref$y.taxa.groups),                                     # 15: Taxon1_Group2_Par
    pref$y.groups.groups,                                           # 16: Group1_Group2_Par
    paste0(pref$y.taxa, pref$y.reaches, pref$y.habitats),           # 17: Taxon1_Reach_Habitat_Par
    paste0(pref$y.groups, pref$y.reaches, pref$y.habitats),         # 18: Group1_Reach_Habitat_Par
    paste0(pref$y.taxa, pref$y.habitats),                           # 19: Taxon1_Habitat_Par
    paste0(pref$y.groups, pref$y.habitats),                         # 20: Group1_Habitat_Par
    paste0(pref$y.taxa, pref$y.reaches),                            # 21: Taxon1_Reach_Par
    paste0(pref$y.groups, pref$y.reaches),                          # 22: Group1_Reach_Par
    pref$y.taxa,                                                    # 23: Taxon1_Par
    pref$y.groups,                                                  # 24: Group1_Par
    paste0(pref$y.reaches, pref$y.habitats),                        # 25: Reach_Habitat_Par
    pref$y.habitats,                                                # 26: Habitat_Par
    pref$y.reaches,                                                 # 27: Reach_Par
    "")                                                             # 28: Par

  vals.inds <- get.inpind.parval.vec(
    pref.names.cols, par.names, pairs$y.names, par, inp, defaults, query.intvec)

  return(vals.inds)
}


# check existence of required parameters
# --------------------------------------

check.inpind.parval.required <- function(
  inpinds.parvals, required = NA, common.no.val.suffix = "")
{

  inds <- inpinds.parvals$inpinds
  vals <- inpinds.parvals$parvals

  if.vals.1d <- is.null(dim(vals))
  if (if.vals.1d) {
    stopifnot(ncol(inds) == 2)
    par.names <- names(vals)
    # make vals a row-vector (matrix)
    dim(vals) <- c(1, length(vals))
    # append 1 "i" column to inds
    inds <- cbind(inds, rep(1, nrow(inds)))
    inds <- inds[, c(1, 3, 2), drop=FALSE]
  } else {
    par.names <- colnames(vals)
  }
  var.no.val.suffix <- length(common.no.val.suffix) == 0 || !if.vals.1d

  for ( j.req in 1:length(required) )
  {
    if ( !is.na(required[j.req]) )
    {
      j <- ( match(required[j.req], par.names) )
      if ( is.na(j) )
      {
        warning("Parameter \"", required[j.req],
                "\" specified to be required ",
                "but not found in provided parameter names",
                sep = "")
      } else {
        for ( i in 1:nrow(vals) )
        {
          if (
            is.na(vals[i,j]) &&
            ( nrow(inds) == 0 || sum(inds[,2]==i & inds[,3]==j) == 0 )
          )
          {
            warning("Parameter \"", required[j.req],
                    "\" specified to be required ",
                    "but missing value both in input and in parameter vector",
                    if (var.no.val.suffix) paste0(" for state variable #", i)
                    else common.no.val.suffix,
                    sep = "")
          }
        }
      }
    }
  }

  invisible(NULL)
}


