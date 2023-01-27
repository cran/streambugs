c_streambugs_init_debug <- function(streambugs.debug) {
    .C("streambugs_init_debug", PACKAGE = 'streambugs',
        as.integer(streambugs.debug)
    )
}

c_streambugs_create_input_structure <- function(ninp) {
    .C("streambugs_create_input_structure", PACKAGE = 'streambugs',
        as.integer(ninp)
    )
}

c_streambugs_create_input <- function(i, n, x, y) {
    .C("streambugs_create_input", PACKAGE = 'streambugs',
        as.integer(i),
        as.integer(n),
        as.double(x),
        as.double(y)
    )
}

c_streambugs_create_parglobal <- function(n_inp, inpinds, n_par, parvals, parnames) {
    .C("streambugs_create_parglobal", PACKAGE = 'streambugs',
        as.integer  (n_inp),
        as.integer  (inpinds),
        as.integer  (n_par),
        as.double   (parvals),
        as.character(parnames)
    )
}

c_streambugs_create_parglobalenvtraits_structure <- function(ntrait) {
    .C("streambugs_create_parglobalenvtraits_structure", PACKAGE = 'streambugs',
        as.integer(ntrait)
    )
}

c_streambugs_create_parglobalenvtrait <- function(i, name, nind, inds, npar, vals, nams) {
    .C("streambugs_create_parglobalenvtrait", PACKAGE = 'streambugs',
         as.integer  (i),
         as.character(name),
         as.integer  (nind),
         as.integer  (inds),
         as.integer  (npar),
         as.double   (vals),
         as.character(nams)
    )
}

c_streambugs_create_parenvcondreach <- function(nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_parenvcondreach", PACKAGE = 'streambugs',
        as.integer  (nind),
        as.integer  (inds),
        as.integer  (npar),
        as.integer  (nrow),
        as.double   (vals),
        as.character(nams)
    )
}

c_streambugs_create_parenvcondhabitat <- function(nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_parenvcondhabitat", PACKAGE = 'streambugs',
        as.integer  (nind),
        as.integer  (inds),
        as.integer  (npar),
        as.integer  (nrow),
        as.double   (vals),
        as.character(nams)
    )
}

c_streambugs_create_parenvcondhabitatgroups_structure <- function(ngroup) {
    .C("streambugs_create_parenvcondhabitatgroups_structure", PACKAGE = 'streambugs',
        as.integer(ngroup)
    )
}

c_streambugs_create_parenvcondhabitatgroup <- function(i, name, nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_parenvcondhabitatgroup", PACKAGE = 'streambugs',
         as.integer  (i),
         as.character(name),
         as.integer  (nind),
         as.integer  (inds),
         as.integer  (npar),
         as.integer  (nrow),
         as.double   (vals),
         as.character(nams)
    )
}

c_streambugs_create_parinitcond <- function(nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_parinitcond", PACKAGE = 'streambugs',
        as.integer  (nind),
        as.integer  (inds),
        as.integer  (npar),
        as.integer  (nrow),
        as.double   (vals),
        as.character(nams)
    )
}

c_streambugs_create_parinput <- function(nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_parinput", PACKAGE = 'streambugs',
        as.integer  (nind),
        as.integer  (inds),
        as.integer  (npar),
        as.integer  (nrow),
        as.double   (vals),
        as.character(nams)
    )
}

c_streambugs_create_partaxapropdirect <- function(nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_partaxapropdirect", PACKAGE = 'streambugs',
        as.integer  (nind),
        as.integer  (inds),
        as.integer  (npar),
        as.integer  (nrow),
        as.double   (vals),
        as.character(nams)
    )
}

c_streambugs_create_partaxaproptraits_structure <- function(ntrait) {
    .C("streambugs_create_partaxaproptraits_structure", PACKAGE = 'streambugs',
        as.integer(ntrait)
    )
}

c_streambugs_create_partaxaproptrait <- function(i, name, nind, inds, npar, nrow, vals, nams) {
    .C("streambugs_create_partaxaproptrait", PACKAGE = 'streambugs',
         as.integer  (i),
         as.character(name),
         as.integer  (nind),
         as.integer  (inds),
         as.integer  (npar),
         as.integer  (nrow),
         as.double   (vals),
         as.character(nams)
    )
}

c_streambugs_create_processes_structure <- function(ny) {
    .C("streambugs_create_processes_structure", PACKAGE = 'streambugs',
        as.integer(ny)
    )
}

c_streambugs_create_proctaxon <- function(i, j, procname, ninp, inpinds, npar,
    parnames, parvals, nstoich, stoichnames, stoichinds, stoichvals) {
    .C("streambugs_create_proctaxon", PACKAGE = 'streambugs',
        as.integer(i),
        as.integer(j),
        as.character(procname),
        as.integer(ninp),
        as.integer(inpinds),
        as.integer(npar),
        as.character(parnames),
        as.double(parvals),
        as.integer(nstoich),
        as.character(stoichnames),
        as.integer(stoichinds),
        as.double(stoichvals)
    )
}

c_streambugs_create_procweb <- function(i, j, procname, ninp, inpinds, npar, parnames, parvals) {
    .C("streambugs_create_procweb", PACKAGE = 'streambugs',
        as.integer(i),
        as.integer(j),
        as.character(procname),
        as.integer(ninp),
        as.integer(inpinds),
        as.integer(npar),
        as.character(parnames),
        as.double(parvals)
    )
}

c_streambugs_create_procwebtaxon <- function(i, j, k, procname, ninp, inpinds,
    npar, parnames, parvals, nstoich, stoichnames, stoichinds, stoichvals) {
    .C("streambugs_create_procwebtaxon", PACKAGE = 'streambugs',
        as.integer(i),
        as.integer(j),
        as.integer(k),
        as.character(procname),
        as.integer(ninp),
        as.integer(inpinds),
        as.integer(npar),
        as.character(parnames),
        as.double(parvals),
        as.integer(nstoich),
        as.character(stoichnames),
        as.integer(stoichinds),
        as.double(stoichvals)
    )
}

c_streambugs_create_fA_structure <- function(nreach) {
    .C("streambugs_create_fA_structure", PACKAGE = 'streambugs',
        as.integer(nreach)
    )
}

c_streambugs_create_fA <- function(i, nreachind, reachind, nfsthabind, fsthabind) {
    .C("streambugs_create_fA", PACKAGE = 'streambugs',
       as.integer(i),
       as.integer(nreachind),
       as.integer(reachind),
       as.integer(nfsthabind),
       as.integer(fsthabind)
    )
}

c_streambugs_delete_inputs <- function()
    .C("streambugs_delete_inputs", PACKAGE = 'streambugs')

c_streambugs_delete_parglobal <- function()
    .C("streambugs_delete_parglobal", PACKAGE = 'streambugs')

c_streambugs_delete_parglobalenvtraits <- function()
    .C("streambugs_delete_parglobalenvtraits", PACKAGE = 'streambugs')

c_streambugs_delete_parenvcondreach <- function()
    .C("streambugs_delete_parenvcondreach", PACKAGE = 'streambugs')

c_streambugs_delete_parenvcondhabitat <- function()
    .C("streambugs_delete_parenvcondhabitat", PACKAGE = 'streambugs')

c_streambugs_delete_parenvcondhabitatgroup <- function()
    .C("streambugs_delete_parenvcondhabitatgroup", PACKAGE = 'streambugs')

c_streambugs_delete_parinitcond <- function()
    .C("streambugs_delete_parinitcond", PACKAGE = 'streambugs')

c_streambugs_delete_parinput <- function()
    .C("streambugs_delete_parinput", PACKAGE = 'streambugs')

c_streambugs_delete_partaxapropdirect <- function()
    .C("streambugs_delete_partaxapropdirect", PACKAGE = 'streambugs')

c_streambugs_delete_partaxaproptraits <- function()
    .C("streambugs_delete_partaxaproptraits", PACKAGE = 'streambugs')

c_streambugs_delete_processes <- function()
    .C("streambugs_delete_processes", PACKAGE = 'streambugs')

c_streambugs_delete_fA <- function()
    .C("streambugs_delete_fA", PACKAGE = 'streambugs')
