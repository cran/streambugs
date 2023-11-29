# load streambugs library:
# ========================
library(streambugs)

# handler to suppress specific warnings from the streambugs library
.streambugs.suppress.warning.handler <- function(w) {
    if(any(grepl("(no (envcond|parameters) for .*|limitation by (\\w)+( \\w+)? not considered)", w))) {
        invokeRestart("muffleWarning")
    }
}

model <- streambugs.example.model.extended()


# run C implementation:
# =====================

# re-compile C routines of streambugs
# after having done modifications of C code:

#compile.streambugs()

# initialize global parameter and parameter name vectors for use in C
# and calculate results:

withCallingHandlers(
  res.C <- run.streambugs(y.names  = model$y.names,
                          times    = model$times,
                          par      = model$par,
                          inp      = model$inp,
                          C        = TRUE),
  warning = .streambugs.suppress.warning.handler
)

# res.C <- run.streambugs(y.names  = model$y.names,
#                         times    = model$times,
#                         par      = model$par,
#                         inp      = model$inp,
#                         C        = TRUE,
#                         file.def = "output/streambugs_extended_sysdef.dat",
#                         file.res = "output/streambugs_extended_res_C.tsv")

res.C$args = c(model, C = TRUE)
