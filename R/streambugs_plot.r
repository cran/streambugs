# function to generate standard plot of streambugs results
# --------------------------------------------------------

# Note: S3-inheritance arguments consistency with generic plot function as the
#       simulation run result `x` has class "streambugs" set. In turn, calling
#       `plot(res, par, inp, ...)` will invoke this method

#' Plot the results of streambugs ODE run
#'
#' Plot time series of all streambugs ODE state variables, for each reach,
#' habitat and group, resulting from the
#' \code{\link{run.streambugs}} function call.
#'
#' @param x matrix with results derived by
#'    \code{\link{run.streambugs}}
#' @param y same as \code{par} in \code{\link{run.streambugs}}
#' @param inp same as \code{inp} in \code{\link{run.streambugs}}
#' @param ... additional argument for the \code{\link[graphics]{plot}} function
#'    call
#'
#' @examples
#' m <- streambugs.example.model.toy()
#' r <- run.streambugs(y.names=m$y.names, times=m$times,  par=m$par, inp=m$inp, C=TRUE)
#' plot(x=r$res, y=m$par, inp=m$inp)
#'
#' @export
plot.streambugs <- function(x,y,inp=NA,...)
{
  res = x; par = y
  sys.def <- streambugs.get.sys.def(y.names=colnames(res)[-1],par=par,inp=inp)
  y.names <- sys.def$y.names

  par.envcond.w <- get.inpind.parval.envcond.reach(
    par.names = c("w"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    required  = c("w"))

  par.envcond.fA <- get.inpind.parval.envcond.habitat(
    par.names = c("fA"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    defaults  = c(fA=1))

  # evaluate time-dependent inputs and merge them with parameters
  # (overwriting parameters if one exists with the same name):
  # -------------------------------------------------------------

  for ( j in 1:nrow(res) )
  {
    # update (time-dependent) parameters:

    inpvals <- interpolate.inputs(inp,res[j,1])

    w  <- streambugs.update.envcond.reach(par.envcond.w,inpvals)[,"w"]
    fA <- streambugs.update.envcond.habitat(par.envcond.fA,inpvals,y.names$ind.fA)[,"fA"]

    # convert results in mass per unit length into mass per unit
    # surface area:

    res[j,-1] <- res[j,-1]/(w*fA)
  }

  # determin maxima of converted results:

  y.max <- rep(NA,length(y.names$groups))
  for ( i in 1:length(y.names$groups) )
  {
    y.max[i] <- 1.1*max(res[,1+which(y.names$y.groups==y.names$groups[i])])
  }
  names(y.max) <- y.names$groups

  # get time:

  t    <- res[,1]

  # plot output by groups:

  par.def <- graphics::par(no.readonly=TRUE)
  graphics::par(mfrow=c(length(y.names$habitats),length(y.names$groups)))
  for ( reach in y.names$reaches )
  {
    for ( habitat in y.names$habitats )
    {
      for ( group in y.names$groups )
      {
        plot(numeric(0),numeric(0),type="n",
             xlim=c(min(t),max(t)),ylim=c(0,y.max[group]),
             xlab="time [a]",ylab="biomass [gDM/m2]",
             main=paste(reach,habitat,group), ...)
        ind <- which(y.names$y.reaches  == reach &
                       y.names$y.habitats == habitat &
                       y.names$y.groups   == group )
        if ( length(ind) > 0 )
        {
          for ( k in 1:length(ind) )
          {
            lines(t,res[,1+ind[k]],lty=k)
          }
          legend(x="topleft",legend=y.names$y.taxa[ind],lty=1:length(ind))
        }
      }
    }
  }
  graphics::par(par.def)

  #write.table(res,paste("output/res_gDMperm2_",name.run,".dat",sep=""),
  #            sep="\t",row.names=F,col.names=T)

}


#' Plot the streambugs "foodweb" graph.
#'
#' Plot the "foodweb" graph depicting interactions between ODE variables in a streambugs
#' model.
#'
#' @param y.names same as \code{y.names} in \code{\link{run.streambugs}}
#' @param par same as \code{par} in \code{\link{run.streambugs}}
#' @param file optional name of a PDF file to plot to
#' @param cex same as \code{cex} in \code{\link[graphics]{par}}, consumed by here by
#'     multiple text generating functions
#' @param font same as \code{font} in \code{\link[graphics]{par}}, consumed here by
#'     \code{\link[graphics]{text}} as a font type for taxa names
#' @param title optional title for the plot
#' @param lwd same as \code{lwd} in \code{\link[graphics]{par}}, consumed here by
#'     \code{\link[graphics]{lines}} as a line width for the "food web" edges
#' @param bg same as \code{bg} in \code{\link[graphics]{par}}
#' @param lcol same as \code{col} in \code{\link[graphics]{par}}, consumed here by
#'     \code{\link[graphics]{lines}} as a line color for the "food web" edges
#' @param ncrit number of inverts in one line at which they are shifted up and down,
#'     alternating
#' @param lcrit number of letters/characters of state names which are plottet at level 2
#' @param survivals vector with the entries \code{"survived"} or \code{"extinct"} for
#'     all state variables with names matching names of taxa state variables and
#'     \code{"SusPOM"}
#' @param observed vector with the entries \code{"never"}/\code{"notobserved"},
#'     \code{"observed"},\code{"sometimes"}, \code{"always"}, or \code{NA}/\code{"NA"}
#'     for all state variables, with names matching names of taxa state variables and
#'     \code{"SusPOM"}
#' @param texts if to plot as "food web"nodes as texts with taxa names; otherwise, plot
#'     points
#' @param pointcol if to color text or point nodes using Eawag coloring scheme
#'     (eawagfarben); otherwise plot all in the same color
#' @param ... additional arguments for the \code{\link[grDevices]{pdf}} graphics device,
#'     relevant only if \code{file} argument was given.
#'
#' @examples
#' model <- streambugs.example.model.toy()
#' foodweb.plot(model$y.names, model$par, cex=1.1, title="complete foodweb", ncrit=8,
#'     lcrit=7, lwd=2, bg="white", lcol="blue", font=2)
#'
#' @export
foodweb.plot <- function(y.names,par,file=NA,cex=1,font=1,title="",
                         lwd=1,bg=colors()[1],lcol=colors()[555],ncrit=8,lcrit=20,
                         survivals=NA,observed=NA,texts=TRUE,pointcol=FALSE,...) {

  #analyse structure
  # ------------------

  # decode states:

  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

  state_names <- c("SusPOM",y.names$taxa)

  # get levels and food of each state variable:

  stoich.Cons  <- get.par.stoich.web("Cons", par) #xxx

  level         <-  rep(1,length(state_names))
  names(level)  <-  state_names
  foods         <-  list()

  # identify foods
  # loop over state variables:

  for ( i in 1:length(state_names) )
  {
    state_name <- state_names[i]

    foods[[state_name]] <- names(stoich.Cons[[state_name]])

  }  #loop over state variables

  # set levels of consumers

  count = 1
  repeat
  {
    #loop over state variables
    for (i in 1:length(foods) )
    {

      state_name <- names(foods[i])

      #update own level

      level[state_name] <- max( (max( level[paste(foods[[state_name]])] ) +1) ,level[state_name] )

      #update level of my predators

      predators<-NA
      for (  j in 1:length(foods))
      {
        pind <- which(foods[[j]]==state_name)
        if (sum(pind) >0) {predators <- c( predators, names(foods[j]) )}
      }
      predators    <- unique(predators[!is.na(predators)] )


      if (length(predators)>0)
      {
        for (p in 1:length(predators))
        {
          predator <- predators[p]
          level[predator] <- max( (level[state_name]+1),level[predator] )
        }
      }
    }

    count <- count+1

    if(count>10){break}

  } #repeat loop



  # ( level )
  # ( foods  )

  structure           <- matrix(NA,nrow=length(state_names),ncol=4)
  colnames(structure) <- c("level","width","xpos","ypos")
  rownames(structure) <- state_names
  structure           <- as.data.frame(structure)
  structure[,"level"] <- level
  structure[,"width"] <- 1

  ind.topdown         <- order(structure$level,decreasing=F)
  names(ind.topdown)  <- state_names[ind.topdown]
  ind.bottomup        <- order(structure$level,decreasing=T)
  names(ind.bottomup) <- state_names[ind.bottomup]

  nr_levels <- NA           #nr of states at each level
  for (l in 1:max(level) )
  {
    nr_levels[l] <-  length(structure[structure[,"level"]==l,"level"])
  }

  oldlevel <- 0
  offset   <- 0
  for ( stat in ind.bottomup )
  {
    #stat <- ind.bottomup[2]

    newlevel <- structure$level[stat]
    n        <- nr_levels[newlevel]

    if ( newlevel != oldlevel )
    {
      offset   <- 0
      structure$xpos[stat] <- offset + 1/(n+1)
      offset   <- structure$xpos[stat]
      oldlevel <- newlevel
    } else
    {
      structure$xpos[stat] <- offset + 1/(n+1)
      offset <- structure$xpos[stat]
    }
  }

  structure$ypos <- structure$level /(max(level)+1)

  ypos <- unique(structure$ypos)
  for (i in 1:length(ypos))
  {
    yposi <- ypos[i]
    indy <- which(structure$ypos==yposi)
    if( length(indy)> ncrit)
    {
      indx <- grep(paste(sort(structure$xpos[indy]),collapse="|"),structure$xpos[indy])

      indh <- which(indx%%2==0)
      indl <- which(indx%%2==1)
      structure$ypos[indy][indh] <- yposi*1.1
      structure$ypos[indy][indl] <- yposi*0.9
    }
  }


  # open pdf file:

  if ( !is.na(file) ) pdf(file,...)

  # initialize plot:

  # par.def <- graphics::par(no.readonly=TRUE)
  graphics::par(mar=c(0,0,0,0),bg=bg)
  plot(numeric(0),numeric(0),xlim=c(1,0),ylim=c(0,1),xlab="",ylab="",axes=F)# ,...


  #check survivals
  if ( length(survivals) == 1 && is.na(survivals) )
  {
    no.survivals <- TRUE
    survivals <- rep(NA,length(state_names))
    names(survivals) <- state_names
  } else {
    no.survivals <- FALSE
    if ( !all(names(survivals) %in% state_names) )
      stop("Names of the survivals vector must match names of the state variables")
    if ( !all(survivals %in% c("survived", "extinct")) )
      stop("Values of the survivals vector must be either \"survived\", or \"extinct\"")
  }

  #check observed
  if ( length(observed) == 1 && is.na(observed))
  {
    no.observed <- TRUE
    observed <- rep(NA,length(state_names))
    names(observed) <- state_names
  } else {
    no.observed <- FALSE
    if ( !all(names(observed) %in% state_names) )
      stop("Names of the observed vector must match names of the state variables")
    if ( !all(observed %in% c(NA, "always", "never", "observed","notobserved", "sometimes")) )
      stop("Values of the observed vector must be either NA/\"NA\", \"always\", \"never\"/\"observed\",\"notobserved\", or \"sometimes\"")
  }


  # draw food connections

  for ( s in 1:length(foods) )
  {
    food <- foods[[s]]
    cons <- names(foods[s])
    if (!is.na(food)[1])
    {
      for (f in 1:length(food))
      {
        from <- cons
        to   <- food[f]
        x <- numeric(2)
        x[1]<- structure[paste(from),"xpos"]
        x[2]<- structure[paste(to),"xpos"]
        y   <- numeric(2)
        y[1]<- structure[paste(from),"ypos"]
        y[2]<- structure[paste(to),"ypos"]

        indcons <-  which(from==names(survivals))
        indconsm <- which(from==names(observed))

        indfood <- which(to==names(survivals))
        indfoodm <- which(to==names(observed))

        if( sum(indfoodm)==0 ) warning (to, " not in observed")
        if( sum(indfood)==0 )  warning (to, " not in survivals")

        cond.1 <- (survivals[indcons]=="extinct"|
                     survivals[indfood]=="extinct"|
                     observed[indconsm]=="never"  |
                     observed[indfoodm]=="never")

        cond.2 <- ( observed[indconsm]=="sometimes" |
                      observed[indfoodm]=="sometimes" |
                      observed[indconsm]=="NA"        |
                      observed[indfoodm]=="NA" )

        cond.3 <- !no.observed # length(observed)>1 & !is.na(observed)


        clf <- lcol

        if( !is.na(cond.1) & cond.1 ) { clf=NA } else {
          if( !is.na(cond.2) & cond.2 & cond.3) {clf=gray(0.5) } }

        lines(x,y,col=clf,lwd=lwd)

      }
    }
  }


  # plot text or points

  for ( l in 1 : max(level) )
  {
    label = rownames(structure)[structure$level == l]

    cl <- rep(1,length(label))

    for (m in 1:length(label))
    {

      if (!pointcol)
      {
        #cl[m] <- gray(0.5)

        if(!no.survivals)
        {
          inds <- which(label[m]==names(survivals))
          if( !is.na(survivals[inds]) )
          {
            if(survivals[inds]=="extinct")  { cl[m] <- gray(0.6) }
          }
        }

        if(!no.observed)
        {
          indm <- which(label[m]==names(observed))
          if( !is.na( observed[indm]) )
          {
            if(observed[indm]=="notobserved" | observed[indm]=="never")
            {cl[m] <- gray(0.7)} else {
              if(observed[indm]== "sometimes" )
              {cl[m] <- gray(0.45)} else {
                if (observed[indm]== "always" | observed[indm]=="observed")
                {cl[m] <- 1} else {
                  if (is.na(observed[indm]) | observed[indm]== "NA")
                  { cl[m] <- gray(0.5) }

                }
              }
            }
          }
        }

      } else
      {
        if (pointcol)
        {

          if(!no.survivals)
          {
            inds <- which(label[m]==names(survivals))
            if( !is.na(survivals[inds]) )
            {
              ifelse(survivals[inds]=="extinct",
                     cl[m] <- rgb(240,123,0,maxColorValue = 255),
                     cl[m] <- rgb(0,173,221,maxColorValue = 255))
            }
          }

          if(!no.observed)
          {
            indm <- which(label[m]==names(observed))
            if( !is.na( observed[indm] ) )
            {
              if(observed[indm]=="notobserved" | observed[indm]=="never")
              {cl[m] <- rgb(240,123,0,maxColorValue = 255)} else {
                if(observed[indm]== "sometimes" )
                {cl[m] <- rgb(188,208,45,maxColorValue = 255)} else {
                  if (observed[indm]== "always" | observed[indm]=="observed")
                  {cl[m] <- rgb(0,173,221,maxColorValue = 255)} else {
                    if (observed[indm]=="NA")
                    { cl[m] <-gray(0.5)}
                  }
                }
              }
            }
          }

        }
      }
    } #end for m

    if(texts)
    {
      if(l==2) {
        label <- substr(rownames(structure)[structure$level == l],1,lcrit)
        label <- ifelse(nchar(rownames(structure)[structure$level == l])>lcrit,
                        paste0(label,"."),label)
      }

      if( length(label)!=length(unique(label))  ) {warning("lcrit: ",lcrit," to small,labels not unique")}

      text(x  = structure$xpos[structure$level == l],
           y  = structure$ypos[structure$level == l],
           labels = label,
           col    = cl,
           font   = font,
           cex    = cex
      )

    }  else
    {
      points(x   = structure$xpos[structure$level == l],
             y   = structure$ypos[structure$level == l],
             pch = 19,
             col = cl,
             cex = cex*2
      )
    }

  } #end for l

  mtext(title,3,line=-2,cex=cex,font=2)

  # reset plot parameters:
  #graphics::par(par.def)

  # close pdf file:

  if ( !is.na(file) ) dev.off()

  # return(list(structure=structure[,c(1,4)])) #foods=foods,products=products

} # end function

