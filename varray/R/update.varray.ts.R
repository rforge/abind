update.varray.ts <- function(va, data, comp.name=va$comp.name, dateblock='%Y', dates.by='bizdays', holidays='NYSEC',
                             vmode='single', along=va$along, dimorder=va$dimorder,
                             env.name=va$env.name, envir=NULL, naidxok=va$naidxok, character.only=FALSE) {
    # Update a varray that stores time series matrix data.
    # Labels on the binding 'along' dimension are dates, stored in order
    # Might need to create new component arrays in the varray.
    # Might need to create new columns in the varray -- though first
    # try to use up ones that have NA labels.
    # Return a varray, which is an in-memory object that refers to its components,
    # and has complete dimnames.
    # There might be some modification of underlying objects.

    # Work out what dates in data are new, and which are old
    # For old dates in data, find the objects to assign into,
    # increase their dims if necessary, and assign values into
    # them.

    # get 'va' as NULL or the varray
    if (is.character(substitute(va)) || character.only) {
        va.name <- va
        if (!is.character(va.name))
            stop('va must be supplied as character data when character.only=TRUE')
        if (exists(va.name)) {
            va <- get(va.name)
        } else {
            va <- NULL
        }
    } else {
        va.name <- '**unknown**'
        force(va)
    }
    # and 'adn' and 'ad'; the dimnames & dims of 'va'
    if (!is.null(va)) {
        if (!inherits(va, 'varray'))
            stop('va is not a varray object')
        adn <- dimnames(va)
        ad <- dim(va)
        for (i in seq(along=va$info)[-1]) {
            if (length(va$info[[i]]$dim) != length(va$info[[1]]$dim))
                stop('chunk ', i, ' has different dimensionality: ',
                     length(va$info[[i]]$dim), ' vs ', length(va$info[[1]]$dim), ' in first')
        }
    } else {
        adn <- rep(list(character(0)), length(dim(data)))
        ad <- rep(0, length(dim(data)))
    }
    if (along < 1 || along > length(non.null(va$info[[1]]$dim, dim(data))))
        stop('along must be in 1..', length(non.null(va$info[[1]]$dim), dim(data)))

    if (is.null(comp.name))
        comp.name <- paste(va.name, dateblock, sep='.')

    # get 'envir' and 'env.name'
    if (!is.null(envir)) {
        env.name <- fixGlobalEnvName(environmentName(envir))
    } else if (is.null(env.name)) {
        envir <- .GlobalEnv
        env.name <- fixGlobalEnvName(environmentName(envir))
    } else {
        envir <- as.environment(env.name)
    }

    # 'ddn' and 'dd' are dimnames and dims of data (the new data)
    ddn <- dimnames(data)
    dd <- dim(data)
    if (!is.null(va) && length(va$info[[1]]$dim) != length(dd))
        stop('component 1 has different dimensionality than data: ',
                     length(va$info[[1]]$dim), ' vs ', length(dd))

    if (is.null(along))
        stop("must specify 'along'")
    if (is.null(dimorder))
        dimorder <- seq(length=length(dd))
    if (is.null(naidxok))
        naidxok <- NA
    rdimorder <- order(dimorder)

    # find existing sub-components
    ex.ai <- match(ddn[[along]], adn[[along]])
    # existing sub-component names
    ex.scn <- as.character(sapply(va$info, '[[', 'name'))
    # If we need any new slices, create and/or expand existing subcomponents
    new.slices <- which(is.na(ex.ai))
    if (length(new.slices) || is.null(va)) {
        # new.slices is integer: the slices in 'data' that need new subcomponents
        new.scn <- format(dateParse(ddn[[along]][new.slices]), format=comp.name)
        all.scn.u <- unique(c(ex.scn, new.scn))
        new.scn.u <- setdiff(unique(new.scn), ex.scn)
        new.slices.sci <- match(new.scn, all.scn.u)
        comp.dn.changed <- rep(FALSE, length(all.scn.u))
        if (is.null(va)) {
            va <- structure(list(dim=NULL, dimnames=NULL, along=along,
                                 info=rep(list(list(name=NULL, dim=NULL, dimnames=NULL, env.name=NULL,
                                          sample=NULL, naidxok=NULL, map=NULL)), length(all.scn.u)),
                                 along.idx=NULL, dimorder=dimorder, naidxok=naidxok,
                                 comp.name=comp.name, env.name=env.name), class='varray')
        } else {
            va$info <- c(va$info, rep(list(list(name=NULL, dim=NULL, dimnames=NULL, env.name=NULL,
                                                sample=NULL, naidxok=NULL, map=NULL)), length(new.scn.u)))
        }
        # if any new sub-components are needed, create them
        for (this.new.scn in new.scn.u) {
            this.data <- asub(data, new.slices[new.scn==this.new.scn], dims=along, drop=FALSE)
            this.data.dn <- dimnames(this.data)
            # find dimnames that don't have all NA values
            for (i in seq(length(dim(data))[-along]))
                this.data.dn[[i]] <- this.data.dn[[i]][apply(this.data, i, function(x) !all(is.na(x)))]
            # do we need to subset this.data down to non-NA data?
            if (!isTRUE(all.equal(dimnames(this.data), this.data.dn)))
                this.data <- asub(this.data, this.data.dn, dims=seq(length=length(dim(this.data))), drop=FALSE)
            # do we need to reverse-permute the data?
            if (all(dimorder==seq(len=length(dim)))) {
                this.datar <- this.data
            } else {
                this.datar <- aperm(this.data, rdimorder)
            }
            j <- match(this.new.scn, all.scn.u)
            va$info[[j]]$name <- this.new.scn
            va$info[[j]]$dim <- dim(this.datar)
            va$info[[j]]$dimnames <- dimnames(this.datar)
            va$info[[j]]$env.name <- env.name
            va$info[[j]]$naidxok <- naidxok
            # will fix 'map' at the end
            assign(this.new.scn, envir=envir, value=this.datar)
        }
        # if we need to add slices to any existing sub-components, do that
    } else {
        comp.dn.changed <- rep(FALSE, length(va$info))
        new.scn <- character(0)
    }
    # data.ii is the slices of 'data' that are to be found in existing components of 'va'
    data.ii <- which(!is.na(ex.ai))
    if (length(data.ii)) {
        # work out which existing components of 'va' we need to work with
        exist.comp.i <- va$along.idx[ex.ai[data.ii]]
        for (this.i in unique(exist.comp.i)) {
            # working with component 'this.i' of 'va', and data slices 'data.i'
            data.i <- data.ii[this.i == exist.comp.i]
            # add in the indices that are new in the 'along' dimension but belong to this component
            data.i <- union(data.i, which(is.na(ex.ai))[new.scn == va$info[[this.i]]$name])
            this.data <- asub(data, data.i, dims=along, drop=FALSE)
            this.data.dn <- dimnames(this.data)
            # find dimnames that don't have all NA values
            for (i in seq(length(dim(data))[-along]))
                this.data.dn[[i]] <- this.data.dn[[i]][apply(this.data, i, function(x) !all(is.na(x)))]
            # do we need to subset this.data down to non-NA data?
            if (!isTRUE(all.equal(dimnames(this.data), this.data.dn)))
                this.data <- asub(this.data, this.data.dn, dims=seq(length=length(dim(this.data))), drop=FALSE)
            # do we need to reverse-permute the data?
            if (all(dimorder==seq(len=length(dim)))) {
                this.datar <- this.data
                this.datar.dn <- this.data.dn
            } else {
                this.datar <- aperm(this.data, rdimorder)
                this.datar.dn <- this.data.dn[rdimorder]
            }
            need.expand <- mapply(va$info[[this.i]]$dimnames, this.datar.dn, FUN=function(old, new) !all(is.element(new, old)))
            env <- as.environment(va$info[[this.i]]$env.name)
            comp.data <- get(va$info[[this.i]]$name, envir=env, inherits=FALSE)
            if (any(need.expand)) {
                comp.dn.changed[this.i] <- TRUE
                new.dn <- mapply(va$info[[this.i]]$dimnames, this.datar.dn, FUN=function(old, new) union(old, new))
                comp.data <- conform(comp.data, new.dn)
                va$info[[this.i]]$dimnames <- new.dn
                va$info[[this.i]]$dim <- sapply(new.dn, length)
            }
            # load values into the component
            afill(comp.data) <- this.datar
            # and save the component
            assign(va$info[[this.i]]$name, envir=env, value=comp.data, inherits=FALSE)
        }
    }
    dn <- lapply(seq(len=length(va$info[[1]]$dim)), function(i)
                 unique(unlist(lapply(va$info, function(x) x$dimnames[[i]]))))
    d <- sapply(dn, length)
    va$along.idx <- integer(d[along])
    naidxok <- all(sapply(va$info, '[[', 'naidxok'))
    if (is.null(dimorder))
        dimorder <- seq(length(d))
    else
        if (!identical(sort(as.numeric(dimorder)), as.numeric(seq(length(d)))))
            stop('dimorder must be 1:length(d) in some permutation')
    # convert d,dn to user dimorder
    if (!all(dimorder == seq(length(d)))) {
        d <- d[dimorder]
        dn <- dn[dimorder]
    }
    rdimorder <- order(dimorder)
    alongd <- rdimorder[along]
    # fix 'map' in all info components
    # eventualy, record which components were changed, and only update those
    for (i in seq(to=1, from=length(va$info))) {
        va$along.idx[match(va$info[[i]]$dimnames[[alongd]], dn[[along]])] <- i
        va$info[[i]]$map <- lapply(seq(along=dn), function(j) match(dn[[j]], va$info[[i]]$dimnames[[rdimorder[j]]]))[rdimorder]
    }
    va$dim <- d
    va$dimnames <- dn
    va
}
