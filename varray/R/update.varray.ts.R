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
    } else {
        adn <- rep(list(character(0)), length(dim(data)))
        ad <- rep(0, length(dim(data)))
    }

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

    if (is.null(along))
        stop("must specify 'along'")
    if (is.null(dimorder))
        dimorder <- seq(length=length(dd))
    if (is.null(naidxok))
        naidxok <- NA

    # find existing sub-components
    ex.ai <- match(ddn[[along]], adn[[along]])
    # existing sub-component names
    ex.scn <- as.character(sapply(va$info, '[[', 'name'))
    # If we need any new slices, create and/or expand existing subcomponents
    if (length(new.slices <- which(is.na(ex.ai))) || is.null(va)) {
        # new.slices is integer the new slices we need in 'data'
        new.scn <- format(dateParse(ddn[[along]][new.slices]), format=comp.name)
        all.scn.u <- unique(c(ex.scn, new.scn))
        new.scn.u <- setdiff(unique(new.scn), ex.scn)
        new.slices.sci <- match(new.scn, all.scn.u)
        va2 <- structure(list(dim=NULL, dimnames=NULL, along=along,
                              info=rep(list(name=NULL, dim=NULL, dimnames=NULL, env.name=NULL,
                                       sample=NULL, naidxok=NULL, map=NULL), length(all.scn.u)),
                              along.idx=NULL, dimorder=NULL, naidxok=NULL,
                              comp.name=NULL, env.name=NULL), class='varray')
        # if any new sub-components are needed, create them
        for (this.new.scn in new.scn.u) {
            this.data <- asub(data, new.scn==this.new.scn, dims=along, drop=FALSE)
            this.data.dn <- dimnames(this.data)
            # find dimnames that don't have all NA values
            for (i in seq(length(dim(data))[-along]))
                this.data.dn[[i]] <- this.data.dn[[i]][apply(this.data, i, function(x) !all(is.na(x)))]
            # do we need to subset this.data down to non-NA data?
            if (!isTRUE(all.equal(dimnames(this.data), this.data.dn)))
                this.data <- asub(this.data, this.data.dn, dims=seq(length=length(dim(this.data))), drop=FALSE)
            # do we need to reverse-permute the data?
            rdimorder <- order(dimorder)
            if (all(dimorder==seq(len=length(dim)))) {
                datar <- data
            } else {
                datar <- aperm(data, rdimorder)
            }
            j <- match(this.new.scn, all.scn.u)
            va2$info[[j]]$name <- this.new.scn
            va2$info[[j]]$dim <- dim(datar)
            va2$info[[j]]$dimnames <- dimnames(datar)
            va2$info[[j]]$env.name <- env.name
            va2$info[[j]]$naidxok <- NA
            # will fix 'map' at the end
            assign(this.new.scn, envir=envir, value=datar)
        }
        # if we need to add slices to any existing sub-components, do that
    }
    browser()

}
