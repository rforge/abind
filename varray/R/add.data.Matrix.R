# Not a generic because the first argument is a character string, so methods
# do not have a class of object to dispatch on.
# add.data <- function(object, data, need.dimnames=list(NULL, NULL), keep.ordered=TRUE, ...) UseMethod('add.data')

add.data.Matrix <- function(object, data, need.dimnames=list(NULL, NULL), keep.ordered=TRUE, ...) {
    # have ... args to satisfy the generic add.data()
    if (length(list(...)))
        warning('additional arguments ignored: ', paste(names(list(...)), collapse=', '))
    va.name <- object
    if (!is.character(va.name))
        stop('object must be supplied as character data naming matrix to be updated')
    if (!exists(va.name)) {
        x <- Matrix(data, sparse=TRUE)
    } else {
        x <- get(va.name)
    }
    dn <- list(unique(c(rownames(x), rownames(data), need.dimnames[[1]])),
               unique(c(colnames(x), colnames(data), need.dimnames[[2]])))
    keep.ordered <- rep(keep.ordered, length.out=length(dim(x)))
    if (any(keep.ordered))
        dn[keep.ordered] <- lapply(dn[keep.ordered], sort, na.last=NA)
    sample <- replace(x[1,1], 1, 0L)
    if (!isTRUE(all.equal(rownames(x), dn[[1]]))) {
        new <- setdiff(dn[[1]], rownames(x))
        if (length(new))
            x <- rBind(x, array(sample, dimnames=list(new, colnames(x)), dim=c(length(new), ncol(x))))
        if (!isTRUE(all.equal(rownames(x), dn[[1]])))
            x <- x[dn[[1]],,drop=FALSE]
    }
    if (!isTRUE(all.equal(colnames(x), dn[[2]]))) {
        new <- setdiff(dn[[2]], colnames(x))
        if (length(new))
            x <- cBind(x, array(sample, dimnames=list(rownames(x), new), dim=c(nrow(x), length(new))))
        if (!isTRUE(all.equal(colnames(x), dn[[2]])))
            x <- x[,dn[[2]],drop=FALSE]
    }
    ii <- cbind(rep(match(rownames(data), rownames(x)), ncol(data)),
                rep(match(colnames(data), colnames(x)), each=nrow(data)))
    x[ii] <- as.vector(data)
    assign(va.name, value=x, pos=1)
}

add.data.matrix <- function(object, data, need.dimnames=list(NULL, NULL), keep.ordered=TRUE, ...) {
    # have ... args to satisfy the generic add.data()
    if (length(list(...)))
        warning('additional arguments ignored: ', paste(names(list(...)), collapse=', '))
    va.name <- object
    if (!is.character(va.name))
        stop('object must be supplied as character data naming matrix to be updated')
    if (!exists(va.name)) {
        if (!is.matrix(data))
            stop('data must be a matrix')
        x <- data
    } else {
        x <- get(va.name)
        if (!is.matrix(x))
            stop(va.name, ' already exists but is not a matrix')
    }
    dn <- list(unique(c(rownames(x), rownames(data), need.dimnames[[1]])),
               unique(c(colnames(x), colnames(data), need.dimnames[[2]])))
    keep.ordered <- rep(keep.ordered, length.out=length(dim(x)))
    if (any(keep.ordered))
        dn[keep.ordered] <- lapply(dn[keep.ordered], sort, na.last=NA)
    sample <- replace(x[1,1], 1, 0L)
    if (!isTRUE(all.equal(rownames(x), dn[[1]]))) {
        new <- setdiff(dn[[1]], rownames(x))
        if (length(new))
            x <- rbind(x, array(sample, dimnames=list(new, colnames(x)), dim=c(length(new), ncol(x))))
        if (!isTRUE(all.equal(rownames(x), dn[[1]])))
            x <- x[dn[[1]],,drop=FALSE]
    }
    if (!isTRUE(all.equal(colnames(x), dn[[2]]))) {
        new <- setdiff(dn[[2]], colnames(x))
        if (length(new))
            x <- cbind(x, array(sample, dimnames=list(rownames(x), new), dim=c(nrow(x), length(new))))
        if (!isTRUE(all.equal(colnames(x), dn[[2]])))
            x <- x[,dn[[2]],drop=FALSE]
    }
    ii <- cbind(rep(match(rownames(data), rownames(x)), ncol(data)),
                rep(match(colnames(data), colnames(x)), each=nrow(data)))
    x[ii] <- as.vector(data)
    assign(va.name, value=x, pos=1)
}
