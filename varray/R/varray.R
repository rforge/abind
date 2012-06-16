varray <- function(..., along=1, dimorder=NULL, env.name=TRUE, envir=NULL, naidxok=NA, dimnames=NULL) {
    # Can call like this:
    #    varray(a, b, c)
    # or varray('a', 'b', 'c')
    # or varray(c('a', 'b', 'c'))
    # Still need to figure out if need to record where to find the objects (i.e., which environments)
    dotv <- list(...)
    fixGlobalEnvName <- function(name) {
        if (name=='R_GlobalEnv') '.GlobalEnv'
        else if (name=='') stop('cannot use an unnamed environment')
        else name
    }
    if (is.character(envir))
        envir <- as.environment(envir)
    if (length(dotv)==1 && is.character(dotv[[1]])) {
        # single argument that is a character vector of names of objects
        info <- lapply(dotv[[1]], FUN=function(arg) {
            name <- arg
            if (!is.null(envir)) {
                x <- get(name, envir=envir)
                x.env.name <- fixGlobalEnvName(environmentName(envir))
            } else {
                x.env.name <- fixGlobalEnvName(find(name))
                x <- get(name)
            }
            if (!env.name)
                x.env.name <- NULL
            if (length(dim(x))==0)
                stop("object '", name, "' is not an array (has no dim)")
            if (is.null(dimnames(x)))
                stop('argument ', name, ' does not refer to an array')
            sample <- asub(x, rep(list(1), length(dim(x))))
            list(name=name, dim=dim(x), dimnames=dimnames(x),
                 env.name=x.env.name, sample=sample,
                 naidxok=if (is.na(naidxok)) inherits(x, 'ff') else naidxok)
        })
    } else {
        # multiple arguments, could be names or objects
        dotargs <- match.call(expand.dots=FALSE)$...
        info <- mapply(list(...), dotargs, SIMPLIFY=FALSE, FUN=function(x, arg) {
            if (is.name(arg)) {
                name <- as.character(arg)
                x.env.name <- fixGlobalEnvName(find(name))
            } else if (is.character(arg) && length(arg)==1) {
                name <- arg
                if (!is.null(envir)) {
                    x <- get(name, envir=envir)
                    x.env.name <- fixGlobalEnvName(environmentName(envir))
                } else {
                    x.env.name <- fixGlobalEnvName(find(name))
                    x <- get(name)
                }
            } else {
                # could consider storing value inside the varray
                stop('arguments must be character data or names')
            }
            if (!env.name)
                x.env.name <- NULL
            if (length(dim(x))==0)
                stop("object '", name, "' is not an array (has no dim)")
            if (is.null(dimnames(x)))
                stop('argument ', name, ' does not refer to an array')
            sample <- asub(x, rep(list(1), length(dim(x))))
            list(name=name, dim=dim(x), dimnames=dimnames(x),
                 env.name=x.env.name, sample=sample,
                 naidxok=if (is.na(naidxok)) inherits(x, 'ff') else naidxok)
        })
    }
    if (along < 1 || along > length(info[[1]]$dim))
        stop('along must be in 1..', length(info[[1]]$dim))
    for (i in seq(along=info)[-1]) {
        if (length(info[[i]]$dim) != length(info[[1]]$dim))
            stop('chunk ', i, ' has different dimensionality: ',
                 length(info[[i]]$dim), ' vs ', length(info[[1]]$dim))
    }
    dn <- lapply(seq(len=length(info[[1]]$dim)), function(i)
                 unique(unlist(lapply(info, function(x) x$dimnames[[i]]))))
    if (!is.null(dimnames))
        for (in in seq(along=dn))
            if (!is.null(dimnames[[i]]))
                dn[[i]] <- dimnames[[i]]
    d <- sapply(dn, length)
    along.idx <- integer(d[along])
    naidxok <- all(sapply(info, '[[', 'naidxok'))
    if (is.null(dimorder))
        dimorder <- seq(length(d))
    if (!identical(sort(as.numeric(dimorder)), as.numeric(seq(length(d)))))
        stop('dimorder must be 1:length(d) in some permutation')
    if (!all(dimorder == seq(length(d)))) {
        d <- d[dimorder]
        dn <- dn[dimorder]
    }
    for (i in seq(to=1, from=length(info))) {
        along.idx[match(info[[i]]$dimnames[[along]], dn[[along]])] <- i
        info[[i]]$map <- lapply(seq(along=dn), function(j) match(dn[[j]], info[[i]]$dimnames[[j]]))
    }
    structure(list(dim=d, dimnames=dn, along=along, info=info,
                   along.idx=along.idx, dimorder=dimorder, naidxok=naidxok),
              class='varray')
}

as.array.varray <- function(x, ...) {
    y <- abind(along=x$along, conform(lapply(x$info,
               function(comp) {
                   if (!is.null(comp$value))
                       comp$value
                   else if (is.null(comp$env.name))
                       get(comp$name, pos=1)
                   else
                       get(comp$name, envir=as.environment(comp$env.name))
               }),
               along=seq(len=length(x$dim))[-x$along]))
    if (!all(x$dimorder == seq(length(x$dim))))
        y <- aperm(y, order(x$dimorder))
    y
}

print.varray <- function(x, ...) {
    dot3 <- function(n) if (n<=4) seq(len=n) else c(1,2,NA,n)
    dnsum <- function(i, d, dn) {
        # summarize dimnames
        if (is.na(i)) {
            return("...")
        } else if (is.null(dn[[i]])) {
            s <- "NULL"
        } else {
            j <- dn[[i]][dot3(d[i])]
            if (any(is.na(j)))
                j[is.na(j)] <- '...'
            s <- paste(j, collapse=', ')
        }
        return(paste('[', i, '] ', s, sep=''))
    }
    cat(paste(x$dim, collapse=' x '), 'virtual', if (length(x$dim)==2) 'matrix' else 'array')
    if (is.null(x$dimnames)) {
        cat(' with NULL dimnames\n')
    } else {
        cat(' with dimnames:\n')
        cat(paste('  ', sapply(dot3(length(x$dim)), dnsum, d=x$dim, dn=x$dimnames), '\n'), sep='')
    }
    if (!all(x$dimorder == seq(length(x$dim))))
        cat('using dims in order', paste(x$dimorder, collapse=', '), '\n')
    cat('with', length(x$info), 'sub-arrays:\n')
    for (k in dot3(length(x$info))) {
        if (is.na(k)) {
            cat('...\n')
        } else {
            cat('sub-array ', k, ': ', sep='')
            if (!is.null(x$info[[k]]$name)) {
                cat("'", x$info[[k]]$name, "' ", sep='')
                if (is.null(x$info[[k]]$env.name)) {
                    found.env <- find(x$info[[k]]$name)
                    if (length(found.env)==0)
                        found.env <- '(not found in any env) '
                    else
                        found.env <- paste("(found in env '", found.env, "') ", sep='')
                } else {
                    env.actual <- try(as.environment(x$info[[k]]$env.name))
                    if (is(env.actual, 'try-error'))
                        found.env <- paste("(no env named '", x$info[[k]]$env.name, "') ", sep='')
                    else if (exists(x$info[[k]]$name, envir=env.actual, inherits=FALSE))
                        found.env <- paste("(in env '", x$info[[k]]$env.name, "') ", sep='')
                    else
                        found.env <- paste("(not found in env '", x$info[[k]]$env.name, "') ", sep='')
                }
                cat(found.env)
            } else {
                cat("(internal value) ")
            }
            cat(paste(x$info[[k]]$dim, collapse=' x '),
                if (length(x$info[[k]]$dim)==2) 'matrix' else 'array')
            if (is.null(x$info[[k]]$dimnames)) {
                cat(' with NULL dimnames\n')
            } else {
                cat(' with dimnames:\n')
                cat(paste('  ', sapply(dot3(length(x$info[[k]]$dim)),
                                       dnsum, d=x$info[[k]]$dim,
                                       dn=x$info[[k]]$dimnames), '\n'), sep='')
            }
        }
    }
}

length.varray <- function(x) prod(x$dim)
dim.varray <- function(x) x$dim
dimnames.varray <- function(x) x$dimnames

mode.varray <- function(x) mode(sapply(x$info, '[[', 'sample'))
storage.mode.varray <- function(x) storage.mode(sapply(x$info, '[[', 'sample'))

"[.varray" <- function(x, ..., drop=TRUE) {
    Nidxs <- nargs() - 1 - (!missing(drop))
    d <- x$dim
    dn <- x$dimnames
    naidxok <- x$naidxok
    if (is.null(naidxok))
        naidxok <- FALSE
    so <- x$dimorder
    if (is.null(so))
        so <- x$storage.order
    # permutation on the stored data to get back to user view
    ap <- order(so)
    # adim and adimnames are the dim and dimnames in the storage order
    ad <- d[so]
    adn <- dn[so]
    integer.max <- .Machine$integer.max
    if (Nidxs==length(d)) {
        # regular indexing
        dotargs <- match.call(expand.dots=FALSE)$...
        # In R, list(...) will stop with an error with empty indices, e.g., x[1,]
        idx.missing <- sapply(dotargs, function(a) is.name(a) && as.character(a)=="")
        ai <- vector("list", length(idx.missing))
        ai[!idx.missing] <- lapply(dotargs[!idx.missing], eval, sys.parent())
        if (length(ai)!=length(d))
            stop("strange ... thought I had ", length(d), " index args, but don't???")
        # put the indices in storage order and convert to integers
        ai <- ai[so]
        idx.missing <- idx.missing[so]
        for (j in seq(along=ai)) {
            if (idx.missing[j]) {
                ai[[j]] <- seq(len=ad[j])
            } else if (is.null(ai[[j]])) {
                # NULL is like numeric(0) in an index
                ai[[j]] <- numeric(0)
            } else if (is.numeric(ai[[j]])) {
                if (storage.mode(ai[[j]])!="integer")
                    k <- as.integer(ai[[j]])
                else
                    k <- ai[[j]]
                k <- k[k!=0]
                if (any(is.true(abs(k) > ad[j])))
                    stop("index out of range on dimension ", ap[j], ", e.g.: ", k[abs(k) > ad[j]][1])
                if (any(is.true(k < 0))) {
                    if (any(is.true(k > 0)))
                        stop("cannot mix negative and positive indices (on dimension ", ap[j], ")")
                    k <- seq(len=ad[j])[k]
                }
                ai[[j]] <- k
            } else if (is.logical(ai[[j]])) {
                if (length(ai[[j]])) {
                    if (ad[j] %% length(ai[[j]]) !=0)
                        stop("dim[j] not a multiple of length(i) for logical index i (on dimension ", ap[j], ")")
                    # need to be careful to preserve NA's in ..1
                    # and create a vector not longer than necessary
                    ai[[j]] <- as.vector(outer(seq(along=ai[[j]])[ai[[j]]],
                                               seq(0,len=ad[j]/length(ai[[j]]), by=length(ai[[j]])), "+"))
                } else {
                    ai[[j]] <- integer(0)
                }
            } else if (is.character(ai[[j]])) {
                if (is.null(adn))
                    stop("cannot use character indexing -- object has no dimension names")
                i <- match(ai[[j]], adn[[j]])
                if (any(is.na(i) & !is.na(ai[[j]])))
                    stop("character indices not in dimnames[[", ap[j], "]], e.g.: ", ai[[j]][is.na(i)][1])
                ai[[j]] <- i
            } else {
                stop("cannot use an index of class ", class(ai[[j]]))
            }
        }
        # Now we have all numeric indices
        # Names of subidx are idx's of submatrix, and values are the indices we want
        subidx <- tapply(seq(along=ai[[x$along]]),
                         factor(x$along.idx[ai[[x$along]]], levels=seq(along=x$info)),
                         FUN=c, simplify=FALSE)
        a <- lapply(which(sapply(subidx, length)>0), function(i) {
            # i is the index of the submatrix we need some data from
            ii <- replace(ai, x$along, list(ai[[x$along]][subidx[[i]]]))
            # map to indices for this submatrix
            jj <- lapply(seq(len=length(x$dim)), function(j)
                return(x$info[[i]]$map[[j]][ii[[j]]]))
            has.nas <- FALSE
            if (!naidxok) {
                # jj1 is the indices without any NA's
                # jj2 is how to find the indices including NA's in jj1
                jj1 <- jj
                jj2 <- jj
                for (k in seq(along=jj)) {
                    if (any(jj.na <- is.na(jj[[k]]))) {
                        jj1[[k]] <- jj[[k]][!jj.na]
                        has.nas <- TRUE
                    }
                    jj2[[k]] <- match(jj[[k]], jj1[[k]])
                }
            }
            if (naidxok || !has.nas) {
                if (!is.null(x$info[[i]]$value)) {
                    y <- do.call('[', c(list(x$info[[i]]$value), jj, list(drop=FALSE)))
                } else {
                    if (is.null(x$info[[i]]$env.name))
                        yy <- get(x$info[[i]]$name, pos=1)
                    else
                        yy <- get(x$info[[i]]$name, envir=as.environment(x$info[[i]]$env.name))
                    y <- do.call('[', c(list(yy), jj, list(drop=FALSE)))
                }
            } else {
                if (!is.null(x$info[[i]]$value)) {
                    y1 <- do.call('[', c(list(x$info[[i]]$value), jj1, list(drop=FALSE)))
                } else {
                    if (is.null(x$info[[i]]$env.name))
                        yy <- get(x$info[[i]]$name, pos=1)
                    else
                        yy <- get(x$info[[i]]$name, envir=as.environment(x$info[[i]]$env.name))
                    y1 <- do.call('[', c(list(yy), jj1, list(drop=FALSE)))
                }
                # we assume that the object that is result of indexing can handle NA's in indices
                # (ordinary R arrays are fine with NA's in indices, it's ff arrays that are not)
                y <- do.call('[', c(list(y1), jj2, list(drop=FALSE)))
            }
            dimnames(y) <- NULL
            y
        })
        if (length(a) > 1) {
            a <- abind(along=x$along, a)
            # are subchunks out-of-order?
            # get numeric indices of the along-dim as they appear in a
            b <- unlist(lapply(which(sapply(subidx, length)>0), function(i)
                ai[[x$along]][subidx[[i]]]), use.names=FALSE)
            # if there are NA's in ai[[x$along]], it's length will not equal b's
            if (any(is.na(ai[[x$along]])) != (length(b)!=length(ai[[x$along]])))
                stop("internal indexing inconsistency: expecting NA's in ai[[x$along]] iff length(b)!=length(ai[[x$along]])")
            if (length(b)!=length(ai[[x$along]]) || !all(b==ai[[x$along]]))
                a <- asub(a, match(ai[[x$along]], b), x$along)
        } else if (length(a) == 1) {
            a <- a[[1]]
            if (any(is.na(ai[[x$along]]))) {
                i <- which(sapply(subidx, length)>0)
                if (length(i)!=1)
                    stop("internal indexing inconsistency: expecting only one non-null subidx")
                b <- ai[[x$along]][subidx[[i]]]
                a <- asub(a, match(ai[[x$along]], b), x$along)
            }
        } else {
            a <- array(numeric(0), dim=sapply(ai, length))
        }
        if (!is.null(dn))
            dimnames(a) <- lapply(seq(len=length(d)), function(i) dn[[i]][ai[[i]]])
        if (!identical(drop, FALSE)) {
            if (identical(drop, TRUE))
                drop <- which(sapply(ai, length)==1)
            if (length(drop))
                a <- adrop(a, drop=drop)
        }
    } else if (Nidxs==1 && !missing(..1)) {
        # matrix or vector indexing
        if (!is.matrix(..1) && (is.logical(..1) || is.numeric(..1))) {
            if (is.logical(..1)) {
                if (! (length(..1) %in% cumprod(c(1,x$dim))))
                    stop("length(x) is not the length of a sub-array for logical index i")
                # need to be careful to preserve NA's in ..1
                # and create a vector not longer than necessary
                vi <- as.vector(outer(seq(along=..1)[..1],
                                      seq(0,len=prod(x$dim)/length(..1), by=length(..1)), "+"))
            } else if (is.numeric(..1)) {
                # numeric
                vi <- ..1
                if (any(abs(vi) > x$info$length, na.rm=TRUE))
                    stop("vector indices out of range")
                if (any(vi < 0, na.rm=TRUE)) {
                    # using -ve vector indices results in the allocation of a
                    # vector indexing the entire array here!
                    if (any(is.na(vi) | vi>0))
                        stop("Only 0's may be mixed with -ve indices")
                    vi <- seq(x$info$length)[vi]
                } else if (any(j <- is.true(vi==0))) {
                    vi <- vi[!j]
                }
            }
            mi <- matrix(0L, nrow=length(vi), ncol=length(d))
            ii <- vi-1
            for (j in seq(along=d)) {
                mi[,j] <- as.integer((ii %% d[j]) + 1)
                if (j < length(d))
                    ii <- floor(ii / d[j])
            }
        } else if (is.matrix(..1)) {
            if (ncol(..1)!=length(d))
                stop("a single argument must be a ", length(d), " column matrix")
            if (mode(..1)!="numeric")
                stop("matrix indexing only works with numeric matrices")
            mi <- ..1[,so]
            vi <- NULL
        } else {
            stop("a single argument must be a ", length(d), " column matrix or a vector")
        }
        # now we have 'mi' - a matrix index
        if (any(i <- (rowSums(mi==0, na.rm=TRUE) > 0)))
            mi <- mi[!i,,drop=FALSE]
        if (storage.mode(mi)!="integer")
            storage.mode(mi) <- "integer"
        if (is.true(any(mi < 1 | (mi - matrix(ad, nrow=nrow(mi), ncol=ncol(mi), byrow=TRUE)) > 0)))
            stop("matrix indices out of range")
        subidx <- tapply(seq(along=mi[,x$along]),
                         factor(x$along.idx[mi[,x$along]], levels=seq(along=x$info)),
                         FUN=c, simplify=FALSE)
        a <- lapply(which(sapply(subidx, length)>0), function(i) {
            # i is the index of the submatrix we need some data from
            # jj is the matrix indices just for this submatrix
            jj <- mi[subidx[[i]],,drop=FALSE]
            # map to indices for this submatrix
            for (j in seq(len=length(x$dim)))
                jj[,j] <- x$info[[i]]$map[[j]][jj[,j]]
            if (!is.null(x$info[[i]]$value)) {
                y <- x$info[[i]]$value[jj]
            } else {
                if (is.null(x$info[[i]]$name))
                    yy <- get(x$info[[i]]$name, pos=1)
                else
                    yy <- get(x$info[[i]]$name, envir=as.environment(x$info[[i]]$env.name))
                y <- yy[jj]
            }
            dimnames(y) <- NULL
            y
        })
        if (length(a) > 1) {
            a <- unlist(a, use.names=FALSE)
            # are subchunks out-of-order?
            # get numeric indices of the along-dim as they appear in 'a'
            b <- unlist(lapply(which(sapply(subidx, length)>0), function(i)
                mi[subidx[[i]],x$along,drop=FALSE]), use.names=FALSE)
            # if there are NA's in mi[,x$along], it's length will not equal b's
            if (any(is.na(mi[,x$along])) != (length(b)!=nrow(mi)))
                stop("internal indexing inconsistency: expecting NA's in mi[,x$along] iff length(b)!=nrow(mi)")
            if (length(b)!=length(mi[,x$along]) || !all(b==mi[,x$along])) {
                # have to reconstruct the full numerical vector index
                if (is.null(vi)) {
                    vi <- mi[,1]
                    for (j in seq(ncol(mi))[-1])
                        vi <- vi + ((mi[,j]-1) * prod(c(1,d)[seq(len=j)]))
                }
                ai <- vi[unlist(subidx)]
                a <- a[match(vi, ai)]
            }
        } else if (length(a) == 1) {
            if (any(is.na(mi[,x$along]))) {
                i <- which(sapply(subidx, length)>0)
                if (length(i)!=1)
                    stop("internal indexing inconsistency: expecting only one non-null subidx")
                a <- replace(rep(NA, nrow(mi)), which(!is.na(mi[,x$along])), a[[1]])
            } else {
                a <- a[[1]]
            }
        } else {
            a <- array(numeric(0), dim=sapply(ai, length))
        }
    } else if (Nidxs<=0) {
        # get both [] and [,drop=F]
        stop("empty subscripts not supported for 'varray' objects")
    } else {
        stop("need ", length(d), " or 1 (matrix or vector) index arguments")
    }
    if (!all(x$dimorder == seq(length(x$dim))))
        return(aperm(a, order(x$dimorder)))
    else
        return(a)
}

is.true <- function(x) (x & !is.na(x))

"dimnames<-.varray" <- function(x, value) stop('dimnames for varray are read-only')
"dim<-.varray" <- function(x, value) stop('dim for varray is read-only')
"length<-.varray" <- function(x, value) stop('length for varray is read-only')
"mode<-.varray" <- function(x, value) stop('mode for varray is read-only')
"storage.mode<-.varray" <- function(x, value) stop('storage.mode for varray is read-only')
"[<-.varray" <- function(x, i, j, ..., value) stop('cannot replace parts a varray (varray is read-only -- you must work with the sub-arrays)')

