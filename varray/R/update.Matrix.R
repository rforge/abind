update.Matrix <- function(x.name, data, need.dimnames=list(NULL, NULL)) {
    if (!exists(x.name)) {
        x <- Matrix(data, sparse=TRUE)
    } else {
        x <- get(x.name)
    }
    dn <- list(sort(unique(c(rownames(x), rownames(data), need.dimnames[[1]])), na.last=NA),
               sort(unique(c(colnames(x), colnames(data), need.dimnames[[2]])), na.last=NA))
    sample <- replace(x[1,1], 1, 0L)
    if (!isTRUE(all.equal(rownames(x), dn[[1]]))) {
        new <- setdiff(dn[[1]], rownames(x))
        if (length(new))
            x <- rBind(x, array2(sample, dimnames=list(new, colnames(x))))
        if (!isTRUE(all.equal(rownames(x), dn[[1]])))
            x <- x[dn[[1]],,drop=FALSE]
    }
    if (!isTRUE(all.equal(colnames(x), dn[[2]]))) {
        new <- setdiff(dn[[2]], colnames(x))
        if (length(new))
            x <- cBind(x, array2(sample, dimnames=list(rownames(x), new)))
        if (!isTRUE(all.equal(rownames(x), dn[[1]])))
            x <- x[dn[[1]],,drop=FALSE]
    }
    ii <- cbind(rep(match(rownames(data), rownames(x)), ncol(data)),
                rep(match(colnames(data), colnames(x)), each=nrow(data)))
    x[ii] <- as.vector(data)
    assign(x.name, value=x, pos=1)
}
