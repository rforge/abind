vstorage.mode <- function(x, deep=FALSE) UseMethod('vstorage.mode')
vstorage.mode.varray <- function(x, deep=FALSE) {
    if (deep) {
        return(storage.mode(sapply(x$info, function(y) {
            env <- as.environment(non.null(y$env.name, non.null(x$env.name, 1)))
            comp.data <- non.null(y$value, get(y$name, envir=env, inherits=FALSE))
            vector(storage.mode(comp.data), 1)
        }, simplify=TRUE)))
    } else {
        return(storage.mode(sapply(x$info, '[[', 'sample', simplify=TRUE)))
    }
}
`vstorage.mode<-` <- function(x, value) UseMethod('vstorage.mode<-')
`vstorage.mode<-.varray` <- function(x, value) {
    # should use 'x$umode' here?
    x$info <- lapply(x$info, function(y) {
        if (!is.null(y$value)) {
            if (storage.mode(y$value) != value)
                storage.mode(y$value) <- value
        } else {
            env <- as.environment(non.null(y$env.name, non.null(x$env.name, 1)))
            comp.data <- get(y$name, envir=env, inherits=FALSE)
            if (storage.mode(comp.data) != value) {
                storage.mode(comp.data) <- value
                assign(y$name, envir=env, value=comp.data, inherits=FALSE)
            }
        }
        storage.mode(y$sample) <- value
        y
    })
    x
}
