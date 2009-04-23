nbTs <- function(x, pairwise.deletion = FALSE, as.matrix = FALSE,
                 scaled = FALSE) {
    if(is.list(x)) x <- as.matrix(x)
    nms <- dimnames(x)[[1]]
    n <- dim(x)
    s <- n[2]
    n <- n[1]

    if (!pairwise.deletion) {
        keep <- .C("GlobalDeletionDNA", x, n, s,
                   rep(1L, s), PACKAGE = "ape")[[4]]
        x <- x[,  as.logical(keep)]
        s <- dim(x)[2]
    }

     Ndist <- n*(n - 1)/2

     d <- .C("nb_ts", x, n, s, double(Ndist),
             as.integer(pairwise.deletion),
             as.integer(scaled),
             DUP = FALSE, NAOK = TRUE,
             PACKAGE = "pan")

    d <- d[[4]]
    attr(d, "Size") <- n
    attr(d, "Labels") <- nms
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    if(as.matrix) d <- as.matrix
    d
}
