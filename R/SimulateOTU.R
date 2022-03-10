SimulateOTU <- function(data, nSam=100, parameters, mu = 1000, size = 25) {

    otu.ids <- colnames(data)
    p.est = parameters$pi
    names(p.est) <- names(parameters$pi)
    theta <- parameters$theta
    gplus <- (1 - theta) / theta
    p.est <- p.est[otu.ids]
    g.est <- p.est * gplus
    comm <- matrix(0, nSam, length(g.est))
    rownames(comm) <- 1:nrow(comm)
    colnames(comm) <- names(g.est)
    comm.p <- comm
    nSeq <- rnbinom(nSam, mu = mu, size = size)

    for (i in 1:nSam) {
        comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
        comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
    }
    OTU = comm[, otu.ids]
    return(OTU = OTU)
}
