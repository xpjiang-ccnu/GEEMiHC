GEEMiHC <-
  function(y, id, covs = NULL, otu.tab, tree = NULL, model, Gamma = c(1,3,5,7,9), Lamda = matrix(c(1, rep(0, 8), rep(1/3, 3), rep(0, 6), rep(1/5, 5), rep(0, 4), rep(1/7, 7), rep(0, 2), rep(1/9, 9)), 5, 9, byrow = T), 
           comp = FALSE, CLR = FALSE, opt.ncl = 30, n.perm = 5000) {
    if (length(Gamma) != nrow(Lamda) | max(Gamma) != ncol(Lamda) | length(Gamma) > 100 | length(Gamma) < 1) {
      stop("The number(maximum value) of the Gamma is out of range or not equal to the number of rows(columns) in Lamda. Please enter the correct Gamma or Lamda!")
    }
    
    for (i in 1:nrow(Lamda)) {
      if (sum(Lamda[i,]) != 1){
        stop(paste(paste("The sum of Lamda's", i, sep = " "), paste("th row isn't 1. Please enter the correct Lamda!"), sep = ""))
      }
      if (i < nrow(Lamda) && max(Lamda[i, as.numeric((Gamma[i]+1):max(Gamma))]) != 0){
        stop("Non-zero values occupy the position of the 0 element in Lamda. Please enter the correct Lamda!")
      }
    }
    
    if (!comp) {
      otu.tab <- t(apply(otu.tab, 1, function(x)x/sum(x)))
    } 
    com.otu.tab <- otu.tab
    if (CLR) {
      otu.tab <- t(apply(com.otu.tab, 1, clr))
    } 
    n <- length(y)
    p <- ncol(otu.tab)
    options(warn=-1)
    
    if (is.null(covs)) {
      fit.GEE.AR <- MGEE(y ~ 1, id=id , family = model, corstr = "AR-1")
      f.y.GEE.AR <- fitted.values(fit.GEE.AR)
      r.GEE.AR <- y - f.y.GEE.AR
      
      fit.GEE.EX <- MGEE(y ~ 1, id=id , family = model, corstr = "exchangeable")
      f.y.GEE.EX <- fitted.values(fit.GEE.EX)
      r.GEE.EX <- y - f.y.GEE.EX
      
      fit.GEE.IN <- MGEE(y ~ 1, id=id , family = model, corstr = "independence")
      f.y.GEE.IN <- fitted.values(fit.GEE.IN)
      r.GEE.IN <- y - f.y.GEE.IN
      
    } else {
      fit.GEE.AR <- MGEE(y ~ ., id=id , family = model, corstr = "AR-1", as.data.frame(covs))
      f.y.GEE.AR <- fitted.values(fit.GEE.AR)
      r.GEE.AR <- y - f.y.GEE.AR
      
      fit.GEE.EX <- MGEE(y ~ ., id=id , family = model, corstr = "exchangeable", as.data.frame(covs))
      f.y.GEE.EX <- fitted.values(fit.GEE.EX)
      r.GEE.EX <- y - f.y.GEE.EX
      
      fit.GEE.IN <- MGEE(y ~ ., id=id , family = model, corstr = "independence", as.data.frame(covs))
      f.y.GEE.IN <- fitted.values(fit.GEE.IN)
      r.GEE.IN <- y - f.y.GEE.IN
      
    }
    
    r.s.GEE.AR <- list()
    r.s.GEE.EX <- list()
    r.s.GEE.IN <- list()
    for (j in 1:n.perm) {
    ind_r.s <- shuffle(n)
    r.s.GEE.AR[[j]] <- r.GEE.AR[ind_r.s]
    r.s.GEE.EX[[j]] <- r.GEE.EX[ind_r.s]
    r.s.GEE.IN[[j]] <- r.GEE.IN[ind_r.s]
    }
    
    n.Zs.GEE.AR <- as.numeric(r.GEE.AR) %*% otu.tab
    n.Zs.GEE.EX <- as.numeric(r.GEE.EX) %*% otu.tab
    n.Zs.GEE.IN <- as.numeric(r.GEE.IN) %*% otu.tab
    n.Z0s.GEE.AR <- lapply(r.s.GEE.AR, function(x) as.numeric(x %*% otu.tab))
    n.Z0s.GEE.EX <- lapply(r.s.GEE.EX, function(x) as.numeric(x %*% otu.tab))
    n.Z0s.GEE.IN <- lapply(r.s.GEE.IN, function(x) as.numeric(x %*% otu.tab))
    U0s.GEE.AR <- lapply(apply(sapply(r.s.GEE.AR, function(x) return(x %*% otu.tab)),1,list),unlist)    
    U0s.GEE.EX <- lapply(apply(sapply(r.s.GEE.EX, function(x) return(x %*% otu.tab)),1,list),unlist)
    U0s.GEE.IN <- lapply(apply(sapply(r.s.GEE.IN, function(x) return(x %*% otu.tab)),1,list),unlist)
    se.GEE.AR <- sapply(U0s.GEE.AR, sd)
    se.GEE.EX <- sapply(U0s.GEE.EX, sd)
    se.GEE.IN <- sapply(U0s.GEE.IN, sd)
    Zs.GEE.AR <- n.Zs.GEE.AR/se.GEE.AR
    Zs.GEE.EX <- n.Zs.GEE.EX/se.GEE.EX
    Zs.GEE.IN <- n.Zs.GEE.IN/se.GEE.IN
    Z0s.GEE.AR <- lapply(n.Z0s.GEE.AR, function(x) x/se.GEE.AR)
    Z0s.GEE.EX <- lapply(n.Z0s.GEE.EX, function(x) x/se.GEE.EX)
    Z0s.GEE.IN <- lapply(n.Z0s.GEE.IN, function(x) x/se.GEE.IN)

      if (class(tree) == "phylo") {
        CDis <- cophenetic(tree)	
        for (j in 1:nrow(CDis)) {
          ind.t <- which(CDis[j,] == 0)
          ind.d <- which(ind.t == j)
          ind.r <- ind.t[-ind.d]
          CDis[j,ind.r] <- min(CDis[j,-ind.t])/2
        }	
        
        if (opt.ncl == "gmx") {
          asw <- numeric(p-1)
          for (j in 2:(p-1)) {
            asw[j] <- pam(CDis, j)$silinfo$avg.width
          }
          k.best <- which.max(asw)
          clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
        }
        if (opt.ncl == "fmx") {
          asw <- numeric(p-1)
          j <- 2
          asw[j] <- pam(CDis, j)$silinfo$avg.width
          while (asw[j-1] < asw[j] & j <= p-1) {
            j <- j + 1
            asw[j] <- pam(CDis, j)$silinfo$avg.width
          }		
          k.best <- j-1
          clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
        }
        if (opt.ncl != "gmx" & opt.ncl != "fmx") {
          asw <- numeric(p-1)
          for (j in 2:opt.ncl) {
            asw[j] <- pam(CDis, j)$silinfo$avg.width
          }
          k.best <- which.max(asw)
          clust.pam <- pam(CDis, k.best, cluster.only=TRUE)
        }
        Ws.GEE.AR <- rep(NA, p)
        Ws.GEE.EX <- rep(NA, p)
        Ws.GEE.IN <- rep(NA, p)
        for (j in 1:p) {
          ind.1 <- match(colnames(CDis)[j], names(clust.pam))
          ind.2 <- which(clust.pam == as.numeric(clust.pam[ind.1]))
          ind.3 <- which(ind.2 == ind.1)
          ind <- ind.2[-ind.3]
          if (length(ind) == 0) {
            Ws.GEE.AR[j] <- 1
            Ws.GEE.EX[j] <- 1
            Ws.GEE.IN[j] <- 1
          } else {
            inv.Ds <- 1/(CDis[j,ind])
            abs.Zs.GEE.AR <- abs(Zs.GEE.AR[ind])
            abs.Zs.GEE.EX <- abs(Zs.GEE.EX[ind])
            abs.Zs.GEE.IN <- abs(Zs.GEE.IN[ind])
            Ws.GEE.AR[j] <- sum(inv.Ds*abs.Zs.GEE.AR)/sum(inv.Ds) + 1
            Ws.GEE.EX[j] <- sum(inv.Ds*abs.Zs.GEE.EX)/sum(inv.Ds) + 1
            Ws.GEE.IN[j] <- sum(inv.Ds*abs.Zs.GEE.IN)/sum(inv.Ds) + 1
          }
        }
        Ws.GEE.AR <- Ws.GEE.AR/sum(Ws.GEE.AR)
        Ws.GEE.EX <- Ws.GEE.EX/sum(Ws.GEE.EX)
        Ws.GEE.IN <- Ws.GEE.IN/sum(Ws.GEE.IN)
        
        ahc.o.GEE.AR <- as.list(GEEMiHC.stat(Zs=Zs.GEE.AR, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.AR))
        ahc.o.GEE.EX <- as.list(GEEMiHC.stat(Zs=Zs.GEE.EX, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.EX))
        ahc.o.GEE.IN <- as.list(GEEMiHC.stat(Zs=Zs.GEE.IN, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.IN))
        ahc.p.GEE.AR <- lapply(apply(sapply(Z0s.GEE.AR, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.AR)), 1, list), unlist)
        ahc.p.GEE.EX <- lapply(apply(sapply(Z0s.GEE.EX, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.EX)), 1, list), unlist)
        ahc.p.GEE.IN <- lapply(apply(sapply(Z0s.GEE.IN, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=Ws.GEE.IN)), 1, list), unlist)
        pvs.GEE.AR <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.AR, ahc.p.GEE.AR)
        pvs.GEE.EX <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.EX, ahc.p.GEE.EX)
        pvs.GEE.IN <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.IN, ahc.p.GEE.IN)
        
        ind.pvs.GEE.AR <- 1-pchisq(Zs.GEE.AR^2, df=1)
        ind.pvs.GEE.EX <- 1-pchisq(Zs.GEE.EX^2, df=1)
        ind.pvs.GEE.IN <- 1-pchisq(Zs.GEE.IN^2, df=1)
        simes.pv.GEE.AR <- round(min(length(ind.pvs.GEE.AR)*ind.pvs.GEE.AR/rank(ind.pvs.GEE.AR)),4)
        simes.pv.GEE.EX <- round(min(length(ind.pvs.GEE.EX)*ind.pvs.GEE.EX/rank(ind.pvs.GEE.EX)),4)
        simes.pv.GEE.IN <- round(min(length(ind.pvs.GEE.IN)*ind.pvs.GEE.IN/rank(ind.pvs.GEE.IN)),4)
        simes.pv0s.ori.GEE.AR <- unlist(lapply(Z0s.GEE.AR, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.ori.GEE.EX <- unlist(lapply(Z0s.GEE.EX, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.ori.GEE.IN <- unlist(lapply(Z0s.GEE.IN, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.GEE.AR <- rep(NA, n.perm)
        simes.pv0s.GEE.EX <- rep(NA, n.perm)
        simes.pv0s.GEE.IN <- rep(NA, n.perm)
        
        for (j in 1:n.perm) {	
          simes.pv0s.GEE.AR[j] <- (length(which(simes.pv0s.ori.GEE.AR[-j] < simes.pv0s.ori.GEE.AR[j]))+0.01)/(n.perm+0.01)
          simes.pv0s.GEE.EX[j] <- (length(which(simes.pv0s.ori.GEE.EX[-j] < simes.pv0s.ori.GEE.EX[j]))+0.01)/(n.perm+0.01)
          simes.pv0s.GEE.IN[j] <- (length(which(simes.pv0s.ori.GEE.IN[-j] < simes.pv0s.ori.GEE.IN[j]))+0.01)/(n.perm+0.01)
        }

        simes.pv.GEE.AR <- (length(which(simes.pv0s.GEE.AR < simes.pv.GEE.AR))+0.01)/(n.perm+0.01)
        simes.pv.GEE.EX <- (length(which(simes.pv0s.GEE.EX < simes.pv.GEE.EX))+0.01)/(n.perm+0.01)
        simes.pv.GEE.IN <- (length(which(simes.pv0s.GEE.IN < simes.pv.GEE.IN))+0.01)/(n.perm+0.01)
        
        l.Gamma <- length(Gamma)
        Tu.GEE.AR <- min(pvs.GEE.AR[1:l.Gamma], simes.pv.GEE.AR)
        Tu.GEE.EX <- min(pvs.GEE.EX[1:l.Gamma], simes.pv.GEE.EX)
        Tu.GEE.IN <- min(pvs.GEE.IN[1:l.Gamma], simes.pv.GEE.IN)
        T0u.GEE.AR <- rep(NA, n.perm)
        T0u.GEE.EX <- rep(NA, n.perm)
        T0u.GEE.IN <- rep(NA, n.perm)
        Tw.GEE.AR <- min(pvs.GEE.AR[(l.Gamma+1):(l.Gamma*2)], simes.pv.GEE.AR)
        Tw.GEE.EX <- min(pvs.GEE.EX[(l.Gamma+1):(l.Gamma*2)], simes.pv.GEE.EX)
        Tw.GEE.IN <- min(pvs.GEE.IN[(l.Gamma+1):(l.Gamma*2)], simes.pv.GEE.IN)
        T0w.GEE.AR <- rep(NA, n.perm)
        T0w.GEE.EX <- rep(NA, n.perm)
        T0w.GEE.IN <- rep(NA, n.perm)
        for (j in 1:n.perm) {
          T0u.o.GEE.AR <- lapply(ahc.p.GEE.AR[1:l.Gamma], function(x) x[j])
          T0u.o.GEE.EX <- lapply(ahc.p.GEE.EX[1:l.Gamma], function(x) x[j])
          T0u.o.GEE.IN <- lapply(ahc.p.GEE.IN[1:l.Gamma], function(x) x[j])
          T0u.p.GEE.AR <- lapply(ahc.p.GEE.AR[1:l.Gamma], function(x) x[-j])
          T0u.p.GEE.EX <- lapply(ahc.p.GEE.EX[1:l.Gamma], function(x) x[-j])
          T0u.p.GEE.IN <- lapply(ahc.p.GEE.IN[1:l.Gamma], function(x) x[-j])
          T0u.GEE.AR[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.AR, T0u.p.GEE.AR))
          T0u.GEE.EX[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.EX, T0u.p.GEE.EX))
          T0u.GEE.IN[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.IN, T0u.p.GEE.IN))
        }
        for (j in 1:n.perm) {
          T0w.o.GEE.AR <- lapply(ahc.p.GEE.AR[(l.Gamma+1):(l.Gamma*2)], function(x) x[j])
          T0w.o.GEE.EX <- lapply(ahc.p.GEE.EX[(l.Gamma+1):(l.Gamma*2)], function(x) x[j])
          T0w.o.GEE.IN <- lapply(ahc.p.GEE.IN[(l.Gamma+1):(l.Gamma*2)], function(x) x[j])
          T0w.p.GEE.AR <- lapply(ahc.p.GEE.AR[(l.Gamma+1):(l.Gamma*2)], function(x) x[-j])
          T0w.p.GEE.EX <- lapply(ahc.p.GEE.EX[(l.Gamma+1):(l.Gamma*2)], function(x) x[-j])
          T0w.p.GEE.IN <- lapply(ahc.p.GEE.IN[(l.Gamma+1):(l.Gamma*2)], function(x) x[-j])
          T0w.GEE.AR[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0w.o.GEE.AR, T0w.p.GEE.AR))
          T0w.GEE.EX[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0w.o.GEE.EX, T0w.p.GEE.EX))
          T0w.GEE.IN[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0w.o.GEE.IN, T0w.p.GEE.IN))
        }
        
        T.GEE.AR <- min(pvs.GEE.AR, simes.pv.GEE.AR)
        T0.GEE.AR <- apply(cbind(T0u.GEE.AR, T0w.GEE.AR, simes.pv0s.GEE.AR),1,min)
        pv.opt.GEE.AR <- (length(which(T0.GEE.AR < T.GEE.AR))+0.01)/(n.perm+0.01)
        
        T.GEE.EX <- min(pvs.GEE.EX, simes.pv.GEE.EX)
        T0.GEE.EX <- apply(cbind(T0u.GEE.EX, T0w.GEE.EX, simes.pv0s.GEE.EX),1,min)
        pv.opt.GEE.EX <- (length(which(T0.GEE.EX < T.GEE.EX))+0.01)/(n.perm+0.01)
        
        T.GEE.IN <- min(pvs.GEE.IN, simes.pv.GEE.IN)
        T0.GEE.IN <- apply(cbind(T0u.GEE.IN, T0w.GEE.IN, simes.pv0s.GEE.IN),1,min)
        pv.opt.GEE.IN <- (length(which(T0.GEE.IN < T.GEE.IN))+0.01)/(n.perm+0.01)
        
        T <- min(T.GEE.AR, T.GEE.EX, T.GEE.IN)
        T0 <- apply(cbind(T0.GEE.AR, T0.GEE.EX, T0.GEE.IN),1,min)
        pv.opt <- (length(which(T0 < T))+0.01)/(n.perm+0.01)
        
        names(simes.pv.GEE.AR) <- "Simes.GEE(AR)"
        names(simes.pv.GEE.EX) <- "Simes.GEE(EX)"
        names(simes.pv.GEE.IN) <- "Simes.GEE(IN)"
        item.pvs.GEE.AR <- cbind(pvs.GEE.AR[1:l.Gamma], pvs.GEE.AR[(l.Gamma+1):(l.Gamma*2)])
        item.pvs.GEE.EX <- cbind(pvs.GEE.EX[1:l.Gamma], pvs.GEE.EX[(l.Gamma+1):(l.Gamma*2)])
        item.pvs.GEE.IN <- cbind(pvs.GEE.IN[1:l.Gamma], pvs.GEE.IN[(l.Gamma+1):(l.Gamma*2)])
        colnames(item.pvs.GEE.AR) <- c("uGEEHC(h).AR", "wGEEHC(h).AR")
        colnames(item.pvs.GEE.EX) <- c("uGEEHC(h).EX", "wGEEHC(h).EX")
        colnames(item.pvs.GEE.IN) <- c("uGEEHC(h).IN", "wGEEHC(h).IN")
        rownames(item.pvs.GEE.AR) <- as.character(Gamma)
        rownames(item.pvs.GEE.EX) <- as.character(Gamma)
        rownames(item.pvs.GEE.IN) <- as.character(Gamma)
        
        ada.pvs <- c(pv.opt.GEE.AR, pv.opt.GEE.EX, pv.opt.GEE.IN, pv.opt)
        names(ada.pvs) <- c("GEEMiHC(AR)", "GEEMiHC(EX)", "GEEMiHC(IN)", "aGEEMiHC")
        
        otu.ids <- colnames(otu.tab)
        otu.ids <- sub("New.ReferenceOTU", "N", otu.ids)
        otu.ids <- sub("New.CleanUp.ReferenceOTU", "NC", otu.ids)
        pvs.GEE.AR <- ind.pvs.GEE.AR
        pvs.GEE.EX <- ind.pvs.GEE.EX
        pvs.GEE.IN <- ind.pvs.GEE.IN
        i0.GEE.AR <- which(pvs.GEE.AR < 0.00000001)
        i0.GEE.EX <- which(pvs.GEE.EX < 0.00000001)
        i0.GEE.IN <- which(pvs.GEE.IN < 0.00000001)
        i1.GEE.AR <- which(pvs.GEE.AR > 0.99999999)
        i1.GEE.EX <- which(pvs.GEE.EX > 0.99999999)
        i1.GEE.IN <- which(pvs.GEE.IN > 0.99999999)
        pvs.GEE.AR[i0.GEE.AR] <- 0.00000001
        pvs.GEE.EX[i0.GEE.EX] <- 0.00000001
        pvs.GEE.IN[i0.GEE.IN] <- 0.00000001
        pvs.GEE.AR[i1.GEE.AR] <- 0.99999999
        pvs.GEE.EX[i1.GEE.EX] <- 0.99999999
        pvs.GEE.IN[i1.GEE.IN] <- 0.99999999
        Is.GEE.AR <- order(order(pvs.GEE.AR))
        Is.GEE.EX <- order(order(pvs.GEE.EX))
        Is.GEE.IN <- order(order(pvs.GEE.IN))
        exp.HC.GEE.AR <- (Is.GEE.AR/p)/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        exp.HC.GEE.EX <- (Is.GEE.EX/p)/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        exp.HC.GEE.IN <- (Is.GEE.IN/p)/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        obs.HC.GEE.AR <- pvs.GEE.AR/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        obs.HC.GEE.EX <- pvs.GEE.EX/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        obs.HC.GEE.IN <- pvs.GEE.IN/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        exp.wHC.GEE.AR <- Ws.GEE.AR*(Is.GEE.AR/p)/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        exp.wHC.GEE.EX <- Ws.GEE.EX*(Is.GEE.EX/p)/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        exp.wHC.GEE.IN <- Ws.GEE.IN*(Is.GEE.IN/p)/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        obs.wHC.GEE.AR <- Ws.GEE.AR*pvs.GEE.AR/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        obs.wHC.GEE.EX <- Ws.GEE.EX*pvs.GEE.EX/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        obs.wHC.GEE.IN <- Ws.GEE.IN*pvs.GEE.IN/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        deviation.uGEEHC.AR <- abs(exp.HC.GEE.AR - obs.HC.GEE.AR)
        deviation.uGEEHC.EX <- abs(exp.HC.GEE.EX - obs.HC.GEE.EX)
        deviation.uGEEHC.IN <- abs(exp.HC.GEE.IN - obs.HC.GEE.IN)
        deviation.wGEEHC.AR <- abs(exp.wHC.GEE.AR - obs.wHC.GEE.AR)
        deviation.wGEEHC.EX <- abs(exp.wHC.GEE.EX - obs.wHC.GEE.EX)
        deviation.wGEEHC.IN <- abs(exp.wHC.GEE.IN - obs.wHC.GEE.IN)
        rank.uGEEHC.AR <- otu.ids[order(deviation.uGEEHC.AR, decreasing = T)]
        rank.uGEEHC.EX <- otu.ids[order(deviation.uGEEHC.EX, decreasing = T)]
        rank.uGEEHC.IN <- otu.ids[order(deviation.uGEEHC.IN, decreasing = T)]
        rank.wGEEHC.AR <- otu.ids[order(deviation.wGEEHC.AR, decreasing = T)]
        rank.wGEEHC.EX <- otu.ids[order(deviation.wGEEHC.EX, decreasing = T)]
        rank.wGEEHC.IN <- otu.ids[order(deviation.wGEEHC.IN, decreasing = T)]
        return(list(graph.els = list(otu.ids = otu.ids, exp.GEEHC.AR = exp.HC.GEE.AR, exp.GEEHC.EX = exp.HC.GEE.EX, exp.GEEHC.IN = exp.HC.GEE.IN, obs.GEEHC.AR = obs.HC.GEE.AR, obs.GEEHC.EX = obs.HC.GEE.EX, obs.GEEHC.IN = obs.HC.GEE.IN, exp.wGEEHC.AR = exp.wHC.GEE.AR, exp.wGEEHC.EX = exp.wHC.GEE.EX, exp.wGEEHC.IN = exp.wHC.GEE.IN, obs.wGEEHC.AR = obs.wHC.GEE.AR, obs.wGEEHC.EX = obs.wHC.GEE.EX, obs.wGEEHC.IN = obs.wHC.GEE.IN), opt.nclust = k.best, rank.order = list(uGEEHC.AR = rank.uGEEHC.AR, uGEEHC.EX = rank.uGEEHC.EX, uGEEHC.IN = rank.uGEEHC.IN, wGEEHC.AR = rank.wGEEHC.AR, wGEEHC.EX = rank.wGEEHC.EX, wGEEHC.IN = rank.wGEEHC.IN), simes.pv.GEE.AR = simes.pv.GEE.AR, simes.pv.GEE.EX = simes.pv.GEE.EX, simes.pv.GEE.IN = simes.pv.GEE.IN, ind.pvs.GEEMiHC.AR = item.pvs.GEE.AR, ind.pvs.GEEMiHC.EX = item.pvs.GEE.EX, ind.pvs.GEEMiHC.IN = item.pvs.GEE.IN, ada.pvs=ada.pvs))
      }
      
      if (class(tree) != "phylo") {
        ahc.o.GEE.AR <- as.list(GEEMiHC.stat(Zs=Zs.GEE.AR, Gamma=Gamma, Lamda=Lamda, Ws=NULL))
        ahc.o.GEE.EX <- as.list(GEEMiHC.stat(Zs=Zs.GEE.EX, Gamma=Gamma, Lamda=Lamda, Ws=NULL))
        ahc.o.GEE.IN <- as.list(GEEMiHC.stat(Zs=Zs.GEE.IN, Gamma=Gamma, Lamda=Lamda, Ws=NULL))
        ahc.p.GEE.AR <- lapply(apply(sapply(Z0s.GEE.AR, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=NULL)), 1, list), unlist)
        ahc.p.GEE.EX <- lapply(apply(sapply(Z0s.GEE.EX, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=NULL)), 1, list), unlist)
        ahc.p.GEE.IN <- lapply(apply(sapply(Z0s.GEE.IN, function(x) GEEMiHC.stat(Zs=x, Gamma=Gamma, Lamda=Lamda, Ws=NULL)), 1, list), unlist)
        pvs.GEE.AR <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.AR, ahc.p.GEE.AR)
        pvs.GEE.EX <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.EX, ahc.p.GEE.EX)
        pvs.GEE.IN <- mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), ahc.o.GEE.IN, ahc.p.GEE.IN)
        
        ind.pvs.GEE.AR <- 1-pchisq(Zs.GEE.AR^2, df=1)
        ind.pvs.GEE.EX <- 1-pchisq(Zs.GEE.EX^2, df=1)
        ind.pvs.GEE.IN <- 1-pchisq(Zs.GEE.IN^2, df=1)
        simes.pv.GEE.AR <- min(length(ind.pvs.GEE.AR)*ind.pvs.GEE.AR/rank(ind.pvs.GEE.AR))
        simes.pv.GEE.EX <- min(length(ind.pvs.GEE.EX)*ind.pvs.GEE.EX/rank(ind.pvs.GEE.EX))
        simes.pv.GEE.IN <- min(length(ind.pvs.GEE.IN)*ind.pvs.GEE.IN/rank(ind.pvs.GEE.IN))
        simes.pv0s.ori.GEE.AR <- unlist(lapply(Z0s.GEE.AR, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.ori.GEE.EX <- unlist(lapply(Z0s.GEE.EX, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.ori.GEE.IN <- unlist(lapply(Z0s.GEE.IN, function(x) {xx <- 1-pchisq(x^2, df=1); min(length(xx)*xx/rank(xx))}))
        simes.pv0s.GEE.AR <- rep(NA, n.perm)
        simes.pv0s.GEE.EX <- rep(NA, n.perm)
        simes.pv0s.GEE.IN <- rep(NA, n.perm)
        
        for (j in 1:n.perm) {	
          simes.pv0s.GEE.AR[j] <- (length(which(simes.pv0s.ori.GEE.AR[-j] < simes.pv0s.ori.GEE.AR[j]))+0.01)/(n.perm+0.01)
          simes.pv0s.GEE.EX[j] <- (length(which(simes.pv0s.ori.GEE.EX[-j] < simes.pv0s.ori.GEE.EX[j]))+0.01)/(n.perm+0.01)
          simes.pv0s.GEE.IN[j] <- (length(which(simes.pv0s.ori.GEE.IN[-j] < simes.pv0s.ori.GEE.IN[j]))+0.01)/(n.perm+0.01)
        }

        simes.pv.GEE.AR <- (length(which(simes.pv0s.GEE.AR < simes.pv.GEE.AR))+0.01)/(n.perm+0.01)
        simes.pv.GEE.EX <- (length(which(simes.pv0s.GEE.EX < simes.pv.GEE.EX))+0.01)/(n.perm+0.01)
        simes.pv.GEE.IN <- (length(which(simes.pv0s.GEE.IN < simes.pv.GEE.IN))+0.01)/(n.perm+0.01)
        
        l.Gamma <- length(Gamma)
        Tu.GEE.AR <- min(pvs.GEE.AR, simes.pv.GEE.AR)
        Tu.GEE.EX <- min(pvs.GEE.EX, simes.pv.GEE.EX)
        Tu.GEE.IN <- min(pvs.GEE.IN, simes.pv.GEE.IN)
        T0u.GEE.AR <- rep(NA, n.perm)
        T0u.GEE.EX <- rep(NA, n.perm)
        T0u.GEE.IN <- rep(NA, n.perm)
        for (j in 1:n.perm) {
          T0u.o.GEE.AR <- lapply(ahc.p.GEE.AR, function(x) x[j])
          T0u.o.GEE.EX <- lapply(ahc.p.GEE.EX, function(x) x[j])
          T0u.o.GEE.IN <- lapply(ahc.p.GEE.IN, function(x) x[j])
          T0u.p.GEE.AR <- lapply(ahc.p.GEE.AR, function(x) x[-j])
          T0u.p.GEE.EX <- lapply(ahc.p.GEE.EX, function(x) x[-j])
          T0u.p.GEE.IN <- lapply(ahc.p.GEE.IN, function(x) x[-j])
          T0u.GEE.AR[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.AR, T0u.p.GEE.AR))
          T0u.GEE.EX[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.EX, T0u.p.GEE.EX))
          T0u.GEE.IN[j] <- min(mapply(function(x, y) (length(which(x < y))+0.01)/(n.perm+0.01), T0u.o.GEE.IN, T0u.p.GEE.IN))
        }
        T0u.minp.GEE.AR <- apply(cbind(T0u.GEE.AR, simes.pv0s.GEE.AR),1,min)
        T0u.minp.GEE.EX <- apply(cbind(T0u.GEE.EX, simes.pv0s.GEE.EX),1,min)
        T0u.minp.GEE.IN <- apply(cbind(T0u.GEE.IN, simes.pv0s.GEE.IN),1,min)
        pv.opt.u.GEE.AR <- (length(which(T0u.minp.GEE.AR < Tu.GEE.AR))+0.01)/(n.perm+0.01)
        pv.opt.u.GEE.EX <- (length(which(T0u.minp.GEE.EX < Tu.GEE.EX))+0.01)/(n.perm+0.01)
        pv.opt.u.GEE.IN <- (length(which(T0u.minp.GEE.IN < Tu.GEE.IN))+0.01)/(n.perm+0.01)
        
        T <- min(Tu.GEE.AR, Tu.GEE.EX, Tu.GEE.IN)
        T0 <- apply(cbind(T0u.minp.GEE.AR, T0u.minp.GEE.EX, T0u.minp.GEE.IN),1,min)
        pv.opt <- (length(which(T0 < T))+0.01)/(n.perm+0.01)
        
        names(simes.pv.GEE.AR) <- "Simes.GEE(AR)"
        names(simes.pv.GEE.EX) <- "Simes.GEE(EX)"
        names(simes.pv.GEE.IN) <- "Simes.GEE(IN)"
        item.pvs.GEE.AR <- as.matrix(pvs.GEE.AR)
        item.pvs.GEE.EX <- as.matrix(pvs.GEE.EX)
        item.pvs.GEE.IN <- as.matrix(pvs.GEE.IN)
        colnames(item.pvs.GEE.AR) <- "uGEEHC(h).AR"
        colnames(item.pvs.GEE.EX) <- "uGEEHC(h).EX"
        colnames(item.pvs.GEE.IN) <- "uGEEHC(h).IN"
        rownames(item.pvs.GEE.AR) <- as.character(Gamma)
        rownames(item.pvs.GEE.EX) <- as.character(Gamma)
        rownames(item.pvs.GEE.IN) <- as.character(Gamma)
        
        ada.pvs <- c(pv.opt.u.GEE.AR, pv.opt.u.GEE.EX, pv.opt.u.GEE.IN, pv.opt)
        names(ada.pvs) <- c("GEEMiHC(AR)", "GEEMiHC(EX)", "GEEMiHC(IN)", "aGEEMiHC")

        otu.ids <- colnames(otu.tab)
        otu.ids <- sub("New.ReferenceOTU", "N", otu.ids)
        otu.ids <- sub("New.CleanUp.ReferenceOTU", "NC", otu.ids)
        pvs.GEE.AR <- ind.pvs.GEE.AR
        pvs.GEE.EX <- ind.pvs.GEE.EX
        pvs.GEE.IN <- ind.pvs.GEE.IN
        i0.GEE.AR <- which(pvs.GEE.AR < 0.00000001)
        i0.GEE.EX <- which(pvs.GEE.EX < 0.00000001)
        i0.GEE.IN <- which(pvs.GEE.IN < 0.00000001)
        i1.GEE.AR <- which(pvs.GEE.AR > 0.99999999)
        i1.GEE.EX <- which(pvs.GEE.EX > 0.99999999)
        i1.GEE.IN <- which(pvs.GEE.IN > 0.99999999)
        pvs.GEE.AR[i0.GEE.AR] <- 0.00000001
        pvs.GEE.EX[i0.GEE.EX] <- 0.00000001
        pvs.GEE.IN[i0.GEE.IN] <- 0.00000001
        pvs.GEE.AR[i1.GEE.AR] <- 0.99999999
        pvs.GEE.EX[i1.GEE.EX] <- 0.99999999
        pvs.GEE.IN[i1.GEE.IN] <- 0.99999999
        Is.GEE.AR <- order(order(pvs.GEE.AR))
        Is.GEE.EX <- order(order(pvs.GEE.EX))
        Is.GEE.IN <- order(order(pvs.GEE.IN))
        exp.HC.GEE.AR <- (Is.GEE.AR/p)/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        exp.HC.GEE.EX <- (Is.GEE.EX/p)/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        exp.HC.GEE.IN <- (Is.GEE.IN/p)/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        obs.HC.GEE.AR <- pvs.GEE.AR/sqrt(pvs.GEE.AR*(1-pvs.GEE.AR)/p)
        obs.HC.GEE.EX <- pvs.GEE.EX/sqrt(pvs.GEE.EX*(1-pvs.GEE.EX)/p)
        obs.HC.GEE.IN <- pvs.GEE.IN/sqrt(pvs.GEE.IN*(1-pvs.GEE.IN)/p)
        deviation.uGEEHC.AR <- abs(exp.HC.GEE.AR - obs.HC.GEE.AR)
        deviation.uGEEHC.EX <- abs(exp.HC.GEE.EX - obs.HC.GEE.EX)
        deviation.uGEEHC.IN <- abs(exp.HC.GEE.IN - obs.HC.GEE.IN)
        rank.uGEEHC.AR <- otu.ids[order(deviation.uGEEHC.AR, decreasing = T)]
        rank.uGEEHC.EX <- otu.ids[order(deviation.uGEEHC.EX, decreasing = T)]
        rank.uGEEHC.IN <- otu.ids[order(deviation.uGEEHC.IN, decreasing = T)]
        return(list(graph.els = list(otu.ids = otu.ids, exp.GEEHC.AR = exp.HC.GEE.AR, exp.GEEHC.EX = exp.HC.GEE.EX, exp.GEEHC.IN = exp.HC.GEE.IN, obs.GEEHC.AR = obs.HC.GEE.AR, obs.GEEHC.EX = obs.HC.GEE.EX, obs.GEEHC.IN = obs.HC.GEE.IN), rank.order = list(uGEEHC.AR = rank.uGEEHC.AR, uGEEHC.EX = rank.uGEEHC.EX, uGEEHC.IN = rank.uGEEHC.IN), simes.pv.GEE.AR = simes.pv.GEE.AR, simes.pv.GEE.EX = simes.pv.GEE.EX, simes.pv.GEE.IN = simes.pv.GEE.IN, ind.pvs.GEEMiHC.AR = item.pvs.GEE.AR, ind.pvs.GEEMiHC.EX = item.pvs.GEE.EX, ind.pvs.GEEMiHC.IN = item.pvs.GEE.IN, ada.pvs = ada.pvs))
      }
    
  }
