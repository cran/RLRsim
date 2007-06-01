`extract.lmeDesign` <-
function(m)
{
    require(mgcv)
    start.level = 1
    data<-m$data
    grps <- getGroups(m)
    n <- length(grps)
    if (is.null(m$modelStruct$varStruct))
        w <- rep(m$sigma, n)
    else {
        w <- 1/varWeights(m$modelStruct$varStruct)
        group.name <- names(m$groups)
        order.txt <- paste("ind<-order(data[[\"", group.name[1],
            "\"]]", sep = "")
        if (length(m$groups) > 1)
            for (i in 2:length(m$groups)) order.txt <- paste(order.txt,
                ",data[[\"", group.name[i], "\"]]", sep = "")
        order.txt <- paste(order.txt, ")")
        eval(parse(text = order.txt))
        w[ind] <- w
        w <- w * m$sigma
    }
    if (is.null(m$modelStruct$corStruct))
        V <- diag(n)
    else {
        c.m <- corMatrix(m$modelStruct$corStruct)
        if (!is.list(c.m))
            V <- c.m
        else {
            V <- matrix(0, n, n)
            gr.name <- names(c.m)
            n.g <- length(c.m)
            j0 <- 1
            ind <- ii <- 1:n
            for (i in 1:n.g) {
                j1 <- j0 + nrow(c.m[[i]]) - 1
                V[j0:j1, j0:j1] <- c.m[[i]]
                ind[j0:j1] <- ii[grps == gr.name[i]]
                j0 <- j1 + 1
            }
            V[ind, ] <- V
            V[, ind] <- V
        }
    }
    V <- as.vector(w) * t(as.vector(w) * V)
    X <- list()
    grp.dims <- m$dims$ncol
    Zt <- model.matrix(m$modelStruct$reStruct, data)
    cov <- as.matrix(m$modelStruct$reStruct)
    i.col <- 1
    n.levels <- length(m$groups)
    Z <- matrix(0, n, 0)
    if (start.level <= n.levels) {
        for (i in 1:(n.levels - start.level + 1)) {
            if(length(levels(m$groups[[n.levels-i+1]]))!=1)
            {
            X[[1]] <- model.matrix(~m$groups[[n.levels - i +
                1]] - 1, contrasts.arg = c("contr.treatment",
                "contr.treatment"))
            }
            else X[[1]]<-matrix(1)
            X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] -
                1)])
            i.col <- i.col + grp.dims[i]
            Z <- cbind(tensor.prod.model.matrix(X),Z)
        }
        Vr <- matrix(0, ncol(Z), ncol(Z))
        start <- 1
        for (i in 1:(n.levels - start.level + 1)) {
            k <- n.levels - i + 1
            for (j in 1:m$dims$ngrps[i]) {
                stop <- start + ncol(cov[[k]]) - 1
                 Vr[ncol(Z)+1-(stop:start),ncol(Z)+1-(stop:start)] <- cov[[k]]
                 start <- stop + 1
            }
        }

        Vlambda <- V/m$sigma^2 + Z %*% Vr %*% t(Z)
    }
X<-model.matrix(formula(m$call$fixed),data)
y<-as.vector(matrix(m$residuals,nc=NCOL(m$residuals))[,NCOL(m$residuals)] +matrix(m$fitted,nc=NCOL(m$fitted))[,NCOL(m$fitted)])
return(list(
             Vlambda=Vlambda, #Cov(y)/Var(Error)
             V=V, #Cov(Error)
             Vr=Vr, #Cov(RanEf)/Var(Error)
             X=X,
             Z=Z,
             sigmasq=m$sigma^2,
             lambda=unique(diag(Vr)),
             y=y
           )
      )
}

