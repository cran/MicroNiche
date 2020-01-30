levins.overlap <-
function(df, q = 1.65){
  #Perform lognormal species rank model and plot LOQ
  Taxa <- df[,1]
  Counts <- df[,-1]
  Counts <- cbind.data.frame(Taxa, Counts)
  melt.Counts <- melt(Counts, id.vars = "Taxa")
  cast.Counts <- dcast(melt.Counts, Taxa ~ Taxa, fun.aggregate = sum)
  tmpdf <- cast.Counts[,-1]
  tmpdf <- cbind.data.frame(colnames(tmpdf), as.numeric(colSums(tmpdf)))
  colnames(tmpdf) <- c("Taxa", "Sum")
  tmp5 <- tmpdf[order(tmpdf$Sum, decreasing = T),]
  Ranks <- c(1:nrow(tmpdf))
  tmp6 <- cbind.data.frame(tmp5, Ranks)
  plot(log(Sum) ~ Ranks, data = tmp6, pch = 16, cex = 1.5, type = "b", ylim = c(0, max(log(tmpdf$Sum))+1), ylab = "Log Abundance", xlab = "Taxon Rank", cex.lab = 1.6, las = 1)
  legend("topright", pch = c(16, 1, 1), col = c("black", "red", "blue"), legend = c("Data", "Model", "LOQ"), bty = "n", cex = 1.5)
  SR <- function(S, a, R){
    S*exp(-a^2*R^2)
  }
  a = sqrt((log(max(tmpdf$Sum))/min(tmpdf$Sum))/max(Ranks)^2) 
  res <- SR(log(max(tmpdf$Sum)), a, Ranks)
  points(res, col = "red", type = "b", cex = 1.5)
  maxLOQ <- rep(q*sd(res), length(Ranks))
  LOQlim <- q*sd(res)
  points(maxLOQ, col = "blue", type = "b")
  rownames(Counts) <- Taxa
  Counts <- Counts[,-1]
  cleanCounts <- c(1:ncol(Counts))
  checkTaxa <- vector()
  #Remove taxa below lOQ
  for(y in 1:nrow(Counts)){
    query <- log(rowSums(Counts[y,]))
    curr.row <- Counts[y,]
    curr.nom <- as.character(Taxa[y])
    if(query < LOQlim){
      LOQnom <- paste(curr.nom, "*", sep = "")
      checkTaxa <- rbind(checkTaxa, LOQnom)
    }else{
      checkTaxa <- rbind(checkTaxa, curr.nom)
    }
  }
  #Transform Counts data frame to relative proportions
  rownames(Counts) <- checkTaxa
  cleanCounts <- Counts[-1,]
  Counts.ra <- (sweep(cleanCounts, 1,rowSums(cleanCounts), '/'))
  Z <- data.frame("LO1,2" = rep(NA, nrow(Counts.ra)^2), "LO2,1" = rep(NA, nrow(Counts.ra)^2), "value" = rep(NA, nrow(Counts.ra)^2))
  tmp.res <- vector()
  tmp.nom1 <- vector()
  tmp.nom2 <- vector()
  for(o in 1:1){
    for(i in 1:nrow(Counts.ra)){
      for(j in 1:nrow(Counts.ra)){
        sp1 <- Counts.ra[i,]
        sp2 <- Counts.ra[j,]
        noms <- rownames(Counts.ra)
        f1 <- function(x,y) x*y
        W <- mapply(f1, sp1, sp2)
        K <- apply(as.matrix(sp1), 1, function(j) j^2)
        WK <- cbind(sum(W),sum(K))
        res1 <- min(WK)/max(WK)
        tmp.nom1 <- rbind(tmp.nom1, noms[i])
        tmp.nom2 <- rbind(tmp.nom2, noms[j])
        tmp.res <- rbind(tmp.res, res1)
      }
    }
    Z$LO1.2 <- tmp.nom1
    Z$LO2.1 <- tmp.nom2
    Z$value <- tmp.res
    cast.df <- dcast(Z, LO1.2 ~ LO2.1, value.var = "value")
    return(cast.df)
  }
}
