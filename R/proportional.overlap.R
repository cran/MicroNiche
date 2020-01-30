proportional.overlap <-
function(df, sampleInfo, envInfo, q = 1.65){
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
  env.prop <- envInfo/sum(envInfo)
  Z <- data.frame("PO1,2" = rep(NA, nrow(Counts.ra)^2), "PO2,1" = rep(NA, nrow(Counts.ra)^2), "value" = rep(NA, nrow(Counts.ra)^2))
  tmp.res <- vector()
  tmp.nom1 <- vector()
  tmp.nom2 <- vector()
  for(o in 1:1){
    for(i in 1:nrow(Counts.ra)){
      for(j in 1:nrow(Counts.ra)){
        sp1 <- Counts.ra[i,]
        sp2 <- Counts.ra[j,]
        noms <- rownames(Counts.ra)
        f1 <- function(x,y) x-y
        subdf <- cbind(sampleInfo, as.numeric(sp1))
        colnames(subdf) <- c("SI", "RA")
        agg <- aggregate(as.numeric(subdf[,2]), by=list(Category=subdf[,1]), FUN=sum)
        envsubdf <- cbind(sampleInfo, as.numeric(env.prop))
        colnames(envsubdf) <- c("SI", "EI")
        envagg <- aggregate(as.numeric(envsubdf[,2]), by=list(Category=envsubdf[,1]), FUN=sum)
        tmp.mat.sp1 <- as.matrix(cbind(as.matrix(agg$x), as.matrix(envagg$x)))
        tmp.res.sp1 <- vector()
        for(q in 1:nrow(tmp.mat.sp1)){
          NB2 <- mapply(f1, max(tmp.mat.sp1[q,]), min(tmp.mat.sp1[q,]))
          tmp.res.sp1 <- c(tmp.res.sp1,NB2)
        }
        res.sp1 <- 1 - 0.5*sum(tmp.res.sp1)
        subdf <- cbind(sampleInfo, as.numeric(sp2))
        colnames(subdf) <- c("SI", "RA")
        agg <- aggregate(as.numeric(subdf[,2]), by=list(Category=subdf[,1]), FUN=sum)
        envsubdf <- cbind(sampleInfo, as.numeric(env.prop))
        colnames(envsubdf) <- c("SI", "EI")
        envagg <- aggregate(as.numeric(envsubdf[,2]), by=list(Category=envsubdf[,1]), FUN=sum)
        tmp.mat.sp2 <- as.matrix(cbind(as.matrix(agg$x), as.matrix(envagg$x)))
        tmp.res.sp2 <- vector()
        for(q in 1:nrow(tmp.mat.sp2)){
          NB3 <- mapply(f1, max(tmp.mat.sp2[q,]), min(tmp.mat.sp2[q,]))
          tmp.res.sp2 <- c(tmp.res.sp2,NB3)
        }
        res.sp2 <- 1 - 0.5*sum(tmp.res.sp2)
        tmp.vec <- c(res.sp1, res.sp2)
        PO.res <- 1-((max(tmp.vec)-min(tmp.vec))/sum(tmp.vec))
        tmp.nom1 <- rbind(tmp.nom1, noms[i])
        tmp.nom2 <- rbind(tmp.nom2, noms[j])
        tmp.res <- rbind(tmp.res, PO.res)
      }
    }
    Z$LO1.2 <- tmp.nom1
    Z$LO2.1 <- tmp.nom2
    Z$value <- tmp.res
    cast.df <- dcast(Z, LO1.2 ~ LO2.1, value.var = "value")
    return(cast.df)
  }
}
