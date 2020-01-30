hurlberts.Bn <-
function(df, R, sampleInfo, envInfo, q = 1.65){
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
  #Transform Counts data frame and envInfo to relative proportions
  rownames(Counts) <- Taxa
  Counts <- Counts[,-1]
  Counts.ra <- (sweep(Counts, 1,rowSums(Counts), '/'))
  env.prop <- envInfo/sum(envInfo)
  #Perform null model testing for Bn
  Z <- data.frame("Null.Bn" = rep(NA, 999))
  tmp2 <- data.frame("Bn" = rep(0,0), "P.val" = rep(0,0), "Adj.P" = rep(0,0))
  for(k in 1:999){
    set.seed(k)
    x<-runif(R)
    RA.proportions <- x/sum(x)
    NB1 <- apply(as.matrix(RA.proportions), 1, function(y) y^2)
    envsubdf <- cbind(sampleInfo, as.numeric(env.prop))
    colnames(envsubdf) <- c("SI", "EI")
    envagg <- aggregate(as.numeric(envsubdf[,2]), by=list(Category=envsubdf[,1]), FUN=sum)
    tmp.mat <- as.matrix(cbind(NB1, as.matrix(envagg$x)))
    foo <- function(z,e) z/e
    NB2 <- mapply(foo, tmp.mat[,1],tmp.mat[,2])
    Z$Null.Bn[k] <- 1/sum(NB2)
  }
  qts <- quantile(Z$Null.Bn, probs = c(0.05, 0.95))
  myplot <- ggplot(Z, aes(x = Z$Null.Bn)) +
    geom_histogram(alpha=0.6, position="identity", bins = 250) +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=22)) +
    labs(x = expression(B[n]), y = "Frequency") +
    geom_vline(xintercept= c(as.numeric(qts[1]), as.numeric(qts[2])), linetype="dashed", color = "red", size=1) 
  print(myplot)
  #Calculate Bn for each taxa
  for(i in 1:nrow(Counts.ra)){
    query <- Counts.ra[i,]
    subdf <- cbind(sampleInfo, as.numeric(query))
    colnames(subdf) <- c("SI", "RA")
    agg <- aggregate(as.numeric(subdf[,2]), by=list(Category=subdf[,1]), FUN=sum)
    W <- apply(as.matrix(agg$x), 2, function(j) j^2)
    envsubdf <- cbind(sampleInfo, as.numeric(env.prop))
    colnames(envsubdf) <- c("SI", "EI")
    envagg <- aggregate(as.numeric(envsubdf[,2]), by=list(Category=envsubdf[,1]), FUN=sum)
    tmp.mat <- as.matrix(cbind(W, as.matrix(envagg$x)))
    foo <- function(z,e) z/e
    res1 <- mapply(foo, tmp.mat[,1],tmp.mat[,2])
    res2 <- 1/sum(res1)
    X <- mean(Z$Null.Bn)
    z <- (res2 - X)/(3/(sqrt(999)))
    pval <- 2*pnorm(-abs(z))
    tmp <- data.frame(Bn = res2, P.val = pval)
    rownames(tmp) <- rownames(Counts[i,])
    tmp2 <- rbind(tmp2, tmp)
  }
  P.adj <- p.adjust(tmp2$P.val, method = "BH")
  tmp2 <- cbind(tmp2, P.adj)
  #Flag the taxa that are below the LOQ
  Below.LOQ <- data.frame("Below.LOQ" = rep("N", nrow(tmp2)), stringsAsFactors = F)
  for(y in 1:nrow(tmpdf)){
    for(w in 1:nrow(tmp2)){
      Nom1 <- as.character(tmpdf$Taxa)[y]
      Nom2 <- as.character(rownames(tmp2))[w]
      checking <-  match(Nom1, Nom2, nomatch = 0)
      if(checking > 0){
        query <- log(tmpdf$Sum[y])
        if(query < LOQlim){
          Below.LOQ$Below.LOQ[w] <- "Y"
        }
      }
    }
  }
  tmp2 <- cbind(tmp2, Below.LOQ)
  return(tmp2)
}
