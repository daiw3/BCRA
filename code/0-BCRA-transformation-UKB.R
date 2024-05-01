TARV_transform <- function(geno, annotation, Z, direction = c("both", "dominant", "recessive"))
{
  # convert genotype from SnpMatrix to matrix and correct 0-1-2 coding
  # snpsum <- col.summary(geno)
  # geno <- t(as(geno, "numeric"))
  # coding.test <- snpsum[, 8] - snpsum[, 6]
  # geno[which(coding.test>0), ] <- 2 - geno[which(coding.test>0), ]
  
  direction <- match.arg(direction)
  if(nrow(geno) != length(annotation)) stop("Error: nrow(geno) != length(annotation)\n")
  if(nrow(geno) != length(Z)) stop("Error: nrow(geno) != length(Z)\n")
  gid <- unique(annotation)
  if(direction %in% c("dominant", "both")) {
    X.pos <- matrix(NA, ncol(geno), length(gid))
    colnames(X.pos) <- paste(gid, "+", sep = "")
  }
  if(direction %in% c("recessive", "both")) {
    X.neg <- matrix(NA, ncol(geno), length(gid))
    colnames(X.neg) <- paste(gid, "-", sep = "")
  }
  for(i in 1:length(gid)) {
    g <- gid[i]
    gindex <- annotation == g
    z <- Z[gindex]
    ord <- order(z, decreasing = T)
    idx <- which(gindex)[ord]
    if(direction %in% c("dominant", "both")) {
      X.pos[, i] <- apply(rbind(geno[idx, ] > 0, TRUE), 2, function(x) which(x)[1])
    }
    if(direction %in% c("recessive", "both")) {
      X.neg[, i] <- apply(rbind(geno[idx, ] > 1, TRUE), 2, function(x) which(x)[1])
    }
  }
  if(direction == "both") {
    X <- cbind(X.pos, X.neg)
  } else if(direction == "dominant") {
    X <- X.pos
  } else if(direction == "recessive") {
    X <- X.neg
  }
  return(X)
}
