
day_cols = c("#fcfdbf", "#feb078", "#f1605d", "#b73779", "#721f81", "#2c115f", "#000004")
names(day_cols) = c(0,3,4,5,6,7,8)

day_cols2 = c("#febb81", "#f8765c", "#d3436e", "#982d80", "#5f187f", "#221150", "#000004")
names(day_cols2) = c(0,3,4,5,6,7,8)
## calculate MVS score and color by hong's color scheme

## ADDITIONAL FUNCTIONS FOR ANALYSIS
geomMean <- function(x, na.rm = T){
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm)
    x <- x[!is.na(x)]
  if (any(x < 0))
    stop("'x' contains negative value(s)")
  ### direct method causes overflow errors, use log method instead
  ### return(prod(x)^(1/length(x)))
  return(exp(sum(log(x[x>0]))/length(x)))
}
getGeneScores <- function(geneMtx, pos, neg, makePos = TRUE, out.missing=TRUE){
  if(out.missing){
    missingpos = pos[!(pos %in% rownames(geneMtx))]
    missingneg = neg[!(neg %in% rownames(geneMtx))]
    missing = c(missingpos,missingneg)
    if(length(missing)>0){
      cat("Missing these genes:",missing,"\n")
    }
  }
  pos = pos[pos %in% rownames(geneMtx)]
  neg = neg[neg %in% rownames(geneMtx)]
  if(makePos){ #If I need to make everything positive
    if(any(geneMtx < 0, na.rm=T)){
      geneMtx <- geneMtx - min(geneMtx, na.rm=T)
    }
  }
  if(length(neg)>=1 && length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T) - apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(pos)>=1){
    scores <- apply(geneMtx[pos, , drop=F], 2, geomMean, na.rm=T)
  }else if(length(neg)>=1){
    scores <- -1*apply(geneMtx[neg, , drop=F], 2, geomMean, na.rm=T)
  }
  if(any(is.nan(scores))){ #this means all upgenes or all downgenes were missing from some samples
    nanIndex = which(is.nan(scores))
    for(i in nanIndex){
      my.score = c(geomMean(geneMtx[pos,i,drop=F],na.rm=T),geomMean(geneMtx[neg,i,drop=F],na.rm=T))
      my.score = my.score[!is.nan(my.score)]
      if(length(my.score) == 1){
        scores[i] = my.score
      }else if(length(my.score) == 2){
        scores[i] = my.score[1] - my.score[2]
      }
    }
  }
  return(scores)
}
