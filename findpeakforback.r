library("quantmod")
whole <- read.delim("D:/matetelo50match/telomate-split/20.cov", header=FALSE, stringsAsFactors=FALSE)
tres <- 150
for (j in 1:20) {
  peakfor <- findPeaks(whole[1:length(whole[, j +2]), j + 2],0)
  peakbak <- findPeaks(whole[length(whole[, j + 2]):1, j + 2],0)
  peakinfofor <- data.frame(character(), integer(), integer(), stringsAsFactors=FALSE)
  peakinfobak <- data.frame(character(), integer(), integer(), stringsAsFactors=FALSE)
  for (i in 1:length(peakfor)) peakinfofor[i,] <- list(whole$V1[peakfor[i] - 1], whole$V2[peakfor[i] - 1], whole[peakfor[i] - 1, j + 2])
  for (i in 1:length(peakbak)) peakinfobak[i,] <- list(whole$V1[length(whole[, j + 2]) - peakbak[i] + 2], whole$V2[length(whole[, j + 2]) - peakbak[i] + 2], whole[length(whole[, j + 2]) - peakbak[i] + 2, j + 2])
  if ( ! (is.na(peakinfofor[,1]) && is.na(peakinfobak[,1]))) {
    peakgroupfor <- data.frame(character(), integer(), integer(), integer(), stringsAsFactors=FALSE)
    peakgroupbak <- data.frame(character(), integer(), integer(), integer(), stringsAsFactors=FALSE)
    k <- 1
    peakgroupfor[k,] <- list(peakinfofor[1,1], peakinfofor[1,2], peakinfofor[1,3], 1)
    for (i in 1:length(peakinfofor[,1])) {
      if (peakinfofor[i,2] != peakgroupfor[k,2]) {
        if (abs(peakinfofor[i,2] - peakgroupfor[k,2]) < tres){
          peakgroupfor[k,4] + 1 -> peakgroupfor[k,4]
          if (peakinfofor[i,3] >= peakgroupfor[k,3]){
            peakgroupfor[k,1] <- peakinfofor[i,1]
            peakgroupfor[k,2] <- peakinfofor[i,2]
            peakgroupfor[k,3] <- peakinfofor[i,3]
          }
        } else {
          k+1 -> k
          peakgroupfor[k,] <- list(peakinfofor[i,1], peakinfofor[i,2], peakinfofor[i,3], 1)
        }
      }
    }
    k <- 1
    peakgroupbak[k,] <- list(peakinfobak[1,1], peakinfobak[1,2], peakinfobak[1,3], 1)
    for (i in 1:length(peakinfobak[,1])) {
      if (peakinfobak[i,2] != peakgroupbak[k,2]) {
        if (abs(peakinfobak[i,2] - peakgroupbak[k,2]) < tres){
          peakgroupbak[k,4] + 1 -> peakgroupbak[k,4]
          if (peakinfobak[i,3] >= peakgroupbak[k,3]){
            peakgroupbak[k,1] <- peakinfobak[i,1]
            peakgroupbak[k,2] <- peakinfobak[i,2]
            peakgroupbak[k,3] <- peakinfobak[i,3]
          }
        } else {
          k+1 -> k
          peakgroupbak[k,] <- list(peakinfobak[i,1], peakinfobak[i,2], peakinfobak[i,3], 1)
        }
      }
    }
    peakgroup <- data.frame(character(), integer(), integer(), integer(), character(), stringsAsFactors=FALSE)
    k <- 1
    for (i in 1:length(peakgroupfor[,1])){
      for (m in 1:length(peakgroupbak[,1])){
        if (abs(peakgroupfor[i,2] - peakgroupbak[length(peakgroupbak[,1]) + 1 - m,2]) < tres){
          if (peakgroupfor[i,4] < peakgroupbak[length(peakgroupbak[,1]) + 1 - m,4]){ ##[i,4] = peak count
            peakgroup[k,] <- list(peakgroupfor[i,1], peakgroupfor[i,2], peakgroupfor[i,3], peakgroupfor[i,4], "for")
            k+1 -> k
          } else if (peakgroupfor[i,4] > peakgroupbak[length(peakgroupbak[,1]) + 1 - m,4]) { ##[i,4] = peak count
            peakgroup[k,] <- list(peakgroupbak[length(peakgroupbak[,1]) + 1 - m,1], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,2], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,3], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,4], "bak")
            k+1 -> k
          } else if (peakgroupfor[i,3] > peakgroupbak[length(peakgroupbak[,1]) + 1 - m,3]) { ##[i,3] = cov
            peakgroup[k,] <- list(peakgroupfor[i,1], peakgroupfor[i,2], peakgroupfor[i,3], peakgroupfor[i,4], "for")
            k+1 -> k
          } else if (peakgroupfor[i,3] < peakgroupbak[length(peakgroupbak[,1]) + 1 - m,3]) { ##[i,3] = cov
            peakgroup[k,] <- list(peakgroupbak[length(peakgroupbak[,1]) + 1 - m,1], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,2], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,3], peakgroupbak[length(peakgroupbak[,1]) + 1 - m,4], "bak")
            k+1 -> k
          } else {
            peakgroup[k,] <- list(peakgroupfor[i,1], peakgroupfor[i,2], peakgroupfor[i,3], peakgroupfor[i,4], "for(=bak)")
            k+1 -> k
          }
        }
      }
    }
    peaksingle <- data.frame(character(), integer(), integer(), integer(), character(), stringsAsFactors=FALSE)
    k <- 1    
    if  ( ! is.na(peakgroup[1,1])) {
      for (i in 1:length(peakgroupfor[,1])){
        memo <- 0
        for (m in 1:length(peakgroup[,1])){
         if (abs(peakgroupfor[i,2] - peakgroup[m,2]) < tres) memo+1 -> memo
        }
        if (memo == 0){
         peaksingle[k,] <- list(peakgroupfor[i,1], peakgroupfor[i,2], peakgroupfor[i,3], peakgroupfor[i,4], "for(single)")
         k+1 -> k
        } 
      }
      for (i in 1:length(peakgroupbak[,1])){
        memo <- 0
        for (m in 1:length(peakgroup[,1])){
         if (abs(peakgroupbak[i,2] - peakgroup[m,2]) < tres) memo+1 -> memo
        }
        if (memo == 0){
          peaksingle[k,] <- list(peakgroupbak[i,1], peakgroupbak[i,2], peakgroupbak[i,3], peakgroupbak[i,4], "bak(single)")
          k+1 -> k
        }
      }
      prev <- length(peakgroup[,1])
    } else {
      for (i in 1:length(peakgroupfor[,1])){
        peaksingle[k,] <- list(peakgroupfor[i,1], peakgroupfor[i,2], peakgroupfor[i,3], peakgroupfor[i,4], "for(single)")
        k+1 -> k
      }
      for (i in 1:length(peakgroupbak[,1])){
        peaksingle[k,] <- list(peakgroupbak[i,1], peakgroupbak[i,2], peakgroupbak[i,3], peakgroupbak[i,4], "bak(single)")
        k+1 -> k
      }
      prev <- length(peakgroup[,1])
    }
    if (length(peaksingle[,1]) > 0){
      for (i in 1:length(peaksingle[,1])){
        peakgroup[prev+i,] <- list(peaksingle[i,1], peaksingle[i,2], peaksingle[i,3], peaksingle[i,4], peaksingle[i,5])
      }
    }
    write.table(peakgroup,file = paste("D:/matetelo50match/telomate-split/" ,j , "-peak.cov", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE)
  } else {
    write.table(peakinfofor,file = paste("D:/matetelo50match/telomate-split/" ,j , "-peak.cov", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
