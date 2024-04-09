calculate_statistics <- function(tbl) {
  P <- sum(tbl[, 2])
  N <- sum(tbl[, 1])
  TP <- tbl[2, 2]
  TN <- tbl[1, 1]
  FP <- tbl[2, 1]
  # FN <- tbl[2, 1]
  TPR <- TP / P
  FPR <- FP / N
  # FDR <- FP / (TP + FP)
  ACC <- (TP + TN) / (P + N)
  # data.frame(P = P, N = N, TP = TP, TN = TN, FP = FP, 
  #            FN = FN, TPR = TPR, FPR = FPR, FDR = FDR)
  data.frame(TPR = TPR, FPR = FPR, ACC = ACC)
}


