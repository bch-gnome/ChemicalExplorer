## CORRELATION CIRCOS #########################################################
load("data/GC-MS/group_correlation.RData")
gc_crr_grop <- list(
  "All Samples" = crr_all$D,
  "Case Samples" = crr_cas$D,
  "Control Samples" = crr_ctr$D,
  "Mother Samples" = crr_mth$D,
  "Father Samples" = crr_fth$D
)
rm(crr_all, crr_cas, crr_ctr, crr_mth, crr_fth)
completeMatrix <- function(mtr, nms) {
  dsp <- matrix(0L, ncol=length(nms), nrow=length(nms))
  colnames(dsp) <- rownames(dsp) <- nms
  for(ii in rownames(mtr)) {
    for(jj in colnames(mtr)) {
      dsp[ii, jj] <- mtr[ii, jj]
    }
  }
  dsp
}
getCorGR <- function(set) {
  mat <- melt(gc_crr_grop[[set]])
  mat <- mat[mat$Var1 != mat$Var2, ]
  mat <- mat[order(mat$value, decreasing = TRUE), ]
  mat_f <- rbind(head(mat, n=36), tail(mat, n=36))
  mat_f <- mat_f[mat_f$value != 0, ]
  mat_f <- dcast(mat_f, Var1~Var2, value.var="value", fill=0)
  rownames(mat_f) <- mat_f[, 1]
  mat_f <- mat_f[ , -1]
  mat_f <- completeMatrix(mat_f, gc_data$fd$uid)
  
  network <- graph_from_adjacency_matrix(as.matrix(mat_f),
    weighted = TRUE, mode = "undirected", diag = FALSE)
  V( network )$color <- gc_data$colour[
    gc_data$fd[names( V( network ) ), "Proposed.Group"]
    ]
  V( network )$name <- gc_data$fd[ names( V( network ) ), "Chemical.Name" ]
  network
}