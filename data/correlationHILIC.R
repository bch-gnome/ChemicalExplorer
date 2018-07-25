## DEAFULT LOADING ############################################################
load("data/LC-MS/hilic_corr_1p.RData")
load("data/LC-MS/hilic_annot2.RData")
hilic_hmdb_sel <- hilic_hmdb_f$Name
names(hilic_hmdb_sel) <- hilic_hmdb_f$mid
hilic_kegg_sel <- hilic_kegg_f$Name
names(hilic_kegg_sel) <- hilic_kegg_f$mid
hilic_lipd_sel <- hilic_lipd_f$Name
names(hilic_lipd_sel) <- hilic_lipd_f$mid
rm(hilic_hmdb_f, hilic_lipd_f)
hilic_kegg_f <- hilic_kegg_f[hilic_kegg_f$MatchCategory.1 == "Unique", ]
hilic_kegg_f <- hilic_kegg_f[!duplicated(hilic_kegg_f$KEGGID), ]
phe <- phe[!duplicated(phe$BioSTOR), ]

## DEAFULT LOADING ############################################################
load("data/LC-MS/kegg_annotation.RData")
hilic_kegg <- kilic_kegg; rm(kilic_kegg)
for(elm in names(c18_kegg)) {
  names(c18_kegg[[elm]]$compounds_id) <- c18_kegg[[elm]]$compounds_name
}
for(elm in names(hilic_kegg)) {
  names(hilic_kegg[[elm]]$compounds_id) <- hilic_kegg[[elm]]$compounds_name
}
rm(elm)

## SHORTING NAMES #############################################################
th_char <- 15

len_hmdb <- nchar(hilic_hmdb_sel)
hilic_hmdb_sel[len_hmdb > th_char] <-
  paste(substr(hilic_hmdb_sel[len_hmdb > th_char], 1, th_char), "...")

len_kegg <- nchar(hilic_kegg_sel)
hilic_kegg_sel[len_kegg > th_char] <-
  paste(substr(hilic_kegg_sel[len_kegg > th_char], 1, th_char), "...")

len_lipd <- nchar(hilic_lipd_sel)
hilic_lipd_sel[len_lipd > th_char] <-
  paste(substr(hilic_lipd_sel[len_lipd > th_char], 1, th_char), "...")

rm(len_lipd, len_kegg, len_hmdb)

## PLOT CORRELATION ###########################################################
correlationHILIC <- function(set1, set2) {
  mat <- valH_f
  
  set <- c(set1, set2)
  if("HMDB" %in% set & "KEGG" %in% set) {
    mat <- mat[mat$Var1 %in% names(hilic_hmdb_sel) &
                 mat$Var2 %in% names(hilic_kegg_sel), ]
    mat$Chemical1 <- hilic_hmdb_sel[mat$Var1]
    mat$Chemical2 <- hilic_kegg_sel[mat$Var2]
    K1 <- "HMDB"
    K2 <- "KEGG"
  } else if("HMDB" %in% set & "LIPID" %in% set) {
    mat <- mat[mat$Var1 %in% names(hilic_hmdb_sel) &
                 mat$Var2 %in% names(hilic_lipd_sel), ]
    mat$Chemical1 <- hilic_hmdb_sel[mat$Var1]
    mat$Chemical2 <- hilic_lipd_sel[mat$Var2]
    K1 <- "HMDB"
    K2 <- "LIPID"
  } else {
    mat <- mat[mat$Var1 %in% names(hilic_lipd_sel) &
                 mat$Var2 %in% names(hilic_kegg_sel), ]
    mat$Chemical1 <- hilic_lipd_sel[mat$Var1]
    mat$Chemical2 <- hilic_kegg_sel[mat$Var2]
    K1 <- "LIPID"
    K2 <- "KEGG"
  }
  
  p1 <- ggplot(mat, aes(Chemical1, Chemical2, fill=value)) + 
    geom_tile() +
    theme_bw() +
    ggtitle(paste0("Correlation between ", K1, " and ", K2)) +
    xlab(K1) + ylab(K2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient2(name="Correlation", low="red", mid="white", high="blue")
  
  ggplotly(p1)
}

## PATHWAYS NAMES #############################################################
keggNameHILIC <- function() {
  sel <- as.character(sapply(hilic_kegg, '[[', 'name'))
  sel <- sel[
    sapply(hilic_kegg, function(it) length(it$compounds_name)) != 0
  ]
  sel[order(sel)]
}
lc_samples_grop <- list(
  "All Samples" = phe$BioSTOR,
  "Case Samples" = phe$BioSTOR[phe$Group == "Proband"],
  "Control Samples" = phe$BioSTOR[phe$Group == "Control"],
  "Mother Samples" = phe$BioSTOR[phe$Group == "Mother"],
  "Father Samples" = phe$BioSTOR[phe$Group == "Father"]
)

getKEGGcorHILIC <- function(path, set) {
  message("HILIC: ", set, " - ", path)
  sel1 <- which(sapply(hilic_kegg, '[[', 'name') == path)
  sel <- hilic_kegg[[sel1]]$compounds_name
  mat <- hilic_kegg_f[hilic_kegg_f$KEGGID %in% sel, ]
  rownames(mat) <- mat$KEGGID
  mat <- mat[ , intersect(colnames(mat), lc_samples_grop[[set]])]
  mat <- cor(t(mat))
  mat <- completeMatrix(mat, sel)
  
  mat_r <- melt(mat)
  colnames(mat_r) <- c("KEGG1", "KEGG2", "Correlation")
  mat_r$Chemical1 <- as.character(hilic_kegg[[sel1]]$compounds_id[mat_r$KEGG1])
  mat_r$Chemical2 <- as.character(hilic_kegg[[sel1]]$compounds_id[mat_r$KEGG2])
  
  plt <- ggplot(mat_r, 
                aes_string(x="Chemical1", y="Chemical2", fill="Correlation")) + 
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("") + ylab("") +
    scale_fill_gradient2(name="Correlation", low="red", mid="white", high="blue") +
    ggtitle(c18_kegg[[sel1]]$name)
  ggplotly(plt)
}