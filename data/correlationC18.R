## DEAFULT LOADING ############################################################
load("data/LC-MS/c18_corr_1p.RData")
load("data/LC-MS/c18_annot2.RData")
c18_hmdb_sel <- c18_hmdb_f$Name
names(c18_hmdb_sel) <- c18_hmdb_f$mid
c18_kegg_sel <- c18_kegg_f$Name
names(c18_kegg_sel) <- c18_kegg_f$mid
c18_lipd_sel <- c18_lipd_f$Name
names(c18_lipd_sel) <- c18_lipd_f$mid
rm(c18_hmdb_f, c18_lipd_f)
c18_kegg_f <- c18_kegg_f[c18_kegg_f$MatchCategory.1 == "Unique", ]
c18_kegg_f <- c18_kegg_f[!duplicated(c18_kegg_f$KEGGID), ]
phe <- phe[!duplicated(phe$BioSTOR), ]

# ## DEAFULT LOADING ############################################################
load("data/LC-MS/kegg_annotation.RData")
hilic_kegg <- kilic_kegg; rm(kilic_kegg)
for(elm in names(c18_kegg)) {
  names(c18_kegg[[elm]]$compounds_id) <- c18_kegg[[elm]]$compounds_name
}
for(elm in names(hilic_kegg)) {
  names(hilic_kegg[[elm]]$compounds_id) <- hilic_kegg[[elm]]$compounds_name
}
rm(elm)
# load("data/LC-MS/kegg_annotation_c18.Rdata")

## SHORTING NAMES #############################################################
th_char <- 15

len_hmdb <- nchar(c18_hmdb_sel)
c18_hmdb_sel[len_hmdb > th_char] <-
  paste(substr(c18_hmdb_sel[len_hmdb > th_char], 1, th_char), "...")

len_kegg <- nchar(c18_kegg_sel)
c18_kegg_sel[len_kegg > th_char] <-
  paste(substr(c18_kegg_sel[len_kegg > th_char], 1, th_char), "...")

len_lipd <- nchar(c18_lipd_sel)
c18_lipd_sel[len_lipd > th_char] <-
  paste(substr(c18_lipd_sel[len_lipd > th_char], 1, th_char), "...")

rm(len_lipd, len_kegg, len_hmdb)

## PLOT CORRELATION ###########################################################
correlationC18 <- function(set1, set2) {
  mat <- valC18_f
  
  set <- c(set1, set2)
  if("HMDB" %in% set & "KEGG" %in% set) {
    mat <- mat[mat$Var1 %in% names(c18_hmdb_sel) &
               mat$Var2 %in% names(c18_kegg_sel), ]
    mat$Chemical1 <- c18_hmdb_sel[mat$Var1]
    mat$Chemical2 <- c18_kegg_sel[mat$Var2]
    K1 <- "HMDB"
    K2 <- "KEGG"
  } else if("HMDB" %in% set & "LIPID" %in% set) {
    mat <- mat[mat$Var1 %in% names(c18_hmdb_sel) &
               mat$Var2 %in% names(c18_lipd_sel), ]
    mat$Chemical1 <- c18_hmdb_sel[mat$Var1]
    mat$Chemical2 <- c18_lipd_sel[mat$Var2]
    K1 <- "HMDB"
    K2 <- "LIPID"
  } else {
    mat <- mat[mat$Var1 %in% names(c18_lipd_sel) &
               mat$Var2 %in% names(c18_kegg_sel), ]
    mat$Chemical1 <- c18_lipd_sel[mat$Var1]
    mat$Chemical2 <- c18_kegg_sel[mat$Var2]
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
keggNameC18 <- function() {
  sel <- as.character(sapply(c18_kegg, '[[', 'name'))
  sel <- sel[
    sapply(c18_kegg, function(it) length(it$compounds_name)) != 0
    ]
  sel[order(sel)]
}

# # path <- "alpha-Linolenic acid metabolism"
# path <- "Tryptophan metabolism"
# # set <- "Case Samples"
# set <- "All Samples"
# 
# xx <- c18_kegg_f[
#   c18_kegg_f$KEGGID %in% c("C05648", "C16518"),
# ]
# rownames(xx) <- xx$KEGGID
# xx <- xx[ , 21:ncol(c18_kegg_f)]
# message("All")
# xx2 <- xx[ , intersect(colnames(xx), lc_samples_grop[["All Samples"]])]
# rcorr(t(xx2), type="pearson")
# message("Cases")
# xx2 <- xx[ , intersect(colnames(xx), lc_samples_grop[["Case Samples"]])]
# rcorr(t(xx2), type="pearson")
# message("Controls")
# xx2 <- xx[ , intersect(colnames(xx), lc_samples_grop[["Control Samples"]])]
# rcorr(t(xx2), type="pearson")
# message("Mothers")
# xx2 <- xx[ , intersect(colnames(xx), lc_samples_grop[["Mother Samples"]])]
# rcorr(t(xx2), type="pearson")
# message("Fathers")
# xx2 <- xx[ , intersect(colnames(xx), lc_samples_grop[["Father Samples"]])]
# rcorr(t(xx2), type="pearson")
# 

getKEGGcorC18 <- function(path, set) {
  message("C18: ", set, " - ", path)
  sel1 <- which(sapply(c18_kegg, '[[', 'name') == path)
  sel <- c18_kegg[[sel1]]$compounds_name
  mat <- c18_kegg_f[c18_kegg_f$KEGGID %in% sel, ]
  rownames(mat) <- mat$KEGGID
  mat <- mat[ , intersect(colnames(mat), lc_samples_grop[[set]])]
  mat <- cor(t(mat))

  
  mat_r <- melt(mat)
  colnames(mat_r) <- c("KEGG1", "KEGG2", "Correlation")
  mat_r$Chemical1 <- as.character(c18_kegg[[sel1]]$compounds_id[mat_r$KEGG1])
  mat_r$Chemical2 <- as.character(c18_kegg[[sel1]]$compounds_id[mat_r$KEGG2])
  
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
