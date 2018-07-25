## DEAFULT LOADING ############################################################
load("data/LC-MS/hilic_data.RData")

hilic_data <- list(
  colour = c("HMDB" = "red", "KEGG" = "green", "LIPID" = "blue"),
  measures = hilic_dta,
  table = hilic_fd,
  pd = gc_data$pd,
  fd = hilic_fd
)
rownames(hilic_data$measures) <- hilic_data$measures$uid
hilic_data$measures <- t(hilic_data$measures[, -1])
rm(hilic_fd, hilic_dta)

## PIE CHART ##################################################################
pieHILIC <- function() {
  dta <- as.data.frame(table(hilic_data$fd$Source))
  plot_ly(dta, labels = ~Var1, values = ~Freq, type = 'pie')
}

## NAMES ARRAY ################################################################
namesHILIC <- function(group) {
  sel <- as.character(hilic_data$fd$Name[hilic_data$fd$Source == group])
  sel
}

## MZ-RT ARRAY ################################################################
mztimeHILIC <- function(metab) {
  dta <- hilic_data$fd[hilic_data$fd$Name == metab, c("mz", "time")]
  sel <- apply(dta, 1, function(row) {
    paste0("mz: ", row[1], "; time: ", row[2])
  })
  as.character(sel)
}

## BOXPLOT ####################################################################
boxplotHILIC <- function(metab, mzt) {
  mzt <- strsplit(gsub("mz: ", "", gsub("; time: ", "-", mzt)), "-")[[1]]
  dta <- hilic_data$fd[hilic_data$fd$Name == metab, c("mz", "time", "mid")]
  dta <- dta[dta$mz == mzt[1] & dta$time == mzt[2], "mid"]
  dta <- as.data.frame(log2(hilic_data$measures[ , dta, drop=FALSE]))
  dta$Group <- as.character(hilic_data$pd[rownames(dta), c("Group")])
  colnames(dta) <- c("ppm", "group")
  plot_ly(dta, y = ~ppm, color = ~group, type = "box") %>%
    layout(
      title = paste(metab, "\n(mz: ", mzt[1], "; time: ", mzt[2], ")"),
      xaxis = list(title = ""), yaxis = list(title = "log2(ppm)")
    )
}