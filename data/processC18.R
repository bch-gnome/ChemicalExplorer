## DEAFULT LOADING ############################################################
load("data/LC-MS/c18_data.RData")

c18_data <- list(
  colour = c("HMDB" = "red", "KEGG" = "green", "LIPID" = "blue"),
  measures = c18_dta,
  table = c18_fd,
  pd = gc_data$pd,
  fd = c18_fd
)
rownames(c18_data$measures) <- c18_data$measures$uid
c18_data$measures <- t(c18_data$measures[, -1])
rm(c18_fd, c18_dta)

## PIE CHART ##################################################################
pieC18 <- function() {
  dta <- as.data.frame(table(c18_data$fd$Source))
  plot_ly(dta, labels = ~Var1, values = ~Freq, type = 'pie')
}

## NAMES ARRAY ################################################################
namesC18 <- function(group) {
  sel <- as.character(c18_data$fd$Name[c18_data$fd$Source == group])
  sel
}

## MZ-RT ARRAY ################################################################
mztimeC18 <- function(metab) {
  dta <- c18_data$fd[c18_data$fd$Name == metab, c("mz", "time")]
  sel <- apply(dta, 1, function(row) {
    paste0("mz: ", row[1], "; time: ", row[2])
  })
  as.character(sel)
}

## BOXPLOT ####################################################################
boxplotC18 <- function(metab, mzt) {
  mzt <- strsplit(gsub("mz: ", "", gsub("; time: ", "-", mzt)), "-")[[1]]
  dta <- c18_data$fd[c18_data$fd$Name == metab, c("mz", "time", "mid")]
  dta <- dta[dta$mz == mzt[1] & dta$time == mzt[2], "mid"]
  dta <- as.data.frame(log2(c18_data$measures[ , dta, drop=FALSE]))
  dta$Group <- as.character(c18_data$pd[rownames(dta), c("Group")])
  colnames(dta) <- c("ppm", "group")
  
  plot_ly(dta, y = ~ppm, color = ~group, type = "box") %>%
    layout(
      title = paste(metab, "\n(mz: ", mzt[1], "; time: ", mzt[2], ")"),
      xaxis = list(title = ""), yaxis = list(title = "log2(ppm)")
    )
}