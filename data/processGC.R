## DEAFULT LOADING ############################################################
load("data/GC-MS/colourChemicals.RData")
load("data/GC-MS/Transformed_Filtered_Data.RData")
grd <- read.xlsx("data/GC-MS/CompoundsRationale.xlsx",
                 sheet = 2)
grd$uid <- paste0("x", gsub("-", "_", grd$CAS.Number))
rownames(grd) <- grd$uid
grd <- grd[grd$In.Batch.2 == "Yes", -c(2, 4, 6, 7, 8, 9) ]
pd$Group[pd$Group=="Proband"] <- "Case"
pd$Group <- factor(pd$Group, levels = c("Control", "Father", "Mother", "Case"))
colnames(grd) <- c("Chemical", "Group", "Classification", "Origin", "uid")

gc_data <- list(
  colour = colorExposures,
  measures = dta,
  table = grd,
  pd = pd,
  fd = dsc
)
rm(dsc, dta, grd, pd, tdta, colorExposures)

## LOAD CORRELATIONS ##########################################################
load("data/GC-MS/paired_correlations.RData")
gc_crr <- all_paired
colnames(gc_crr) <- c("Exposure 1", "Exposure 2", "Correlation", 
  "Group", "Metabolites")
gc_crr$Correlation <- as.numeric(as.character(gc_crr$Correlation))
rm(all_paired, eF)


## PIE CHART ##################################################################
pieGC <- function() {
  dta <- as.data.frame(table(gc_data$fd$Proposed.Group))
  plot_ly(dta, labels = ~Var1, values = ~Freq, type = 'pie')
}

## FAMILY ARRAY ###############################################################
familcyGC <- function() {
  fm <- unique(gc_data$fd$Proposed.Group)
  fm <- fm[order(fm)]
  fm
}

## NAMES ARRAY ###############################################################
namesGC <- function(group) {
  sel <- gc_data$fd$Chemical.Name[gc_data$fd$Proposed.Group == group]
  sel
}

## BOXPLOT ####################################################################
boxplotGC <- function(metab) {
  dta <- gc_data$fd$uid[gc_data$fd$Chemical.Name == metab]
  dta <- as.data.frame(gc_data$measures[ , dta, drop = FALSE])
  dta$Group <- as.character(gc_data$pd[rownames(dta), c("Group")])
  colnames(dta) <- c("ppm", "group")
  plot_ly(dta, y = ~ppm, color = ~group, type = "box") %>%
    layout(
      title = metab,
      xaxis = list(title = ""), 
      yaxis = list(title = "ppm")
    )
}

## CORRELATION GRID ###########################################################
correlationGC <- function() {
  plt <- ggplot(gc_crr, 
      aes_string(x="Metabolites", y="Group", fill="Correlation")) + 
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("") + ylab("") +
    scale_fill_gradient2(name="Correlation", low="red", mid="white", high="blue")
  ggplotly(plt)
}
