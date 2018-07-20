rm(list = ls())

source("ba_functions.r")

data1 <- read.csv("ba_2007_dataset.csv", stringsAsFactors = F, header = T)
data1$diff <- data1$diff <- data1$RV - data1$IC

plot(RV ~ IC, data = data1,
     pch = as.character(data1$Subj.), col = "white")

text(x = data1$IC, y = data1$RV, labels = data1$Subj., adj = NULL,
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1, col = NULL, font = NULL)

results_table <- get_ba_result(data1$RV, data1$IC, data1$Subj., plotit = T, verbose = T)
require(knitr)
kable(results_table)
kable(standard_ba(data1$RV, data1$IC, plotit = T, verbose = T))

plot_sd(data1$RV, data1$IC, data1$Subj., plotit = T, verbose = T)
plot_subject_averages(data1$RV, data1$IC, data1$Subj., verbose = T)

