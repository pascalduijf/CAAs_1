#dependencies
library(ggplot2)
library(gridExtra)
library(plotly)
setwd("~/Desktop/Volcano_All_plots")
diff_df <- read.delim(file = "~/Desktop/Volcano_All_plots/MESO_type_cooccurence_stats_2017-08-10/MESO_type_cooccurence_stats_2017-08-10.txt" ,header = TRUE,sep = ',')
colnames(diff_df)
diff_df <- diff_df[c("sp_name", "Frequency_of_real", "p_value","X.log.pval.","obs_cooccur","color")]
#diff_df["group"] <- "NotSignificant"
diff_df[which(diff_df['p_value'] < .05),"group"] <- "Significant"
diff_df[which(diff_df['p_value'] > .05),"group"] <- "Not-Significant"
top_peaks <- diff_df[with(diff_df, order(-obs_cooccur,p_value)),][0:0,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-obs_cooccur,p_value)),][0:0,])
a <- list()
for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
        x = m[["Frequency_of_real"]],
        y = m[["X.log.pval."]],
        text = m[["sp_name"]],
        xref = "x",
        yref = "y",
        showarrow = FALSE,
      arrowhead = 0.1,
       ax = 15,
        ay = -40
    )
}

f <- list(
    family = "Arial, monospace",
    size = 18,
    color = "Black")

x <- list(
    title = "Frequency (%) of co-occurring events",
    titlefont = f
)
y <- list(
    title = "-log10(P Value)",
    titlefont = f
)
p <- plot_ly(data = diff_df, colors=c('orange','blue'),x = diff_df$Frequency_of_real, y = diff_df$X.log.pval., text = diff_df$sp_name, mode = "markers", color = diff_df$color) %>% layout(xaxis = x, yaxis= y,title ="Cooccuring CAL-SCNAs in Mesothelioma ; Total Cases=87") %>% layout(annotations = a)
p
