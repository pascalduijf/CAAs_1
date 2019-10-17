library(cooccur)
library(igraph)
library(qgraph)
library(scales) # rescale()
library(Hmisc)  # %nin%

setwd("/Users/lfink/projects/duijf_chromothripsis/chr_changes/")
datafile <- "/Users/lfink/projects/duijf_chromothripsis/chr_changes/Network-A.txt"
data     <- read.csv(datafile, header=TRUE, sep="\t") # get some basic information before trying to make plots


################ DRAW NETWORKS
## FUNCTION:
# read datafile with events, select a group using the "Origin/Type/Blood_solid" name, provide a name to use in file names and plot titles
draw_chr_network <- function(datafile, groupname, printname){
  #datafile      <- "/Users/lfink/projects/duijf_chromothripsis/chr_changes/Network-A.txt"
  #groupname     <- "Solid"
  #groupname     <- "BRCA"
  #printname     <- groupname
  pvalue        <- 0.05 # min p-value; probably never needed as most p-values close to 0
  #max_events    <- 75 # cap number of events to plot; otherwise network is too big
  max_events    <- 100
  
  
  data            <- read.csv(datafile, header=TRUE, sep="\t")
  if(groupname %in% c("Epithelial", "Ectoderm", "Neural ectoderm")) { 
    data            <- subset(data, Origin==groupname) ###  (Epithelial, etc)
    filetag         <- "origin"
  }
  if(groupname %in% c("Blood", "Solid")) {
    data            <- subset(data, Blood_solid==groupname) ###  (Blood, Solid)
    filetag         <- "bs"    
  }
  if(groupname %nin% c("Epithelial", "Ectoderm", "Neural ectoderm", "Blood", "Solid")) {
    data            <- subset(data, Type==groupname) ###  (ACC, BRCA, etc)
    filetag         <- "type"
  }
  
  # reformat data frame; remove unwanted observations
  samples          <- data$Sample
  data             <- data[,1:92] # remove X chromosome changes
  data             <- subset(data, select=-c(Sample, Type, Origin, Blood_solid))
  tdata            <- t(data)
  colnames(tdata)  <- samples
  n                <- length(samples)
  
  # count number of occurrences of each gain or loss in the desired subset of data
  counts           <- as.data.frame(colSums(data))
  colnames(counts) <- c("counts")
  # in rownames, replace "gain" with "+" and "loss" with "-" 
  rownames(counts) <- gsub("gain", "+", rownames(counts))
  rownames(counts) <- gsub("loss", "-", rownames(counts))
  counts$event     <- rownames(counts)
  
  #obs.v.exp(cooccur.chr)
  #pair.profile(cooccur.chr)
  cooccur.chr      <- cooccur(mat=tdata, type="spp_site", thresh=TRUE, spp_names=TRUE)
  
  plot(cooccur.chr)
  prob.table(cooccur.chr)
  pair.attributes(cooccur.chr) # can see these on the plot
  pair(cooccur.chr, spp = "gain1p") # looking at details for first one "gain1p"
  # can see that it has one significant negative cooccurrence with "loss9p"
  # can see positive cooccurrences with 6 others
  
  # print cooccurrence heatmap to PDF
  pdfname          <- gsub(" ", "_", printname)
  pdffile          <- paste0(pdfname, "_", filetag, "_cooccurence_matrix_", format(Sys.Date(), format="%F"), ".pdf")
  pdf(file=pdffile, height=12, width=12)
  print(plot(cooccur.chr))
  dev.off()
  
  
  ## FORMAT CO-OCCURRENCE RESULTS
  mat             <- cooccur.chr$results
  
  # remove impossible events (loss1p/gain1p)
  mat             <- subset(mat,     gsub("^([[:alpha:]]+)", "", mat$sp1_name) != gsub("^([[:alpha:]]+)",     "", mat$sp2_name))
  # replace "gain" with "+" and "loss" with "-"
  mat$sp1_name    <- gsub("gain", "+", mat$sp1_name)
  mat$sp1_name    <- gsub("loss", "-", mat$sp1_name)
  mat$sp2_name    <- gsub("gain", "+", mat$sp2_name)
  mat$sp2_name    <- gsub("loss", "-", mat$sp2_name)
  
  # create edges based on log2(obs/exp)
  mat$obs_cooccur[which(mat$obs_cooccur == 0)] <- 0.00001 #handle occurrences of 0 (make it 0.00001 for a pseudocount)
  mat$log2oVe     <- log2(mat$obs_cooccur / mat$exp_cooccur)
  
  # assign p_lt or p_gt to p-value based over- or under-expected observations; sort by p-value
  mat$obs_v_exp   <- mat$obs_cooccur / mat$exp_cooccur
  mat$p_value     <- mat$p_lt
  mat$p_value[which(mat$obs_v_exp > 1)] <- mat$p_gt[which(mat$obs_v_exp > 1)]
  mat             <- mat[with(mat, order(p_value)), ]
  mat             <- mat[with(mat, order(p_value,-abs(mat$log2oVe))), ]

  # remove unnecessary columns
  #mat             <- subset(mat, select=c(sp1_name, sp2_name, p_lt, obs_cooccur, log2oVe))  
  mat             <- subset(mat, select=c(sp1_name, sp2_name, log2oVe, p_value, obs_cooccur, exp_cooccur))  
  
  # remove all events with "non-significant" p-values
  mat             <- subset(mat, p_value < pvalue)
  
  # order by most frequenct occurrences
  #mat             <- mat[with(mat, order(-obs_cooccur, p_value)), ]
  mat             <- mat[with(mat, order(-obs_cooccur, -abs(log2oVe))), ]
  

  # write cooccurrence stats to file
  datname         <- gsub(" ", "_", printname)
  datfile         <- paste0(datname, "_", filetag, "_cooccurence_stats_", format(Sys.Date(), format="%F"), ".txt")
  write.table(mat, datfile, col.names = TRUE, row.names = FALSE, quote=FALSE)
  

  ## DRAW NETWORKS
  # We used this line in the plot function -->  comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
  
  # just select the top X significant/big effect events or network gets too tangled
  num_events      <- dim(mat)[1]
  num_include     <- min(max_events, num_events)
  mat             <- mat[1:num_include,] 
  # for printing to plot; show highest p-value used
  max_pval_in_set <- max(mat$p_value)
  #mat             <- subset(mat, select=c(sp1_name, sp2_name, log2oVe, p_value))
  
  #effs <- effect.sizes(mod = cooccur.chr, standardized = TRUE, matrix = FALSE) # DMG example
  #hist(effs$effects)
  #obs.v.exp(cooccur.chr)
   
  #mat             <- subset(mat, select=c(sp1_name, sp2_name, log2oVe, obs_cooccur)) 
  mat             <- subset(mat, select=c(sp1_name, sp2_name, log2oVe)) # DMG but if we go with the log ratios as edges
  Q               <- qgraph(mat)
  nodecolor       <- vector()
  nodecolor[which(gsub("^(\\+)(\\w+)", "\\1", Q$graphAttributes$Nodes$names) == "+")] <- "green4"
  nodecolor[which(gsub("^(-)(\\w+)",   "\\1", Q$graphAttributes$Nodes$names) == "-")] <- "red3"
  # get size of each node - based on frequency of that gain or loss in the data set
  nodesize        <- counts[which(counts$event %in% Q$graphAttributes$Nodes$names),]
  nodesize        <- nodesize[Q$graphAttributes$Nodes$names,] # order them to match the Nodes
  # scale nodesize so min and max node size is consistent across all plots
  # One thing to think about is the fact that we usually calculate effects sizes as the difference in the expected and 
  # observed co-occurrences and then we usually standardize them to the number of "sites" available
  scalednodesize <- rescale(nodesize$counts, c(3, 6))
  
  # draw network
  pdfname         <- gsub(" ", "_", printname)
  pdffile         <- paste(pdfname, "_", filetag, "_network_", format(Sys.Date(), format="%F"), ".pdf", sep="")
  pdf(file=pdffile, height=12, width=12)
  Q               <- qgraph(mat, label.cex=0.9, color=nodecolor, posCol = "blue", negCol = "darkorange", arrows=FALSE,
                            label.color="white",
                            vsize=scalednodesize 
                            )
  title(paste0("Top ", max_events, " Co-occurring Chromosome Gains and Losses in ", printname, "; n=", n, "; max p-value=", max_pval_in_set), cex=1.5, line = 2.5)
  dev.off()
}
###


# DRAW NETWORKS FOR EACH CANCER (RESET FUNCTION PARAMETERS ACCORDINGLY)
cancer_types    <- unique(data$Type)
cancer_types    <- as.character(cancer_types)
for (i in 1:length(cancer_types)) {
#for (i in 1:2) {
  draw_chr_network(datafile, cancer_types[i], cancer_types[i])
}

# DRAW NETWORKS FOR EACH GROUP (RESET FUNCTION PARAMETERS ACCORDINGLY)
#data            <- read.csv(datafile, header=TRUE, sep="\t")
cancer_groups   <- c("Epithelial", "Ectoderm", "Neural ectoderm")
cancer_groupspr <- gsub(" ", "_", cancer_groups)
for (i in 1:length(cancer_groups)) {
#for (i in 1:2) {
  draw_chr_network(datafile, cancer_groups[i], cancer_groupspr[i])
}

# DRAW NETWORKS FOR EACH BLOOD_SOLID (RESET FUNCTION PARAMETERS ACCORDINGLY)
#  for solid: Of 3828 species pair combinations, 3811 pairs (99.56 %) were removed from the analysis because expected co-occurrence was < 1 and 17 pairs were analyzed
#data            <- read.csv(datafile, header=TRUE, sep="\t")
#cancer_groups   <- c("Blood", "Solid")
draw_chr_network(datafile, "Blood", "Blood")
draw_chr_network(datafile, "Solid", "Solid")
