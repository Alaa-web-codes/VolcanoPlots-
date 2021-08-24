library(ggplot2)
library(dplyr)

#===================
#===================
#**Modify**
DIR <- "data/" #where your input file is
file <- "Protein data CT and ejection .csv" #name of input file

#test reading
#data <- read.csv(paste0(DIR,file), stringsAsFactors = FALSE)
#head(data)


#**Define columns for genes and two conditions**
#*
gene_col <- 1 #column where my proteins are

cond1_colstart <-2 #column where condition 1 starts
cond1_colend <- 37 #column where condition 1 ends

cond2_colstart <- 38 #column where condition 2 starts
cond2_colend <- 50 #column where condition 2 ends
#===================
#===================

#**Import data**
data <- read.csv(paste0(DIR,file), stringsAsFactors = FALSE)
View(data) #show data

#**Calculate p-value**
pval <- vector(mode="numeric")

for (i in 1:nrow(data)) {
  
  cond1 <- unlist(data[i,cond1_colstart:cond1_colend])
  cond2 <- unlist(data[i,cond2_colstart:cond2_colend])
  
  result <- t.test(cond1,cond2)$p.value
  
  pval <- append(pval,result)
  
}


#**Calculate foldchange**
mean_cond1 <- vector(mode="numeric")
mean_cond2 <- vector(mode="numeric")

for (i in 1:nrow(data)) {
  
  cond1 <- unlist(data[i,cond1_colstart:cond1_colend])
  cond2 <- unlist(data[i,cond2_colstart:cond2_colend])
  
  result1 <- mean(cond1)
  result2 <- mean(cond2)
  
  mean_cond1 <- append(mean_cond1,result1)
  mean_cond2 <- append(mean_cond2,result2)
  
}

foldchange <- mean_cond1/mean_cond2 #cond1 = celltherapy #cond2= Rejection


#**Log scale**

neg_log10_pval <- -log(pval,10)
log10_foldchange <- log(foldchange,10)


#**Combine to dataframe**
results <- cbind(data[,gene_col,drop=FALSE],pval, foldchange, neg_log10_pval,log10_foldchange)
View(results) #show results


#**Volcano Plot**
#============
#**Modify**
y_label <- "-log10 (p-value) \n" # \n - this adds a break, to have more space between axis label and plot
x_label <- "\n log10 fold change"
length_x_axes=c(-1.3,1.3) #x-axes limits - check your data for limits - results - log10_foldchange
length_y_axes=c(0,12) #y-axes limits - check your data for limits - results - neg_log10_pval
sig_level=0.05 #significance level - for horizontal line in plot
FC_cutoff <- 3 #foldchange cutoff - for vertical lines in plot
#===========

horizontal_line=-log10(sig_level) #calculate value for horizontal line
FC <- log(FC_cutoff,10) #calculate value for vertical lines

VolcanoPlot <- ggplot(data=results, 
                      aes_string(x=log10_foldchange,
                                 y=neg_log10_pval)) +
  geom_hline(yintercept = horizontal_line,linetype="dashed",color="darkgrey",size=0.5)+ # horizontal line
  geom_vline(xintercept= FC, linetype="dashed",color="darkgrey",size=0.5)+ #vertical line
  geom_vline(xintercept= -FC, linetype="dashed",color="darkgrey",size=0.5)+  #vertical line
  geom_point(shape=21, fill="lightgrey", color="black",size=2)+  
  xlab(x_label)+ #x axes label
  ylab(y_label)+ #y axes label
  xlim(length_x_axes)+ # x axes limits
  scale_y_continuous(expand = c(0, 0), limits = length_y_axes)+ #expand forces the y axes to start at 0, limits refers to your y axes limits
  theme(axis.line.x = element_line(color="black", size = 1), #x axes line
        axis.line.y = element_line(color="black", size = 1))+ #y axes line
  theme(axis.text=element_text(size=10), # size of axes texts
        axis.title=element_text(size=12)) # size of axes titles

plot(VolcanoPlot) #show Volcano Plot


#**Export**
#export result dataframe --------------------------------------
#--------------------------------------------------------------
write.csv(results, paste0(DIR,"results1.csv"),row.names = FALSE) 

#export plot as jpg
#Note that you can also export plots manually via RStudio - also adjust "size" elements in the ggplot according to your chosen height/width
#1. Make jpeg
jpeg(paste0(DIR,"VolcanoPlot.jpg"),width = 3.5, height = 3.5, units = 'in', res = 300)
# 2. Create the plot
plot(VolcanoPlot)
# 3. Close the file
dev.off()
