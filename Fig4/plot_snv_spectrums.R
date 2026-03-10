# Modified from code written by Sarah Moody, Iñigo Martincorena and Tim Coorens
library(RColorBrewer)
library(stringr)

mutation_matrix_LCM <- read.table("../data/LCM_SBS96_matrix_corrected_merged.txt", header = TRUE, sep = "\t",check.names = F)
mutation_matrix_bulk <- read.table("../data/bulk_SBS96_matrix_corrected.txt", header = TRUE, sep = "\t",check.names = F)
mutation_matrix_cancer <- read.table("../data/RCC_SBS96_matrix.txt", header = TRUE, sep = "\t",check.names = F)
mutation_matrix_blood <- read.table("../data/blood_SBS96_matrix_corrected.txt", header = TRUE, sep = "\t",check.names = F)

colnames(mutation_matrix_LCM) <- gsub(
  "_Proximal_tubules$",
  " - Proximal tubules",
  colnames(mutation_matrix_LCM)
)
colnames(mutation_matrix_LCM) <- gsub(
  "_Distal_tubules$",
  " - Distal tubules",
  colnames(mutation_matrix_LCM)
)
colnames(mutation_matrix_LCM) <- gsub(
  "_Medulla$",
  " - Medulla",
  colnames(mutation_matrix_LCM)
)
colnames(mutation_matrix_LCM) <- gsub(
  "_Glomeruli$",
  " - Glomeruli",
  colnames(mutation_matrix_LCM)
)

colnames(mutation_matrix_cancer) <- str_c(colnames(mutation_matrix_cancer), " - Tumor")
colnames(mutation_matrix_blood) <- str_c(colnames(mutation_matrix_blood), " - Blood")
colnames(mutation_matrix_bulk) <- str_c(colnames(mutation_matrix_bulk), " - Bulk cortex")

mutation_matrix <- cbind(mutation_matrix_LCM,mutation_matrix_cancer[,-1], mutation_matrix_blood[,-1])
write.csv(mutation_matrix[,grepl("MutationType|PD46868|PD42887|PD67546",colnames(mutation_matrix))],
          "../../data_S1/Fig4_SNV_spectrums1.csv", row.names = FALSE)

# ------------plotting mutational spectra for different individuals----------
# choose one of the four individuals
# PD46868
mutation_matrix<- mutation_matrix[,grepl("MutationType|PD46868",colnames(mutation_matrix))]
# PD42887
mutation_matrix<- mutation_matrix[,grepl("MutationType|PD42887",colnames(mutation_matrix))]
# PD67546
mutation_matrix<- mutation_matrix[,grepl("MutationType|PD67546",colnames(mutation_matrix))]
# PD55516
mutation_matrix <- cbind(mutation_matrix_bulk,mutation_matrix_cancer[,-1], mutation_matrix_blood[,-1])
mutation_matrix <- mutation_matrix[,grepl("MutationType|PD55516",colnames(mutation_matrix))]
write.csv(mutation_matrix,
          "../../data_S1/Fig4_SNV_spectrums2.csv", row.names = FALSE)
#Prepare and reorder the matrix
mutation_matrix$Ref <- substr(mutation_matrix$MutationType, 3, 3)
mutation_matrix$Mut <- substr(mutation_matrix$MutationType, 5, 5)
mutation_matrix <- mutation_matrix[order(mutation_matrix$Ref, mutation_matrix$Mut),]
SNV_data <- mutation_matrix[,sapply(mutation_matrix, is.numeric)]


plot_SNV_spectrums = function(SNV_data){
  # Specify Context Type
  sig_cat = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sig_cat,each=16),rep(ctx_vec,times=6),sep=",")
  snv_context = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  # Specify Vectors for plot colours
  col_vec_num <- rep(16,6)
  col_vec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  
  # Convert to matrix
  sig_plot <- as.matrix(SNV_data)
  
  # Set up Signature Names and title
  plot_title <- title
  sample_names <- colnames(sig_plot)
  
  # Set up Counts
  sample_counts <- colSums(sig_plot)
  
  # Plot signatures
  for (i in 1:ncol(sig_plot)) {
    # set up output files
    pdf(paste0(sample_names[i],"_SNV_spectrum.pdf"),width=10,height=3.5)
    par(xaxs='i', xpd = FALSE, mgp = c(3, 0.2, 0))
    #  define maxy
    max_prob <- sig_plot[,i];maxy = 300 ## Set y scale for each individual. PD42887:1500, PD46868:150, PD67546:270, PD55516:270
    #  call empty plot
    b = barplot(sig_plot[,i], col = NA, border = NA, axes = FALSE, las = 2, 
                ylim=c(0,1.5*maxy), names.arg = NA, cex.lab = 1.5, 
                cex.names = 0.6, cex.axis = 1.5, space = 1, ylab = "Mutations", 
                xaxs='i')  
    # add gridlines
    abline_pos = pretty(0:(1.5*maxy), n = 3)
    abline(h = abline_pos[2], col = 'grey90')
    abline(h = abline_pos[3], col = 'grey90')
    abline(h = abline_pos[4], col = 'grey90')
    # add axis
    axis(2, at = pretty(0:(1.5*maxy), n = 3), col = NA, las = 1, cex.axis = 1.5)
    # call columns
    b = barplot(as.numeric(sig_plot[,i]), axes = FALSE, col=col_vec, add = T, border = NA, space = 1)
    # add box surronding plot
    box(lty = 1, col = 'grey90')
    # add mutation count
    title(paste0(sample_names[i]), line = 2.5, adj = 0.5, cex.main = 2)
    title(paste0("Total SNV: ", sample_counts[i]), line = -1, adj = 0.99, cex.main = 1.5)
    # add rectangles and annotations on top of the plot
    par(xpd = NA)
    for (j in 1:length(sig_cat)) {
      xpos = b[c(sum(col_vec_num[1:j])-col_vec_num[j]+1,sum(col_vec_num[1:j]))]
      rect(xpos[1]-0.5, maxy*1.5, xpos[2]+0.5, maxy*1.6, border=NA, col=unique(col_vec)[j])
      text(x=mean(xpos), pos=3, y=maxy*1.6, label=sig_cat[j], cex = 1.5)
    } 
    par(xpd = FALSE)
    dev.off()
    
  }
}

setwd("../figure/Fig4/SBS12")
plot_SNV_spectrums(SNV_data)
