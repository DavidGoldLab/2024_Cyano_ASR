library(phytools)
library(corHMM)

# Load the data
# Copied and renamed tree ./1_Tree_Building/CyanosEnvironment.fa.contree > 2_Bacteria_213.tree
# Converted 2_Bacteria_213_Ancestral_States.xlsx into txt (replaced all "-" characters with "_" to match tree)
tree <- read.tree(file = "2-Bacteria_213.tree")
character_matrix <- read.table("2-Bacteria_213_Ancestral_States.txt", header = TRUE, row.names = 1, sep = "\t")

#####
# QC
#####

# Check if the number of tips in the tree and rows in the discrete trait matrix match
if (Ntip(tree) != nrow(character_matrix)) {
  stop("Number of tips in the tree and rows in the discrete trait matrix do not match.")
}

# Check if tree IDs match IDs in text file
tree_ids <- tree$tip.label
text_ids <- rownames(character_matrix)
matching_ids <- intersect(tree_ids, text_ids)
missing_in_text <- setdiff(tree_ids, text_ids)
missing_in_tree <- setdiff(text_ids, tree_ids)
cat("\nIDs in tree file but not in text file:\n")
print(missing_in_text)
cat("\nIDs in text file but not in tree file:\n")
print(missing_in_tree)

#####################################################
# Make PDF with node labels (for interpreting output)
#####################################################

# Find nodes of interest
pdf("2-Bacteria_213.tree_with_nodes.pdf", width = 8, height = 40)
# Plot the tree
plotTree(tree, ftype="i", fsize=0.7, lwd=1)
# Add node labels
nodelabels(bg="white", cex=0.5)
# Close the PDF device
dev.off()

############################
# Perform ASR (phytools)
############################

# List of column names
column_names <- c("Air","brackish","Cultivated","Food","freshwater",
	"host_associated","hydrothermal","marine","Plant_Associated",
	"soil","unknown","wastewater","other")

# ASR for ER model

for (column_name in column_names) {
  # Extract the specified column as a numeric vector while preserving row names
  modified_vector <- as.numeric(character_matrix[, column_name])
  # Check for non-numeric elements or missing values
  if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
    stop("Vector contains non-numeric or missing values.")
  }
  # Assign row names to the modified vector
  names(modified_vector) <- rownames(character_matrix)
  # Perform stochastic character mapping (ER Model)
  sim_result <- make.simmap(tree, modified_vector, model = "ER", nsim = 1000)
  # Export density map as PDF
  pdf_name <- paste(column_name, "-ER_Model-Density_Map.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  density_map <- densityMap(sim_result, fsize = 0.5)
  dev.off()
  # Export node summary as PDF
  pdf_name <- paste(column_name, "-ER_Model-Node_Summary.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  summary <- plot(summary(sim_result), fsize = 0.5)
  dev.off()
  # summarize the results
  summary <- describe.simmap(sim_result)
  # See posterior probabilities of nodes
  txt_name <- paste(column_name, "-ER_Model-Summary.txt", sep = "")
  capture.output(print(summary[["ace"]]), file = txt_name)
}

# Repeat for ARD model

for (column_name in column_names) {
  # Extract the specified column as a numeric vector while preserving row names
  modified_vector <- as.numeric(character_matrix[, column_name])
  # Check for non-numeric elements or missing values
  if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
    stop("Vector contains non-numeric or missing values.")
  }
  # Assign row names to the modified vector
  names(modified_vector) <- rownames(character_matrix)
  # Perform stochastic character mapping (ER Model)
  sim_result <- make.simmap(tree, modified_vector, model = "ARD", pi="estimated", nsim = 1000)
  # Export density map as PDF
  pdf_name <- paste(column_name, "-ARD_Model-Density_Map.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  density_map <- densityMap(sim_result, fsize = 0.5)
  dev.off()
  # Export node summary as PDF
  pdf_name <- paste(column_name, "-ARD_Model-Node_Summary.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  summary <- plot(summary(sim_result), fsize = 0.5)
  dev.off()
  # summarize the results
  summary <- describe.simmap(sim_result)
  # See posterior probabilities of nodes
  txt_name <- paste(column_name, "-ARD_Model-Summary.txt", sep = "")
  capture.output(print(summary[["ace"]]), file = txt_name)
}

###############################################
# Perform ASR with Hidden Rate Models (CorHMM)
###############################################

# ASR HRM for ER model

for (column_name in column_names) {
  # Extract the specified column as a vector while preserving row names
  modified_vector <- character_matrix[, column_name, drop = FALSE]
  # Convert the specified column to a numeric vector
  modified_vector <- as.numeric(character_matrix[, column_name])
  # Check for non-numeric elements or missing values
  if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
    stop("Vector contains non-numeric or missing values.")
  }
  # Create a data frame with species names and trait values
  data <- data.frame(species = rownames(character_matrix), trait = modified_vector)
  # Fit the ARD Model with Hidden Rates
  result <- corHMM(tree, data, rate.cat = 2, model = "ER")
  # Visualize the ancestral state reconstruction
  pdf_name <- paste(column_name, "-ER-HRM_Model-Node_Summary.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  #  Plot density map
  plotRECON(tree, result$states, piecolors=NULL, cex=0.5, pie.cex=0.25, height=11, width=8.5, show.tip.label=TRUE)
  # Close the PDF device
  dev.off()
  # See posterior probabilities of nodes
  txt_name <- paste(column_name, "-ER-HRM_Model-Summary.txt", sep = "")
  capture.output(print(result$states), file = txt_name)
}

# ASR HRM for ARD model

for (column_name in column_names) {
  # Extract the specified column as a vector while preserving row names
  modified_vector <- character_matrix[, column_name, drop = FALSE]
  # Convert the specified column to a numeric vector
  modified_vector <- as.numeric(character_matrix[, column_name])
  # Check for non-numeric elements or missing values
  if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
    stop("Vector contains non-numeric or missing values.")
  }
  # Create a data frame with species names and trait values
  data <- data.frame(species = rownames(character_matrix), trait = modified_vector)
  # Fit the ARD Model with Hidden Rates
  result <- corHMM(tree, data, rate.cat = 2, model = "ARD")
  # Visualize the ancestral state reconstruction
  pdf_name <- paste(column_name, "-ARD-HRM_Model-Node_Summary.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  #  Plot density map
  plotRECON(tree, result$states, piecolors=NULL, cex=0.5, pie.cex=0.25, height=11, width=8.5, show.tip.label=TRUE)
  # Close the PDF device
  dev.off()
  # See posterior probabilities of nodes
  txt_name <- paste(column_name, "-ARD-HRM_Model-Summary.txt", sep = "")
  capture.output(print(result$states), file = txt_name)
}


#######################
# Model fit comparison
#######################

## Compare ER and ARD, with and without hidden rate models)

for (column_name in column_names) {
  # Extract the specified column as a numeric vector while preserving row names
  modified_vector <- as.numeric(character_matrix[, column_name])
  # Check for non-numeric elements or missing values
  if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
    stop("Vector contains non-numeric or missing values.")
  }
  # Assign row names to the modified vector
  names(modified_vector) <- rownames(character_matrix)
  ## ASR using phytools::ancr (Model = ER)
  er_model<-fitMk(tree,modified_vector,model="ER")
  er_ancr<-ancr(er_model)
  er_ancr
  ## ASR using phytools::ancr (Model = HRM ER)
  er_hrm<-fitHRM(tree,modified_vector,model="ER", parallel=TRUE)
  er_hrm_ancr<-ancr(er_hrm)
  er_hrm_ancr
  ## ASR using phytools::ancr (Model = ARD)
  ard_model<-fitMk(tree,modified_vector,model="ARD")
  ard_ancr<-ancr(ard_model)
  ard_ancr
  ## ASR using phytools::ancr (Model = HRM ARD)
  ard_hrm<-fitHRM(tree,modified_vector,model="ARD", parallel=TRUE)
  ard_hrm_ancr<-ancr(ard_hrm)
  ard_hrm_ancr
  # Compare models
  anova_result <- anova(er_model,er_hrm,ard_model,ard_hrm)
  # Capture the output
  anova_output <- capture.output(anova_result)
  # Specify the file name
  file_name <- paste(column_name, "-ANOVA-Model_Test.txt", sep = "")
  # Write the output to the text file
  writeLines(anova_output, con = file_name)
  # Confirm the file has been written
  cat("ANOVA results have been written to", file_name)
}

###############
# Session info
###############
sessionInfo()

# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.7
# 
# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: US/Pacific
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] phytools_2.3-0 maps_3.4.2     ape_5.8       
# 
# loaded via a namespace (and not attached):
#  [1] doParallel_1.0.17       cli_3.6.3               nlme_3.1-164           
#  [4] rlang_1.1.4             generics_0.1.3          clusterGeneration_1.3.8
#  [7] scatterplot3d_0.3-44    quadprog_1.5-8          grid_4.4.1             
# [10] expm_1.0-0              MASS_7.3-60.2           lifecycle_1.0.4        
# [13] foreach_1.5.2           numDeriv_2016.8-1.1     compiler_4.4.1         
# [16] igraph_2.0.3            codetools_0.2-20        coda_0.19-4.1          
# [19] pkgconfig_2.0.3         Rcpp_1.0.13             DEoptim_2.2-8          
# [22] fastmatch_1.1-4         optimParallel_1.0-2     phangorn_2.12.1        
# [25] lattice_0.22-6          digest_0.6.37           mnormt_2.1.1           
# [28] parallel_4.4.1          magrittr_2.0.3          Matrix_1.7-0           
# [31] tools_4.4.1             iterators_1.0.14        combinat_0.0-8         

