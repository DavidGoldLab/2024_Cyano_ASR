library(phytools)
packageVersion("phytools")
# 2.0.13

# Load the data
# Copied and renamed tree ./1_Tree_Building/CyanosEnvironment.fa.contree > 2_Bacteria_213.tree
# Converted 2_Bacteria_213_Ancestral_States.xlsx into txt (replaced all "-" characters with "_" to match tree)
tree <- read.tree(file = "2_Bacteria_213.tree")
character_matrix <- read.table("2_Bacteria_213_Ancestral_States.txt", header = TRUE, row.names = 1, sep = "\t")

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


# ########################
# # Perform ASR (1 sample)
# ########################
# 
# # Specify the column name you want to include
# column_name <- "Air"
# # Extract the specified column as a vector while preserving row names
# modified_vector <- character_matrix[, column_name, drop = FALSE]
# # Extract the specified column as a numeric vector while preserving row names
# modified_vector <- as.numeric(character_matrix[, column_name])
# # Check for non-numeric elements or missing values
# if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
#   stop("Vector contains non-numeric or missing values.")
# }
# # Assign row names to the modified vector
# names(modified_vector) <- rownames(character_matrix)
# # Perform stochastic character mapping
# sim_result <- make.simmap(tree, modified_vector, method = "ER", nsim = 1000)
# # Create PDF File name
# pdf_name <- paste("density_map_simulation_", column_name, ".pdf", sep = "")
# pdf(pdf_name, width = 8, height = 10)
# #  Plot density map
# density_map <- densityMap(sim_result,fsize=0.5)
# # Close the PDF device
# dev.off()

############################
# Perform ASR (all samples)
############################

# List of column names
column_names <- c("Air","brackish","Cultivated","Food","freshwater",
	"host_associated","hydrothermal","marine","Plant_Associated",
	"soil","unknown","wastewater","other")

# Loop over each column
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
  pdf_name <- paste("density_map_simulation-", column_name, "-ER_Model.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  density_map <- densityMap(sim_result, fsize = 0.5)
  dev.off()
  # Export summary results as PDF
  pdf_name <- paste("node_summary-", column_name, "-ER_Model.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  summary <- plot(summary(sim_result), fsize = 0.5)
  dev.off()
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
  # Perform stochastic character mapping (ARD Model)
  sim_result <- make.simmap(tree, modified_vector, model = "ARD", pi="estimated", nsim = 1000)
  # Export density map as PDF
  pdf_name <- paste("density_map_simulation-", column_name, "-ARD_Model.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  density_map <- densityMap(sim_result, fsize = 0.5)
  dev.off()
  # Export summary results as PDF
  pdf_name <- paste("node_summary-", column_name, "-ARD_Model.pdf", sep = "")
  pdf(pdf_name, width = 8, height = 10)
  summary <- plot(summary(sim_result), fsize = 0.5)
  dev.off()
}

###############################################
# Get node summaries (used for select analyses)
###############################################
# summarize the results
summary <- describe.simmap(sim_result)
# See posterior probabilities of nodes
summary[["ace"]]

#######################
# Model fit comparison
#######################

## Compare ER and ARD, with and without hidden rate models)

## Prepare dataset ##
# Specify the column name you want to include
column_name <- "Freshwater"
# Extract the specified column as a vector while preserving row names
modified_vector <- character_matrix[, column_name, drop = FALSE]
# Extract the specified column as a numeric vector while preserving row names
modified_vector <- as.numeric(character_matrix[, column_name])
# Check for non-numeric elements or missing values
if (any(is.na(modified_vector)) || any(!is.finite(modified_vector))) {
  stop("Vector contains non-numeric or missing values.")
}
# Assign row names to the modified vector
names(modified_vector) <- rownames(character_matrix)

## ASR using phytools::ancr (Model = ER)
Freshwater.er_model<-fitMk(tree,modified_vector,model="ER")
Freshwater.er_ancr<-ancr(Freshwater.er_model)
Freshwater.er_ancr

## ASR using phytools::ancr (Model = HRM ER)
Freshwater.er_hrm<-fitHRM(tree,modified_vector,model="ER", parallel=TRUE)
Freshwater.er_hrm_ancr<-ancr(Freshwater.er_hrm)
Freshwater.er_hrm_ancr

## ASR using phytools::ancr (Model = ARD)
Freshwater.ard_model<-fitMk(tree,modified_vector,model="ARD")
Freshwater.ard_ancr<-ancr(Freshwater.ard_model)
Freshwater.ard_ancr

## ASR using phytools::ancr (Model = HRM ARD)
Freshwater.ard_hrm<-fitHRM(tree,modified_vector,model="ARD", parallel=TRUE)
Freshwater.ard_hrm_ancr<-ancr(Freshwater.ard_hrm)
Freshwater.ard_hrm_ancr

# Compare models
anova(Freshwater.er_model,Freshwater.er_hrm,Freshwater.ard_model,Freshwater.ard_hrm)

## Results ##
#                         log(L) d.f.      AIC     weight
# Freshwater.er_model  -54.79540    1 111.5908 0.09512106
# Freshwater.er_hrm    -50.68841    3 107.3768 0.78222047
# Freshwater.ard_model -53.81283    2 111.6257 0.09347721
# Freshwater.ard_hrm   -48.97702    8 113.9540 0.02918126

# Session info
sessionInfo
