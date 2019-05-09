# R script to perform misc. data analyses for Mavrich & Hatfull, mBio, 2019.
# Travis Mavrich.
# Note: this code merges and analyzes various input data sets 
# prepared by other tools including Python and Excel.
# Distinct analyses and code blocks separated by "###".

### 1. Prepare environment.
### 2. Define functions.
### 3. Import datasets.
### 4. Average immunity data.
### 5. Compare variability between replicates.
### 6. Misc. general analyses of averaged immunity data.
### 7. Compare reciprocal infection assays.
### 8. Compare lysogen and cloned-repressor strains.
### 9. Compare L5 and L5-derivative phages.
### 10. Compare escape mutants and parent phages.
### 11. Compute immunity profile correlations.
### 12. Analyze misc. phage genome metrics.
### 13. Compare stoperator site data.


### 1. Prepare environment.

# Setup dependencies.

# The melt function of reshape2 package is needed to convert a matrix
# to unique-pair table.
# install.packages("reshape2")
library(reshape2)

# Stringdist is needed to compute hamming distance between strings.
# install.packages("stringdist")
library(stringdist)

# Set working directory variables specific to local directory structure.
DIR_INPUT = "~/scratch/immunity_analysis/input/"
DIR_OUTPUT = "~/scratch/immunity_analysis/output/"

#Set paths for all input files.
IMMUNITY_DATA_FILENAME = 
  paste(DIR_INPUT,
        "immunity_data.csv",
        sep="")

GENOMIC_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "genomic_distance_data.csv",
        sep="")

PHAGE_METADATA_FILENAME =
  paste(DIR_INPUT,
        "phage_metadata.csv",
        sep="")

REPRESSOR_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "repressor_336_distance_data.csv",
        sep="")

CAS4_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "cas4_311_distance_data.csv",
        sep="")

ENDOVII_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "endovii_306_distance_data.csv",
        sep="")

DNAPOL_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "dnapol_311_distance_data.csv",
        sep="")

PORTAL_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "portal_311_distance_data.csv",
        sep="")

STOPERATOR_PWM_DATA_FILENAME = 
  paste(DIR_INPUT,
        "stoperator_pwm_distances.csv",
        sep="")

RECIPROCAL_INFECTION_TABLE_FILENAME = 
  paste(DIR_INPUT,
        "reciprocal_infection_matrix.csv",
        sep="")

STOPERATOR_SITES_FILENAME = 
  paste(DIR_INPUT,
        "stoperator_site_predictions.csv",
        sep="")


setwd(DIR_OUTPUT)
###
###
###
###
###
###
###
###
###
### 2. Define functions.

# Compute comparison fields.
compute_comparisons <- function(table){

  table$subcluster_compare <-
    ifelse(table$defending_subcluster == table$challenging_subcluster,
           as.character(table$defending_subcluster),
           "different")
  
  table$source_compare <-
    ifelse(table$defending_source == table$challenging_source,
           as.character(table$defending_source),
           "different")
  
  table$temperate_empirical_compare <-
    ifelse(table$defending_temperate_empirical == 
             table$challenging_temperate_empirical,
           as.character(table$defending_temperate_empirical),
           "different")
  
  table$functional_repressor_compare <-
    ifelse(table$defending_repressor_functional == 
             table$challenging_repressor_functional,
           as.character(table$defending_repressor_functional),
           "different")
  
  table$lysogen_type_compare <-
    ifelse(table$defending_lysogen_type == 
             table$challenging_lysogen_type,
           as.character(table$defending_lysogen_type),
           "different")
  
  table$integrase_compare <-
    ifelse(table$defending_pham_integrase == table$challenging_pham_integrase,
           as.character(table$defending_pham_integrase),
           "different")
  
  table$parb_compare <-
    ifelse(table$defending_pham_parb == table$challenging_pham_parb,
           as.character(table$defending_pham_parb),
           "different")
  
  table$repressor_hth_compare <-
    stringdist(as.character(table$defending_repressor_hth),
               as.character(table$challenging_repressor_hth),
               method = "hamming")
  
  table$repressor_length_full_compare <-
    abs(table$defending_repressor_length_full - 
          table$challenging_repressor_length_full)

  table$repressor_length_nterm_compare <-
    abs(table$defending_repressor_length_nterm - 
          table$challenging_repressor_length_nterm)

  table$repressor_length_cterm_compare <-
    abs(table$defending_repressor_length_cterm - 
          table$challenging_repressor_length_cterm)
  
  table$gene_content_clade_compare <-
    ifelse(table$defending_gene_content_clade == 
             table$challenging_gene_content_clade,
           as.character(table$defending_gene_content_clade),
           "different")
  
  table$subcluster_compare <- 
    as.factor(table$subcluster_compare)
  table$source_compare <- 
    as.factor(table$source_compare)
  table$temperate_empirical_compare <-
    as.factor(table$temperate_empirical_compare)
  table$functional_repressor_compare <-
    as.factor(table$functional_repressor_compare)
  table$lysogen_type_compare <-
    as.factor(table$lysogen_type_compare)
  table$integrase_compare <- 
    as.factor(table$integrase_compare)
  table$parb_compare <- 
    as.factor(table$parb_compare)
  table$gene_content_clade_compare <-
    as.factor(table$gene_content_clade_compare)
  
  return(table)
}


match_lys_clone_data <- function(lysogen_data,clone_data){
  
  clone_match_columns <- c("defending_challenging",
                           "defending_phage",
                           "challenging_phage",
                           "averaged_score",
                           "nuc_dist",
                           "gcd",
                           "repressor_cterm_mafft_dist_uncorrected",
                           "stoperator_pwd_dist_euc",
                           "gene_content_clade_compare")
  
  lysY_reduced <- subset(lysogen_data,select = c(clone_match_columns))
  names(lysY_reduced) <- paste('lys_',
                               names(lysY_reduced),
                               sep="")
  
  cloneY_reduced <- subset(clone_data,select = c(clone_match_columns))
  names(cloneY_reduced) <- paste('clone_',
                                 names(cloneY_reduced),
                                 sep="")
  
  lys_clone_compare <- merge(lysY_reduced,cloneY_reduced,
                             by.x='lys_defending_challenging',
                             by.y='clone_defending_challenging')
  
  lys_clone_compare$lys_averaged_score <- 
    as.numeric(as.character(lys_clone_compare$lys_averaged_score))
  
  lys_clone_compare$clone_averaged_score <- 
    as.numeric(as.character(lys_clone_compare$clone_averaged_score))
  
  lys_clone_compare$score_diff <- 
    lys_clone_compare$clone_averaged_score - 
    lys_clone_compare$lys_averaged_score
  
  return(lys_clone_compare)
}


subset_interclade <- function(table,clade_column){
  table_subset <- subset(table,
                         table[,clade_column] == "different")
  return(table_subset)
}


subset_intraclade <- function(table,clade_column,clade){
  table_subset <- subset(table,
                         table[,clade_column] == clade)
  return(table_subset)
}


subset_same <- function(table,phage1,phage2){
  table_subset <- subset(table,
                         as.character(table[,phage1]) ==
                           as.character(table[,phage2]))
  return(table_subset)
}


subset_diff <- function(table,phage1,phage2){
  table_subset <- subset(table,
                         as.character(table[,phage1]) !=
                           as.character(table[,phage2]))
  return(table_subset)
}


compute_binned_freq1 <- function(table,bin_num){
  
  table$defending_phage <- factor(table$defending_phage)
  table$challenging_phage <- factor(table$challenging_phage)
  
  output_table <- data.frame(bin_num,
                             nlevels(table$defending_phage),
                             nlevels(table$challenging_phage),
                             nrow(table),
                             nrow(subset(table,
                                         table$subcluster_compare != 
                                           "different")),
                             nrow(subset(table,
                                         table$subcluster_compare == 
                                           "different")))
  
  names(output_table) <- c("bin",
                           "defending_freq",
                           "challenging_freq",
                           "total_num_assays",
                           "num_intra_subcluster_assays",
                           "num_inter_subcluster_assays")
  
  output_table$defending_percent <- 
    output_table$defending_freq / nlevels(clade2_env$defending_phage)
  
  output_table$challenging_percent <-
    output_table$challenging_freq / nlevels(clade2_env$challenging_phage)
  
  output_table$total_assays_percent <-
    output_table$total_num_assays / nrow(clade2_env)
  
  output_table$intra_subcluster_percent <-
    output_table$num_intra_subcluster_assays / 
    nrow(subset(clade2_env,
                clade2_env$subcluster_compare != "different"))
  
  output_table$inter_subcluster_percent <-
    output_table$num_inter_subcluster_assays / 
    nrow(subset(clade2_env,
                clade2_env$subcluster_compare == "different"))
  
  return(output_table)
}


compute_binned_freq2 <- function(table1,bin_num,table2){
  
  output_table <- data.frame(bin_num,nrow(table1))
  names(output_table) <- c("bin","freq")
  output_table$freq_percent <- output_table$freq / nrow(table2)
  return(output_table)
  
}


# Report number of data points in each table.
print_rows <- function(table_1,table_2,table_3){
  
  total_rows <- nrow(table_1) + nrow(table_2) + nrow(table_3)
  
  print(paste("Number of data points in first table: ",nrow(table_1)))
  print(paste("Number of data points in second table: ",nrow(table_2)))
  print(paste("Number of data points in third table: ",nrow(table_3)))
  print(paste("Number of total data points: ",total_rows))
  
}


# Compute correlation.
compute_correlation <- function(table,field1,field2){

    lm_result <- lm(table[,field1] ~ table[,field2])
  summary(lm_result)
  
}


# Plot to compare infection scores of two groups of phages.
plot_tricolor_scatter1 <- function(table1,
                                  table2,
                                  table3,
                                  x_data,
                                  y_data,
                                  x_range,
                                  y_range,
                                  filename){
  
  # For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  table3 <- subset(table3,select = c(x_data,y_data))
  table3 <- na.omit(table3)
  
  print_rows(table1,table2,table3)

  # Plot data.
  par(mar=c(4,8,8,4))
  plot(table1[,x_data],
       table1[,y_data],
       pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
       xlim=x_range,ylim=y_range,col="red")
  par(new=TRUE)
  plot(table2[,x_data],
       table2[,y_data],
       pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
       xlim=x_range,ylim=y_range,col="grey")
  par(new=TRUE)
  plot(table3[,x_data],
       table3[,y_data],
       pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
       xlim=x_range,ylim=y_range,col="black")
  abline(0,1)
  
  dev.copy(pdf,filename)
  dev.off()
  
}


# Plot to compare infection scores against a distance metric.
plot_tricolor_scatter2 <- function(table1,
                                   table2,
                                   table3,
                                   x_data,
                                   y_data,
                                   x_range,
                                   y_range,
                                   filename,
                                   abline_value = "no"){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  table3 <- subset(table3,select = c(x_data,y_data))
  table3 <- na.omit(table3)

  print_rows(table1,table2,table3)

  # Plot data.
  par(mar=c(4,8,16,4))
  plot(table1[,x_data],
       table1[,y_data],
       xlim=x_range,ylim=y_range,pch=16,cex=2,cex.axis=2,
       ann=FALSE,main=NULL,las=1,col="red")
  par(new=TRUE)
  plot(table2[,x_data],
       table2[,y_data],
       xlim=x_range,ylim=y_range,pch=16,cex=2,cex.axis=2,
       ann=FALSE,main=NULL,las=1,col="grey")
  par(new=TRUE)
  plot(table3[,x_data],
       table3[,y_data],
       xlim=x_range,ylim=y_range,pch=16,cex=2,cex.axis=2,
       ann=FALSE,main=NULL,las=1,col="black")
  
  if(abline_value == "yes"){
    abline(h=0)
  }

  dev.copy(pdf,filename)
  dev.off()
  
  compute_correlation(table1,x_data,y_data)

}


# Plot to compare genome metrics.
plot_bicolor_scatter1 <- function(table1,
                                  table2,
                                  x_data,
                                  y_data,
                                  x_range,
                                  y_range,
                                  filename,
                                  abline_value = "no"){
  
  # For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  total_rows <- nrow(table1) + nrow(table2)
  
  # For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))
  print(paste("Number of total data points: ",total_rows))
  
  # Plot data.
  par(mar=c(4,8,8,4))
  plot(table1[,x_data],
       table1[,y_data],
       xlim=x_range,ylim=y_range,pch=1,cex=1,cex.axis=2,
       ann=FALSE,main=NULL,las=1,col="grey")
  par(new=TRUE)
  plot(table2[,x_data],
       table2[,y_data],
       xlim=x_range,ylim=y_range,pch=1,cex=1,cex.axis=2,
       ann=FALSE,main=NULL,las=1,col="red")
  
  if(abline_value == "yes"){
    abline(0,1)
    }

  dev.copy(pdf,filename)
  dev.off()
}


plot_bargraph1 <- function(table1,value1,value2,y_range,filename){
  
  par(mar=c(4,8,4,4))
  barplot(table1[,value1],
          names.arg=table1[,value2],
          ylim=y_range,
          main=NULL,ann=FALSE,las=1,cex.axis=2,col="black")
  dev.copy(pdf,filename)
  dev.off()
}


plot_bargraph2 <- function(table1,value1,y_range,filename){

  par(mar=c(4,8,8,4))
  barplot(summary(table1[,value1]),
          ylim=y_range,
          main=NULL,ann=FALSE,las=1,cex.axis=2,col="black")
  dev.copy(pdf,filename)
  dev.off()
}


plot_hist1 <- function(table1,value1,num_breaks,x_range,y_range,filename){
  
  par(mar=c(4,8,8,4))
  hist(table1[,value1],
       breaks=num_breaks,
       xlim=x_range,
       ylim=y_range,
       main=NULL,ann=FALSE,las=1,cex.axis=2,col="black")
  dev.copy(pdf,filename)
  dev.off()
}
###
###
###
###
###
###
###
###
###
### 3. Import datasets.

# Import immunity data.
# Data structure:
# 1. immunity_assay_id (unique identifier for each specific assay)
# 2. immunity_set
# 3. date
# 4. notebook
# 5. page
# 6. strain (systematic strain name)
# 7. prophage
# 8. repressor_clone
# 9. strain_type (lysogen or repressor_clone)
# 10. defending_phage
# 11. challenging_phage
# 12. assay_type (multiple_titer or single_titer)
# 13. lawn_notes
# 14. lawn_reliability (0 = unreliable; 3 = reliable)
# 15. tested_titer
# 16. phage_reliability (0 = unreliable; 3 = reliable)
# 17. observed_infection_strength
# 18. observed_turbidity
# 19. observed_plaque_size
# 20. observed_plaques
# 21. Infection score (systematically scored infection phenotype)
immunity_data <- read.csv(IMMUNITY_DATA_FILENAME,sep=",",header=TRUE)
immunity_data$immunity_assay_id <- as.factor(immunity_data$immunity_assay_id)
immunity_data$immunity_set <- as.factor(immunity_data$immunity_set)
immunity_data$notebook <- as.factor(immunity_data$notebook)
immunity_data$page <- as.factor(immunity_data$page)
immunity_data$lawn_reliability <- as.factor(immunity_data$lawn_reliability)
immunity_data$phage_reliability <- as.factor(immunity_data$phage_reliability)
immunity_data$score <- as.factor(immunity_data$score)

# These fields contain NA's:
# observed_infection_strength
# observed_turbidity
# observed_plaque_size
# observed_plaques

# These fields contain 'unspecified', which need to be converted to NA and
# re-factored:
# prophage
# repressor_clone
# tested_titer
# phage_reliability
immunity_data[immunity_data == "Unspecified"] <- NA
immunity_data[immunity_data == "unspecified"] <- NA
immunity_data$prophage <- factor(immunity_data$prophage)
immunity_data$repressor_clone <- factor(immunity_data$repressor_clone)
immunity_data$phage_reliability <- factor(immunity_data$phage_reliability)
immunity_data$score <- factor(immunity_data$score)

# Convert titer to numeric.
immunity_data$tested_titer <-
  as.numeric(as.character(immunity_data$tested_titer))


# Create unique identifiers
immunity_data$defending_challenging <- paste(immunity_data$defending_phage,
                                             "_",
                                             immunity_data$challenging_phage,
                                             sep="")

immunity_data$defending_challenging <-
  as.factor(immunity_data$defending_challenging)

immunity_data$experiment_id <- paste(immunity_data$assay_type,
                                     "_",
                                     immunity_data$strain_type,
                                     "_",
                                     immunity_data$defending_phage,
                                     "_",
                                     immunity_data$challenging_phage,
                                     sep="")

immunity_data$experiment_id <- as.factor(immunity_data$experiment_id)


# Import Mash genomic distance data for all phages in Actino1321 database,
# including all reciprocal data and self-comparison data.
# Data structure:
# 1. phage1_phage2
# 2. modified_mash_distance (whole genome nucleotide distance, nuc_dist)
# 3. pham_pham_dissimilarity (gene content dissimilarity, gcd)
genomic_distance_data <- read.csv(GENOMIC_DISTANCE_DATA_FILENAME,
                                  sep=",",
                                  header=TRUE)


# Change column names for better readability.
names(genomic_distance_data) <- c("phage1_phage2",
                                "nuc_dist",
                                "gcd")


# Import phage metadata.
# Data structure:
# 1. phageid
# 2. host
# 3. cluster
# 4. subcluster
# 5. size
# 6. lysogen_type (extrachromosomal, integration, none)
# 7. pham_integrase (imported as factor)
# 8. pham_para (imported as int)
# 9. source (environment or lab)
# 10. parent phage
# 11. repressor_functional (yes, no, NA; For Cluster A phages, 
#     is the immunity repressor predicted to be functional?)
# 12. temperate_empirical (no, unknown, yes, NA; For Cluster A phages,
#     can a lysogen be generated?)
# 13. repressor_hth
# 14. repressor_length_full (imported as factor)
# 15. repressor_length_nterm (imported as factor)
# 16. repressor_length_cterm (imported as factor)
# 17. pham_parb (imported as int)
# 18. gene_content_clade (clade2 = the "L5 clade")
# 19. coordinate_pleft (coordinate for manual alignment)
# 20. coordinate_repressor (coordinate for manual alignment)
# 21. coordinate_genome_center (coordinate for manual alignment)
phage_metadata <- read.csv(PHAGE_METADATA_FILENAME,sep=",",header=TRUE)
phage_metadata$pham_para <- as.factor(phage_metadata$pham_para)
phage_metadata$pham_parb <- as.factor(phage_metadata$pham_parb)


# Several fields contain variants of 'unspecified', so convert these to 'NA'.
# subcluster = Unspecified
# pham_integrase = ""
# repressor_functional = "not_applicable"
# temperate_empirical = "not_applicable"

phage_metadata[phage_metadata == "Unspecified"] <- NA
phage_metadata[phage_metadata == "unspecified"] <- NA
phage_metadata[phage_metadata == ""] <- NA
phage_metadata[phage_metadata == "not_applicable"] <- NA

phage_metadata$subcluster <- factor(phage_metadata$subcluster)
phage_metadata$pham_integrase <-
  factor(phage_metadata$pham_integrase)
phage_metadata$repressor_functional <-
  factor(phage_metadata$repressor_functional)
phage_metadata$temperate_empirical <-
  factor(phage_metadata$temperate_empirical)
phage_metadata$gene_content_clade <-
  factor(phage_metadata$gene_content_clade)
phage_metadata$repressor_length_full <-
  as.numeric(as.character(phage_metadata$repressor_length_full))
phage_metadata$repressor_length_nterm <-
  as.numeric(as.character(phage_metadata$repressor_length_nterm))
phage_metadata$repressor_length_cterm <-
  as.numeric(as.character(phage_metadata$repressor_length_cterm))
phage_metadata$coordinate_pleft <-
  as.numeric(as.character(phage_metadata$coordinate_pleft))
phage_metadata$coordinate_repressor <-
  as.numeric(as.character(phage_metadata$coordinate_repressor))
phage_metadata$coordinate_genome_center <-
  as.numeric(as.character(phage_metadata$coordinate_genome_center))



# Actino1321 Immunity Repressor protein distance data, including
# reciprocal data and self-comparison data for 336 homologs. There is no
# data for escape mutants or for parent phages that are natural mutants
# (e.g. d29, misswhite, jeffabunny) with no repressor annotated.
# Data structure:
# 1. phage1_phage2
# 2. repressor_full_mafft_phyml_dist
# 3. repressor_nterm_mafft_phyml_dist
# 4. repressor_cterm_mafft_phyml_dist
# 5. repressor_full_mafft_dist_uncorrected
# 6. repressor_nterm_mafft_dist_uncorrected
# 7. repressor_cterm_mafft_dist_uncorrected
repressor_distance_data <-
  read.csv(REPRESSOR_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 Cas4-family protein mafft distance data, including
# reciprocal data and self-comparison data for 311 homologs present 
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# 1. phage1_phage2
# 2. cas4_mafft_dist_uncorrected
cas4_distance_data <-
  read.csv(CAS4_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 EndoVII protein mafft distance data, including
# reciprocal data and self-comparison data for 306 homologs present
# in Cluster A parent phages. There is not data escape mutants.
# Data structure:
# 1. phage1_phage2
# 2. endovii_mafft_dist_uncorrected
endovii_distance_data <-
  read.csv(ENDOVII_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 DNA Polymerase protein mafft distance data, including 
# reciprocal data and self-comparison data for 311 homologs present
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# 1. phage1_phage2
# 2. dnapol_mafft_dist_uncorrected
dnapol_distance_data <-
  read.csv(DNAPOL_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 Portal protein mafft distance data, including
# reciprocal data and self-comparison data for 311 homologs present
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# 1. phage1_phage2
# 2. portal_mafft_dist_uncorrected
portal_distance_data <-
  read.csv(PORTAL_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Position weight matrix distance data, including reciprocal data
# and self-comparison data for all 327 Cluster A phage genomes 
# (including escape mutants) from Actino1321 database.
# Data structure:
# 1. phage1
# 2. phage2
# 3. dist_pearson (Pairwise distance using Pearson metric)
# 4. dist_euc (Pairwise distance using Euclidean metric)
stoperator_pwm_data <- read.csv(STOPERATOR_PWM_DATA_FILENAME,
                                sep=",",
                                header=TRUE)

stoperator_pwm_data$phage1_phage2 <- paste(stoperator_pwm_data$phage1,
                                           "_",
                                           stoperator_pwm_data$phage2,
                                           sep="")
stoperator_pwm_data$phage1_phage2 <-
  as.factor(stoperator_pwm_data$phage1_phage2)

names(stoperator_pwm_data) <- c("phage1",
                                "phage2",
                                "stoperator_pwd_dist_pearson",
                                "stoperator_pwd_dist_euc",
                                "phage1_phage2")

stoperator_pwm_data <- subset(stoperator_pwm_data,
                              select=c("phage1_phage2",
                                       "stoperator_pwd_dist_pearson",
                                       "stoperator_pwd_dist_euc"))



# Match the phage metadata to the defending phage in each immunity assay.
# There are 2 defending phages in the immunity dataset that are not in
# the Actino1321 database = 'l5_gp71-flag' and 'redrock_bxb1'. Retain only
# the immunity data involving phages that are present in Actino1321.
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('defending',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")

main_immunity_data <- merge(immunity_data,
                            phage_metadata_to_match,
                            by.x="defending_phage",
                            by.y="defending_phageid")



# Match the phage metadata to the superinfecting phage in each immunity assay.
# There are several superinfecting phages that are not in the Actino1321
# database = 'd29_mutant','l5_flag1-1','l5_ha1','l5_ha2','l5_ha3','l5_ha3_hazy'.
# Retain only data for phages that are present in Actino1321.
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('challenging',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")

main_immunity_data <- merge(main_immunity_data,
                            phage_metadata_to_match,
                            by.x="challenging_phage",
                            by.y="challenging_phageid")

# At this point, all data in main_immunity_data is derived from phages 
# that are also present in the Actino1321 database.



# Match the genomic distance data. Defending phages and superinfecting phages
# not in Actino1321 will not be matched.
main_immunity_data <- merge(main_immunity_data,
                            genomic_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2")


# Match the gene distance data. Many comparisons will not be matched.
main_immunity_data <- merge(main_immunity_data,
                            repressor_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            cas4_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            endovii_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            dnapol_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            portal_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)


# Match PWM distance data.
main_immunity_data <- merge(main_immunity_data,
                            stoperator_pwm_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)


# Compute comparison fields.
main_immunity_data <- compute_comparisons(main_immunity_data)


# Some immunity assays contain data that is not as reliable as others. 
# Retain only confident data, and discard questionable data.
main_immunity_data_unreduced <- main_immunity_data

main_immunity_data <- subset(main_immunity_data,
                             main_immunity_data$lawn_reliability != 1 &
                               main_immunity_data$phage_reliability != 1)


# Histogram to assess the range of titers used, which can be used to
# compute multiplicity of infection.

par(mar=c(4,8,8,4))
hist(log(main_immunity_data$tested_titer,10),xlim=c(0,10),
     main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",
     breaks=20)
dev.copy(pdf,"tested_titers.pdf")
dev.off()


# At this point, all data in main_immunity_data is derived from phages 
# that are present in the Actino1321 database AND only confident data.
###
###
###
###
###
###
###
###
###
### 4. Average immunity data.

# For most analyses, the averaged non-duplicated immunity data is needed. 
# Averages can be computed by unique experiment_id.
main_immunity_data$experiment_id <-
  factor(main_immunity_data$experiment_id)
main_immunity_data$score <- as.numeric(as.character(main_immunity_data$score))


# Average infection scores for each unique experiment_id identifier.
immunity_average <- aggregate(main_immunity_data[,"score"],
                              list(main_immunity_data$experiment_id),mean)
names(immunity_average) <- c("experiment_id",
                             "averaged_score") 
immunity_average$experiment_id <- factor(immunity_average$experiment_id)


# Compute the range of scores for each unique assay. First compute the minimum
# score and maximum score, then compute the range.
immunity_min <- aggregate(main_immunity_data[,"score"],
                          list(main_immunity_data$experiment_id),min)
names(immunity_min) <- c("experiment_id",
                         "min_score") 
immunity_min$experiment_id <- factor(immunity_min$experiment_id)


immunity_max <- aggregate(main_immunity_data[,"score"],
                          list(main_immunity_data$experiment_id),max)
names(immunity_max) <- c("experiment_id",
                         "max_score") 
immunity_max$experiment_id <- factor(immunity_max$experiment_id)


immunity_average <- merge(immunity_average,
                          immunity_min,
                          by.x="experiment_id",
                          by.y="experiment_id")

immunity_average <- merge(immunity_average,
                          immunity_max,
                          by.x="experiment_id",
                          by.y="experiment_id")

immunity_average$range_score <- immunity_average$max_score -
  immunity_average$min_score
immunity_average$range_score <- as.factor(immunity_average$range_score)


# Create table of immunity data metadata that will be added back to the
# averaged experiment_id data. Only retain the indicated columns.
reduced_immunity_metadata1 <- subset(main_immunity_data,
                                     select=c('experiment_id',
                                              'defending_challenging',
                                              'prophage',
                                              'repressor_clone',
                                              'strain_type',
                                              'defending_phage',
                                              'challenging_phage',
                                              'assay_type'))

reduced_immunity_metadata1 <-
  reduced_immunity_metadata1[!duplicated(reduced_immunity_metadata1),]

reduced_immunity_metadata1$experiment_id <-
  factor(reduced_immunity_metadata1$experiment_id)

reduced_immunity_metadata1$defending_challenging <-
  factor(reduced_immunity_metadata1$defending_challenging)

reduced_immunity_metadata1$prophage <-
  factor(reduced_immunity_metadata1$prophage)

reduced_immunity_metadata1$repressor_clone <-
  factor(reduced_immunity_metadata1$repressor_clone)

reduced_immunity_metadata1$strain_type <-
  factor(reduced_immunity_metadata1$strain_type)

reduced_immunity_metadata1$defending_phage <-
  factor(reduced_immunity_metadata1$defending_phage)

reduced_immunity_metadata1$challenging_phage <-
  factor(reduced_immunity_metadata1$challenging_phage)

reduced_immunity_metadata1$assay_type <-
  factor(reduced_immunity_metadata1$assay_type)

immunity_average <- merge(immunity_average,
                          reduced_immunity_metadata1,
                          by.x='experiment_id',
                          by.y='experiment_id')


# Match the phage metadata.
phage_metadata_to_match <- phage_metadata

names(phage_metadata_to_match) <- paste('defending',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")

immunity_average <- merge(immunity_average,
                          phage_metadata_to_match,
                          by.x="defending_phage",
                          by.y="defending_phageid")


phage_metadata_to_match <- phage_metadata

names(phage_metadata_to_match) <- paste('challenging',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")

immunity_average <- merge(immunity_average,
                          phage_metadata_to_match,
                          by.x="challenging_phage",
                          by.y="challenging_phageid")


# Match the genomic distance data. Defending phages and challenging phages
# not in the Actino1321 database will not be matched.
immunity_average <- merge(immunity_average,
                          genomic_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2")

# Match the gene distance data. Many comparisons will not be matched.
immunity_average <- merge(immunity_average,
                          repressor_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          cas4_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          endovii_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          dnapol_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          portal_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          stoperator_pwm_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)


# Compute comparison fields.
immunity_average <- compute_comparisons(immunity_average)


# Compute frequency of each unique ave_data identifier and append to
# the averaged datasets.
main_immunity_data_freq <-
  as.data.frame(table(main_immunity_data$experiment_id))
names(main_immunity_data_freq) <- c('experiment_id','frequency')

immunity_average <- merge(immunity_average,
                          main_immunity_data_freq,
                          by.x='experiment_id',
                          by.y='experiment_id')

immunity_average$frequency <- as.factor(immunity_average$frequency) 


# Export all averaged data.
# This table was not used for the publication.
# write.table(immunity_average,
#             paste(DIR_OUTPUT,"immunity_data_averaged.csv",sep=""),
#             sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

# Export reduced data for Table S1.
output_fields <- c("strain_type",
                   "defending_phage",
                   "challenging_phage",
                   "nuc_dist",
                   "gcd",
                   "repressor_full_mafft_dist_uncorrected",
                   "stoperator_pwd_dist_euc",
                   "frequency",
                   "min_score",
                   "max_score",
                   "averaged_score")

immunity_average_reduced_for_output <- 
  subset(immunity_average,
         select = output_fields,
         immunity_average$assay_type == "multiple_titer")

output_fields_renamed <- c("strain type",
                           "defending phage",
                           "challenging phage",
                           "whole genome nucleotide distance",
                           "whole genome gene content dissimilarity",
                           "repressor distance",
                           "stoperator motif distance",
                           "replicates",
                           "minimum infection score",
                           "maximum infection score",
                           "average infection score")

names(immunity_average_reduced_for_output) <- output_fields_renamed

# Summary.
nrow(immunity_average_reduced_for_output)
sum(as.numeric(as.character(immunity_average_reduced_for_output$replicates)))

write.table(immunity_average_reduced_for_output,
            paste(DIR_OUTPUT,"Table_S1_averaged_immunity_data.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

# Summary - Data in Table S1 was used to generate heatmaps in Excel
# for Figs. 5a, 6b, 10c, S8a, S8b, S8d, S8e.

###
###
###
###
###
###
###
###
###
### 5. Compare variability between replicates.
# Map averages back to main immunity table, subtract average from original
# values, then plot histogram of differences to show how reliable the
# dataset is.
immunity_ave_to_match <- subset(immunity_average,
                                select = c("experiment_id",
                                           "averaged_score",
                                           "min_score",
                                           "max_score",
                                           "range_score",
                                           "frequency"))

names(immunity_ave_to_match) <- paste('ave_data',
                                      '_',
                                      names(immunity_ave_to_match),
                                      sep="")


# At this point main_immunity_data retains only phages that are present in
# the Actino1321 database and only confident data.
main_immunity_data <- merge(main_immunity_data,
                            immunity_ave_to_match,
                            by.x='experiment_id',
                            by.y='ave_data_experiment_id',
                            all.x=TRUE)

main_immunity_data$ave_data_averaged_score_diff <-
  as.numeric(as.character(main_immunity_data$score)) -
  as.numeric(as.character(main_immunity_data$ave_data_averaged_score))
# Note: if the diff is positive = assay's score is higher than the average.


immunity_average$frequency <-
  as.numeric(as.character(immunity_average$frequency))

immunity_average$frequency2 <-
  ifelse(immunity_average$frequency > 10,
         11,
         immunity_average$frequency)

immunity_average$frequency <- as.factor(immunity_average$frequency)
immunity_average$frequency2 <- as.factor(immunity_average$frequency2)


# QC: Assess how many replicates there are for each unique assay.
plot_bargraph2(immunity_average,
               "frequency",
               c(0,800),
               "immunity_assay_average_frequency.pdf")

plot_bargraph2(immunity_average,
               "frequency2",
               c(0,800),
               "immunity_assay_average_frequency.pdf_adjusted.pdf")


# Summary - number of unique multi-titer assays = 1,050.
immunity_ave_multi <- subset(immunity_average,
                             immunity_average$assay_type == "multiple_titer")

nrow(immunity_ave_multi)


# Summary - number of unique multi-titer assays with > 1 replicate = 67%.
immunity_ave_multi_reps <- subset(immunity_ave_multi,
                                  immunity_ave_multi$frequency != "1")
nrow(immunity_ave_multi_reps)/nrow(immunity_ave_multi)


# Summary - number of unique multi-titer assays with > 1 replicate and
# score range < 2 = 82%.
nrow(subset(immunity_ave_multi_reps,
            immunity_ave_multi_reps$range_score == "0" |
              immunity_ave_multi_reps$range_score == "1")) / 
  nrow(immunity_ave_multi_reps)

plot_bargraph2(immunity_ave_multi_reps,
               "range_score",
               c(0,400),
               "immunity_assay_replicate_range.pdf")

###
###
###
###
###
###
###
###
###
### 6. Misc. general analyses of averaged immunity data.

immunity_ave_subset1 <- subset(
  immunity_average,
  immunity_average$strain_type == 'lysogen' &
    immunity_average$assay_type == 'multiple_titer' &
    immunity_average$source_compare == 'environment' &
    immunity_average$temperate_empirical_compare == 'yes' &
    immunity_average$functional_repressor_compare == 'yes'
)

immunity_ave_subset1$defending_challenging <-
  factor(immunity_ave_subset1$defending_challenging)

immunity_ave_subset1$defending_phage <-
  factor(immunity_ave_subset1$defending_phage)

immunity_ave_subset1$challenging_phage <-
  factor(immunity_ave_subset1$challenging_phage)


# Split data by clade and by identical (same) or different (diff) phages.
immunity_ave_subset1_interclade <- 
  subset_interclade(immunity_ave_subset1,
                    "gene_content_clade_compare")

immunity_ave_subset1_intraclade2 <- 
  subset_intraclade(immunity_ave_subset1,
                    "gene_content_clade_compare",
                    "clade2")

immunity_ave_subset1_intraclade2_same <-
  subset_same(immunity_ave_subset1_intraclade2,
                   "defending_phage",
                   "challenging_phage")

immunity_ave_subset1_intraclade2_diff <-
  subset_diff(immunity_ave_subset1_intraclade2,
                     "defending_phage",
                     "challenging_phage")


# Fig. 5b sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_vs_infection_score.pdf")


# Fig. 5b sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_environment_stop_vs_infection_score.pdf")


# Fig. S5a sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "gcd",
                       "averaged_score",
                       c(0,1),
                       c(0,6),
                       "lysogen_environment_gcd_vs_infection_score.pdf")


# Fig. S5a sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "nuc_dist",
                       "averaged_score",
                       c(0,0.5),
                       c(0,6),
                       "lysogen_environment_nuc_vs_infection_score.pdf")


# Fig. S5d sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "portal_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_portal_vs_infection_score.pdf")


# Fig. S5d sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_pol_vs_infection_score.pdf")


# Fig. S5d sub-panel 3
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_endovii_vs_infection_score.pdf")


# Fig. S5d sub-panel 4
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_cas4_vs_infection_score.pdf")


# Fig. 9b sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_nterm_vs_infection_score.pdf")


# Fig. 9b sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "repressor_hth_compare",
                       "averaged_score",
                       c(0,10),
                       c(0,6),
                       "lysogen_environment_rep_hth_vs_infection_score.pdf")


# Fig. 9b sub-panel 3
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_diff,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_same,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_score",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_cterm_vs_infection_score.pdf")


# Compare differences between prophage maintenance strategy.
int_int <- subset(immunity_ave_subset1,
                  immunity_ave_subset1$lysogen_type_compare == 'integration')

extra_extra <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'extrachromosomal')

int_extra <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'different' &
    immunity_ave_subset1$defending_lysogen_type == 'integration' &
    immunity_ave_subset1$challenging_lysogen_type == 'extrachromosomal')

extra_int <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'different' &
    immunity_ave_subset1$defending_lysogen_type == 'extrachromosomal' &
    immunity_ave_subset1$challenging_lysogen_type == 'integration')

int_int$defending_pham_integrase <- factor(int_int$defending_pham_integrase)
extra_extra$defending_pham_parb <- factor(extra_extra$defending_pham_parb)


int_int_same <- subset(int_int,
                       int_int$integrase_compare != "different")
int_int_diff <- subset(int_int,
                       int_int$integrase_compare == "different")

extra_extra_same <- subset(extra_extra,
                           extra_extra$parb_compare != "different")
extra_extra_diff <- subset(extra_extra,
                           extra_extra$parb_compare == "different")


# Integrating and extrachromosomal.
int_extra_interclade <- subset_interclade(int_extra,
                                          "gene_content_clade_compare")
int_extra_intraclade2 <- subset_intraclade(int_extra,
                                           "gene_content_clade_compare",
                                           "clade2")

int_extra_intraclade2_same <- subset_same(int_extra_intraclade2,
                                                    "defending_phage",
                                                    "challenging_phage")

int_extra_intraclade2_diff <- subset_diff(int_extra_intraclade2,
                                                        "defending_phage",
                                                        "challenging_phage")

# QC
nrow(int_extra)
nrow(int_extra_interclade)
nrow(int_extra_intraclade2)
nrow(int_extra_intraclade2_same)
nrow(int_extra_intraclade2_diff)


# Extrachromosomal and integrating.
extra_int_interclade <- subset_interclade(extra_int,
                                          "gene_content_clade_compare")
extra_int_intraclade2 <- subset_intraclade(extra_int,
                                           "gene_content_clade_compare",
                                           "clade2")

extra_int_intraclade2_same <- subset_same(extra_int_intraclade2,
                                                    "defending_phage",
                                                    "challenging_phage")

extra_int_intraclade2_diff <- subset_diff(extra_int_intraclade2,
                                                        "defending_phage",
                                                        "challenging_phage")

# QC
nrow(extra_int)
nrow(extra_int_interclade)
nrow(extra_int_intraclade2)
nrow(extra_int_intraclade2_same)
nrow(extra_int_intraclade2_diff)


# Both integrating and same integrase pham.
int_int_same_interclade <- subset_interclade(int_int_same,
                                             "gene_content_clade_compare")
int_int_same_intraclade2 <- subset_intraclade(int_int_same,
                                              "gene_content_clade_compare",
                                              "clade2")

int_int_same_intraclade2_same <-
  subset_same(int_int_same_intraclade2,
                   "defending_phage",
                   "challenging_phage")

int_int_same_intraclade2_diff <- 
  subset_diff(
    int_int_same_intraclade2,
    "defending_phage",
    "challenging_phage")

# QC
nrow(int_int_same)
nrow(int_int_same_interclade)
nrow(int_int_same_intraclade2)
nrow(int_int_same_intraclade2_same)
nrow(int_int_same_intraclade2_diff)


# Both integrating but different integrase phams.
int_int_diff_interclade <- subset_interclade(int_int_diff,
                                             "gene_content_clade_compare")
int_int_diff_intraclade2 <- subset_intraclade(int_int_diff,
                                              "gene_content_clade_compare",
                                              "clade2")

int_int_diff_intraclade2_same <-
  subset_same(int_int_diff_intraclade2,
                   "defending_phage",
                   "challenging_phage")

int_int_diff_intraclade2_diff <-
  subset_diff(int_int_diff_intraclade2,
                     "defending_phage",
                     "challenging_phage")

# QC
nrow(int_int_diff)
nrow(int_int_diff_interclade)
nrow(int_int_diff_intraclade2)
nrow(int_int_diff_intraclade2_same)
nrow(int_int_diff_intraclade2_diff)


# Both extrachromosomal and same ParB pham.
extra_extra_same_interclade <- subset_interclade(extra_extra_same,
                                                 "gene_content_clade_compare")
extra_extra_same_intraclade2 <- subset_intraclade(extra_extra_same,
                                                  "gene_content_clade_compare",
                                                  "clade2")

extra_extra_same_intraclade2_same <-
  subset_same(extra_extra_same_intraclade2,
                   "defending_phage",
                   "challenging_phage")

extra_extra_same_intraclade2_diff <-
  subset_diff(extra_extra_same_intraclade2,
                     "defending_phage",
                     "challenging_phage")

# QC
nrow(extra_extra_same)
nrow(extra_extra_same_interclade)
nrow(extra_extra_same_intraclade2)
nrow(extra_extra_same_intraclade2_same)
nrow(extra_extra_same_intraclade2_diff)


# Both extrachromosomal but different ParB phams.
extra_extra_diff_interclade <- subset_interclade(extra_extra_diff,
                                                 "gene_content_clade_compare")
extra_extra_diff_intraclade2 <- subset_intraclade(extra_extra_diff,
                                                  "gene_content_clade_compare",
                                                  "clade2")

extra_extra_diff_intraclade2_same <-
  subset_same(extra_extra_diff_intraclade2,
                   "defending_phage",
                   "challenging_phage")

extra_extra_diff_intraclade2_diff <-
  subset_diff(extra_extra_diff_intraclade2,
                     "defending_phage",
                     "challenging_phage")

# QC
nrow(extra_extra_diff)
nrow(extra_extra_diff_interclade)
nrow(extra_extra_diff_intraclade2)
nrow(extra_extra_diff_intraclade2_same)
nrow(extra_extra_diff_intraclade2_diff)


# Fig. S5b sub-panel 1
plot_tricolor_scatter2(int_int_same_intraclade2_diff,
                       int_int_same_interclade,
                       int_int_same_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_same_stop_vs_infection_score.pdf")


# Fig. S5b sub-panel 2
plot_tricolor_scatter2(int_int_diff_intraclade2_diff,
                       int_int_diff_interclade,
                       int_int_diff_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_diff_stop_vs_infection_score.pdf")


# Fig. S5b sub-panel 3
plot_tricolor_scatter2(int_extra_intraclade2_diff,
                       int_extra_interclade,
                       int_extra_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_parb_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 1
plot_tricolor_scatter2(extra_extra_same_intraclade2_diff,
                       extra_extra_same_interclade,
                       extra_extra_same_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_same_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 2
plot_tricolor_scatter2(extra_extra_diff_intraclade2_diff,
                       extra_extra_diff_interclade,
                       extra_extra_diff_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_diff_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 3
plot_tricolor_scatter2(extra_int_intraclade2_diff,
                       extra_int_interclade,
                       extra_int_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_int_stop_vs_infection_score.pdf")


# Binned Data Stats.
# Unlike analyses above, this subsetted data includes only L5-clade phages,
# and regardless of whether they are environmental or have a
# functional repressor.

immunity_ave_subset2 <- subset(
  immunity_average,
  immunity_average$strain_type == 'lysogen' &
    immunity_average$assay_type == 'multiple_titer')

immunity_ave_subset2$defending_lysogen_type <-
  factor(immunity_ave_subset2$defending_lysogen_type)

immunity_ave_subset2$challenging_lysogen_type <-
  factor(immunity_ave_subset2$challenging_lysogen_type)

immunity_ave_subset2$challenging_phage <-
  factor(immunity_ave_subset2$challenging_phage)

immunity_ave_subset2$defending_phage <-
  factor(immunity_ave_subset2$defending_phage)

clade2_env <- subset(
  immunity_ave_subset2,
  immunity_ave_subset2$gene_content_clade_compare == "clade2" &
    immunity_ave_subset2$source_compare == "environment"
)

clade2_env <- subset(clade2_env,
                     as.character(clade2_env$defending_phage) != 
                       as.character(clade2_env$challenging_phage))

clade2_env$defending_phage <- factor(clade2_env$defending_phage)
clade2_env$challenging_phage <- factor(clade2_env$challenging_phage)

clade2_env_bin0 <- subset(clade2_env,
                          clade2_env$averaged_score <= 0.5)
clade2_env_bin1 <- subset(clade2_env,
                          clade2_env$averaged_score > 0.5 &
                            clade2_env$averaged_score <= 1.5)
clade2_env_bin2 <- subset(clade2_env,
                          clade2_env$averaged_score > 1.5 &
                            clade2_env$averaged_score <= 2.5)
clade2_env_bin3 <- subset(clade2_env,
                          clade2_env$averaged_score > 2.5 &
                            clade2_env$averaged_score <= 3.5)
clade2_env_bin4 <- subset(clade2_env,
                          clade2_env$averaged_score > 3.5 &
                            clade2_env$averaged_score <= 4.5)
clade2_env_bin5 <- subset(clade2_env,
                          clade2_env$averaged_score > 4.5 &
                            clade2_env$averaged_score <= 5.5)
clade2_env_bin6 <- subset(clade2_env,
                          clade2_env$averaged_score > 5.5)

clade2_bin0_freq <- compute_binned_freq1(clade2_env_bin0,"bin0")
clade2_bin1_freq <- compute_binned_freq1(clade2_env_bin1,"bin1")
clade2_bin2_freq <- compute_binned_freq1(clade2_env_bin2,"bin2")
clade2_bin3_freq <- compute_binned_freq1(clade2_env_bin3,"bin3")
clade2_bin4_freq <- compute_binned_freq1(clade2_env_bin4,"bin4")
clade2_bin5_freq <- compute_binned_freq1(clade2_env_bin5,"bin5")
clade2_bin6_freq <- compute_binned_freq1(clade2_env_bin6,"bin6")

clade2_binned_frequency <- rbind(clade2_bin0_freq,
                                 clade2_bin1_freq,
                                 clade2_bin2_freq,
                                 clade2_bin3_freq,
                                 clade2_bin4_freq,
                                 clade2_bin5_freq,
                                 clade2_bin6_freq)

clade2_binned_frequency$bin <- factor(clade2_binned_frequency$bin,
                                      c("bin6",
                                        "bin5",
                                        "bin4",
                                        "bin3",
                                        "bin2",
                                        "bin1",
                                        "bin0"))


# Fig. S4a sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "defending_percent",
               "bin",
               c(0,1),
               "infection_score_percent_defending.pdf")


# Fig. S4a sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "challenging_percent",
               "bin",
               c(0,1),
               "infection_score_percent_challenging.pdf")


# Fig. S4a sub-panel 3
plot_bargraph1(clade2_binned_frequency,
               "total_assays_percent",
               "bin",
               c(0,0.3),
               "infection_score_percent_total.pdf")


# Fig. S4b sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "inter_subcluster_percent",
               "bin",
               c(0,0.5),
               "infection_score_percent_intersubcluster.pdf")


# Fig. S4b sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "intra_subcluster_percent",
               "bin",
               c(0,0.5),
               "infection_score_percent_intrasubcluster.pdf")
###
###
###
###
###
###
###
###
###
### 7. Compare reciprocal infection assays.
# Using averaged data, compute differences in reciprocal immunity data.
# Vectored_data = data that is specific to the immunity assay
# and depends on which phage is infecting and defending:
vectored_column_names <- c("experiment_id",
                           "defending_challenging",
                           "defending_phage",
                           "challenging_phage",
                           "averaged_score",
                           "assay_type",
                           "frequency"
)

reciprocal_data <- subset(immunity_average,
                          immunity_average$strain_type == 'lysogen' & 
                            immunity_average$assay_type == 'multiple_titer')

vector2_data <- subset(reciprocal_data,select = vectored_column_names)
names(vector2_data) <- paste('vector2_',names(vector2_data),sep="")

for (column_name in vectored_column_names){
  
  new_column_name <- paste('vector1_',column_name,sep="")
  
  print(new_column_name)
  
  names(reciprocal_data)[names(reciprocal_data) == column_name] <-
    paste('vector1_',column_name,sep="")
}

reciprocal_data$vector2_defending_challenging <- 
  paste(reciprocal_data$vector1_challenging_phage,
        "_",
        reciprocal_data$vector1_defending_phage,
        sep="")

reciprocal_data <- merge(reciprocal_data,
                         vector2_data,
                         by.x='vector2_defending_challenging',
                         by.y='vector2_defending_challenging')

# Compute difference in infection profiles.
reciprocal_data$averaged_score_diff <-
  abs(reciprocal_data$vector1_averaged_score - 
        reciprocal_data$vector2_averaged_score)

# Now that all data has been matched, the data is duplicated. (e.g. vector1
# will contain l5_alma and alma_l5, so both vector2's will be present as well).
# To solve this, remove all vector1 rows in which the defending_phage and
# challenging_phage are not alphabetically ordered. This does NOT retain
# self-comparisons (e.g. alma_alma).
reciprocal_data$vector1_alpha_ordered <- 
  as.character(reciprocal_data$vector1_defending_phage) < 
  as.character(reciprocal_data$vector1_challenging_phage)

# Now retain only the unique pairwise comparisons.
reciprocal_data_alpha_ordered <- 
  subset(reciprocal_data,
         reciprocal_data$vector1_alpha_ordered == TRUE)

# Summary
nrow(reciprocal_data_alpha_ordered)

# Output the reciprocal dataset.
# This table was not used for the publication.
# write.table(reciprocal_data_alpha_ordered,
#             paste(DIR_OUTPUT,"reciprocal_immunity_data.csv",sep=""),
#             sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


# Bin reciprocal data for histogram.

# Reduce the reciprocal dataset to only environmental phages.
reciprocal_unique_envY <- 
  subset(reciprocal_data,
         reciprocal_data$vector1_alpha_ordered == TRUE &
           reciprocal_data$source_compare == 'environment')


# Bin data.
reciprocal_bin0 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff <= 0.5)
reciprocal_bin1 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 0.5 &
                            reciprocal_unique_envY$averaged_score_diff <= 1.5)
reciprocal_bin2 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 1.5 &
                            reciprocal_unique_envY$averaged_score_diff <= 2.5)
reciprocal_bin3 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 2.5 &
                            reciprocal_unique_envY$averaged_score_diff <= 3.5)
reciprocal_bin4 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 3.5 &
                            reciprocal_unique_envY$averaged_score_diff <= 4.5)
reciprocal_bin5 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 4.5 &
                            reciprocal_unique_envY$averaged_score_diff <= 5.5)
reciprocal_bin6 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_score_diff > 5.5)

# Compute frequency.
reciprocal_bin0_freq <- compute_binned_freq2(reciprocal_bin0,
                                             "bin0",
                                             reciprocal_unique_envY)
reciprocal_bin1_freq <- compute_binned_freq2(reciprocal_bin1,
                                             "bin1",
                                             reciprocal_unique_envY)
reciprocal_bin2_freq <- compute_binned_freq2(reciprocal_bin2,
                                             "bin2",
                                             reciprocal_unique_envY)
reciprocal_bin3_freq <- compute_binned_freq2(reciprocal_bin3,
                                             "bin3",
                                             reciprocal_unique_envY)
reciprocal_bin4_freq <- compute_binned_freq2(reciprocal_bin4,
                                             "bin4",
                                             reciprocal_unique_envY)
reciprocal_bin5_freq <- compute_binned_freq2(reciprocal_bin5,
                                             "bin5",
                                             reciprocal_unique_envY)
reciprocal_bin6_freq <- compute_binned_freq2(reciprocal_bin6,
                                             "bin6",
                                             reciprocal_unique_envY)

reciprocal_binned_freq <- rbind(reciprocal_bin0_freq,
                                reciprocal_bin1_freq,
                                reciprocal_bin2_freq,
                                reciprocal_bin3_freq,
                                reciprocal_bin4_freq,
                                reciprocal_bin5_freq,
                                reciprocal_bin6_freq)

reciprocal_binned_freq$bin <- factor(reciprocal_binned_freq$bin,
                                     c("bin6",
                                       "bin5",
                                       "bin4",
                                       "bin3",
                                       "bin2",
                                       "bin1",
                                       "bin0"))

# QC: Should equal 0.
sum(reciprocal_binned_freq$freq) - nrow(reciprocal_unique_envY)


# Fig. S4c
plot_bargraph1(reciprocal_binned_freq,
               "freq_percent",
               "bin",
               c(0,0.5),
               "reciprocal_infection_scores_percent_total.pdf")


# Subsetted reciprocal_unique_envY contains data for confident assays involving
# lysogen, multi-titer, and environmental phages.Since the data is reciprocal,
# it by definition reflects empirically determined temperate phages with
# functional repressors.
reciprocal_interclade <- 
  subset_interclade(reciprocal_unique_envY,
                    "gene_content_clade_compare")

reciprocal_intraclade2 <- 
  subset_intraclade(reciprocal_unique_envY,
                    "gene_content_clade_compare",
                    "clade2")

reciprocal_intraclade2_same <- 
  subset_same(reciprocal_intraclade2,
                   "vector1_defending_phage",
                   "vector1_challenging_phage")

reciprocal_intraclade2_diff <- 
  subset_diff(reciprocal_intraclade2,
                     "vector1_defending_phage",
                     "vector1_challenging_phage")


# Fig. 5b sub-panel 3
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_vs_infection_score.pdf")


# Fig. 5b sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "stoperator_pwd_dist_euc",
                       "averaged_score_diff",
                       c(0,5),
                       c(0,6),
                       "reciprocal_lys_env_stop_vs_infection_score.pdf")


# Fig. S5a sub-panel 3
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "gcd",
                       "averaged_score_diff",
                       c(0,1),
                       c(0,6),
                       "reciprocal_lys_env_gcd_vs_infection_score.pdf")


# Fig. S5a sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "nuc_dist",
                       "averaged_score_diff",
                       c(0,0.5),
                       c(0,6),
                       "reciprocal_lys_env_nuc_vs_infection_score.pdf")


# Fig. S5d sub-panel 5
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "portal_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_portal_vs_infection_score.pdf")


# Fig. S5d sub-panel 6
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_pol_vs_infection_score.pdf")


# Fig. S5d sub-panel 7
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_endovii_vs_infection_score.pdf")


# Fig. S5d sub-panel 8
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_cas4_vs_infection_score.pdf")


# Fig. 9b sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_nterm_vs_infection_score.pdf")


# Fig. 9b sub-panel 5
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "repressor_hth_compare",
                       "averaged_score_diff",
                       c(0,10),
                       c(0,6),
                       "reciprocal_lys_env_rep_hth_vs_infection_score.pdf")


# Fig. 9b sub-panel 6
plot_tricolor_scatter2(reciprocal_intraclade2_diff,
                       reciprocal_interclade,
                       reciprocal_intraclade2_same,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_score_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_cterm_vs_infection_score.pdf")

###
###
###
###
###
###
###
###
###
### 8. Compare lysogen and cloned-repressor strains.

# Summary - total lysogen-CRS reciprocal assays performed.
lys_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$strain_type == 'lysogen')

clone_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$strain_type == 'repressor_clone')

lys_clone_matched <- match_lys_clone_data(lys_to_match,
                                          clone_to_match)

nrow(lys_clone_matched)


# Environment, functional repressor, and empirically temperate data from
# averaged multiple-titer lysogen and clone data.
env_lys_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$source_compare == 'environment' &
           immunity_average$functional_repressor_compare == 'yes' &
           immunity_average$temperate_empirical_compare == 'yes' &
           immunity_average$strain_type == 'lysogen')

env_clone_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$source_compare == 'environment' &
           immunity_average$functional_repressor_compare == 'yes' &
           immunity_average$temperate_empirical_compare == 'yes' &
           immunity_average$strain_type == 'repressor_clone')

env_lys_clone_matched <- match_lys_clone_data(env_lys_to_match,
                                              env_clone_to_match)

env_lys_clone_interclade <- 
  subset_interclade(env_lys_clone_matched,
                    "lys_gene_content_clade_compare")

env_lys_clone_intraclade2 <- 
  subset_intraclade(env_lys_clone_matched,
                    "lys_gene_content_clade_compare",
                    "clade2")

env_lys_clone_intraclade2_same <- 
  subset_diff(env_lys_clone_intraclade2,
                     "lys_defending_phage",
                     "lys_challenging_phage")

env_lys_clone_intraclade2_diff <- 
  subset_same(env_lys_clone_intraclade2,
              "lys_defending_phage",
              "lys_challenging_phage")


# Escape mutant data from averaged multiple-titer lysogen and clone data.
escape_lys_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$source_compare != 'environment' &
           immunity_average$strain_type == 'lysogen')

escape_clone_to_match <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$source_compare != 'environment' &
           immunity_average$strain_type == 'repressor_clone')

escape_lys_clone_matched <- match_lys_clone_data(escape_lys_to_match,
                                                 escape_clone_to_match)


# Summary
nrow(escape_lys_clone_matched) + nrow(env_lys_clone_matched)

escape_lys_clone_interclade <- 
  subset_interclade(escape_lys_clone_matched,
                    "lys_gene_content_clade_compare")

escape_lys_clone_intraclade2 <- 
  subset_intraclade(escape_lys_clone_matched,
                    "lys_gene_content_clade_compare",
                    "clade2")

escape_lys_clone_intraclade2_same <- 
  subset_same(escape_lys_clone_intraclade2,
                   "lys_defending_phage",
                   "lys_challenging_phage")

escape_lys_clone_intraclade2_diff <- 
  subset_diff(escape_lys_clone_intraclade2,
                     "lys_defending_phage",
                     "lys_challenging_phage")

# QC
nrow(env_lys_clone_matched)
nrow(escape_lys_clone_matched)

nrow(env_lys_clone_intraclade2_same)
nrow(env_lys_clone_intraclade2_diff)
nrow(env_lys_clone_interclade)

nrow(env_lys_clone_intraclade2_same) +
  nrow(env_lys_clone_intraclade2_diff) + 
  nrow(env_lys_clone_interclade)

nrow(escape_lys_clone_intraclade2_diff)
nrow(escape_lys_clone_interclade)
nrow(escape_lys_clone_intraclade2_same)

nrow(escape_lys_clone_intraclade2_diff) + 
  nrow(escape_lys_clone_interclade) + 
  nrow(escape_lys_clone_intraclade2_same)


# Fig. 6c
plot_tricolor_scatter1(env_lys_clone_intraclade2_same,
                       env_lys_clone_interclade,
                       env_lys_clone_intraclade2_diff,
                       "lys_averaged_score","clone_averaged_score",
                       c(0,6),
                       c(0,6),
                       "lys_vs_crs_env_infection_score.pdf")

compute_correlation(env_lys_clone_matched,
                    "lys_averaged_score",
                    "clone_averaged_score")


# Fig. 7d
plot_tricolor_scatter1(escape_lys_clone_intraclade2_diff,
                       escape_lys_clone_interclade,
                       escape_lys_clone_intraclade2_same,
                       "lys_averaged_score","clone_averaged_score",
                       c(0,6),
                       c(0,6),
                       "lys_vs_crs_escape_infection_score.pdf")

compute_correlation(escape_lys_clone_matched,
                    "lys_averaged_score",
                    "clone_averaged_score")


# Fig. 6d
plot_tricolor_scatter2(env_lys_clone_intraclade2_same,
                       env_lys_clone_interclade,
                       env_lys_clone_intraclade2_diff,
                       "lys_stoperator_pwd_dist_euc",
                       "score_diff",
                       c(0,5),
                       c(-4,4),
                       "lys_vs_crs_env_inf_score_stop_vs_score_diff.pdf",
                       "yes")


# Fig. 7e
plot_tricolor_scatter2(escape_lys_clone_intraclade2_diff,
                       escape_lys_clone_interclade,
                       escape_lys_clone_intraclade2_same,
                       "lys_stoperator_pwd_dist_euc",
                       "score_diff",
                       c(0,5),
                       c(-4,4),
                       "lys_vs_crs_escape_inf_score_stop_vs_score_diff.pdf",
                       "yes")

###
###
###
###
###
###
###
###
###
### 9. Compare L5 and L5-derivative phages.
# All averaged multiple-titer lysogen and clone data.
l5_columns <- c("defending_challenging",
                "defending_phage",
                "challenging_phage",
                "averaged_score",
                "nuc_dist",
                "gcd",
                "repressor_cterm_mafft_dist_uncorrected")

immunity_ave_subset3 <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$strain_type == 'lysogen')

# Compare L5-related superinfection patterns.
chal_l5 <- subset(immunity_ave_subset3,
                  immunity_ave_subset3$challenging_phage == 'l5')

chal_phitm41 <- subset(immunity_ave_subset3,
                       immunity_ave_subset3$challenging_phage == 'phitm41',
                       select = c(l5_columns))
chal_phitm1 <- subset(immunity_ave_subset3,
                      immunity_ave_subset3$challenging_phage == 'phitm1',
                      select = c(l5_columns))
chal_phitm4 <- subset(immunity_ave_subset3,
                      immunity_ave_subset3$challenging_phage == 'phitm4',
                      select = c(l5_columns))
chal_phitm6 <- subset(immunity_ave_subset3,
                      immunity_ave_subset3$challenging_phage == 'phitm6',
                      select = c(l5_columns))

names(chal_l5) <- paste('l5_',names(chal_l5),sep="")
names(chal_phitm1) <- paste('phitm1_',names(chal_phitm1),sep="")
names(chal_phitm41) <- paste('phitm41_',names(chal_phitm41),sep="")
names(chal_phitm4) <- paste('phitm4_',names(chal_phitm4),sep="")
names(chal_phitm6) <- paste('phitm6_',names(chal_phitm6),sep="")


chal_l5_assays <- merge(chal_l5,
                        chal_phitm41,
                        by.x='l5_defending_phage',
                        by.y='phitm41_defending_phage',
                        all.x=TRUE,
                        all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,
                        chal_phitm1,
                        by.x='l5_defending_phage',
                        by.y='phitm1_defending_phage',
                        all.x=TRUE,
                        all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,
                        chal_phitm4,
                        by.x='l5_defending_phage',
                        by.y='phitm4_defending_phage',
                        all.x=TRUE,
                        all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,
                        chal_phitm6,
                        by.x='l5_defending_phage',
                        by.y='phitm6_defending_phage',
                        all.x=TRUE,
                        all.y=TRUE)

chal_l5_assays$l5_phitm41_score_diff <-
  chal_l5_assays$phitm41_averaged_score - chal_l5_assays$l5_averaged_score

chal_l5_assays$l5_phitm1_score_diff <-
  chal_l5_assays$phitm1_averaged_score - chal_l5_assays$l5_averaged_score

chal_l5_assays$l5_phitm4_score_diff <-
  chal_l5_assays$phitm4_averaged_score - chal_l5_assays$l5_averaged_score

chal_l5_assays$l5_phitm6_score_diff <-
  chal_l5_assays$phitm6_averaged_score - chal_l5_assays$l5_averaged_score

chal_l5_assays$phitm1_phitm6_score_diff <-
  chal_l5_assays$phitm6_averaged_score - chal_l5_assays$phitm1_averaged_score

chal_l5_assays$phitm1_phitm4_score_diff <-
  chal_l5_assays$phitm4_averaged_score - chal_l5_assays$phitm1_averaged_score


chal_l5_assays_interclade <-
  subset_interclade(chal_l5_assays,
                    "l5_gene_content_clade_compare")

chal_l5_assays_intraclade2 <- 
  subset_intraclade(chal_l5_assays,
                    "l5_gene_content_clade_compare",
                    "clade2")

# Used as a dummy table for plotting. "empty" is not a valid clade comparison.
chal_l5_assays_intraclade2_empty <- 
  subset_intraclade(chal_l5_assays,
                    "l5_gene_content_clade_compare",
                    "empty")

chal_l5_assays_intraclade2_phitm1_same <- 
  subset_same(chal_l5_assays_intraclade2,
                   "phitm1_challenging_phage",
                   "l5_defending_phage")

chal_l5_assays_intraclade2_phitm1_diff <- 
  subset_diff(chal_l5_assays_intraclade2,
                     "phitm1_challenging_phage",
                     "l5_defending_phage")

# QC
nrow(chal_l5_assays_intraclade2)
nrow(subset(chal_l5_assays_intraclade2,
            is.na(chal_l5_assays_intraclade2$phitm1_challenging_phage)))
nrow(chal_l5_assays_intraclade2_phitm1_same)
nrow(chal_l5_assays_intraclade2_phitm1_diff)
nrow(chal_l5_assays_intraclade2_empty)


# Fig. 10d sub-panel 1
plot_tricolor_scatter2(chal_l5_assays_intraclade2_phitm1_diff,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_phitm1_same,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm1_averaged_score",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_phitm1.pdf")


# Fig. 10d sub-panel 2
plot_tricolor_scatter2(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm4_averaged_score",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_phitm4.pdf")


# Fig. S8g sub-panel 1
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_score",
                       "phitm41_averaged_score",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm41.pdf")


# Fig. S8g sub-panel 2
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_score",
                       "phitm6_averaged_score",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm6.pdf")


# Fig. S8g sub-panel 3
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_score",
                       "phitm1_averaged_score",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm1.pdf")


# Compare superinfection patterns against L5-related prophages.
def_l5 <- subset(immunity_ave_subset3,
                 immunity_ave_subset3$defending_phage == 'l5')

def_phitm41 <- subset(immunity_ave_subset3,
                      immunity_ave_subset3$defending_phage == 'phitm41',
                      select = c(l5_columns))
def_phitm1 <- subset(immunity_ave_subset3,
                     immunity_ave_subset3$defending_phage == 'phitm1',
                     select = c(l5_columns))
def_phitm6 <- subset(immunity_ave_subset3,
                     immunity_ave_subset3$defending_phage == 'phitm6',
                     select = c(l5_columns))

names(def_l5) <- paste('l5_',names(def_l5),sep="")
names(def_phitm41) <- paste('phitm41_',names(def_phitm41),sep="")
names(def_phitm1) <- paste('phitm1_',names(def_phitm1),sep="")
names(def_phitm6) <- paste('phitm6_',names(def_phitm6),sep="")


def_l5_assays <- merge(def_l5,
                       def_phitm41,
                       by.x='l5_challenging_phage',
                       by.y='phitm41_challenging_phage',
                       all.x=TRUE,
                       all.y=TRUE)
def_l5_assays <- merge(def_l5_assays,
                       def_phitm1,
                       by.x='l5_challenging_phage',
                       by.y='phitm1_challenging_phage',
                       all.x=TRUE,
                       all.y=TRUE)
def_l5_assays <- merge(def_l5_assays,
                       def_phitm6,
                       by.x='l5_challenging_phage',
                       by.y='phitm6_challenging_phage',
                       all.x=TRUE,
                       all.y=TRUE)


def_l5_assays$l5_phitm41_score_diff <-
  def_l5_assays$phitm41_averaged_score - def_l5_assays$l5_averaged_score

def_l5_assays$l5_phitm1_score_diff <-
  def_l5_assays$phitm1_averaged_score - def_l5_assays$l5_averaged_score

def_l5_assays$l5_phitm6_score_diff <-
  def_l5_assays$phitm6_averaged_score - def_l5_assays$l5_averaged_score

def_l5_assays$phitm1_phitm6_score_diff <-
  def_l5_assays$phitm6_averaged_score - def_l5_assays$phitm1_averaged_score

def_l5_assays_interclade <- 
  subset_interclade(def_l5_assays,
                    "l5_gene_content_clade_compare")

def_l5_assays_intraclade2 <- 
  subset_intraclade(def_l5_assays,
                    "l5_gene_content_clade_compare",
                    "clade2")

# Used as a dummy table for plotting.
# phiTM4 is lytic, so unable to form a lysogen.
def_l5_assays_clade2_empty <- 
  subset(def_l5_assays,
         def_l5_assays$l5_defending_phage == "phiTM4")


# QC
nrow(def_l5_assays)
nrow(def_l5_assays_intraclade2)
nrow(def_l5_assays_interclade)
nrow(def_l5_assays_clade2_empty)

summary(def_l5_assays$l5_challenging_phage)
summary(def_l5_assays_intraclade2$l5_challenging_phage)
summary(def_l5_assays_interclade$l5_challenging_phage)


# Fig. S8h sub-panel 1
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_score",
                       "phitm41_averaged_score",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm41.pdf")


# Fig. S8h sub-panel 2
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_score",
                       "phitm6_averaged_score",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm6.pdf")


# Fig. S8h sub-panel 3
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_score",
                       "phitm1_averaged_score",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm1.pdf")

###
###
###
###
###
###
###
###
###
### 10. Compare escape mutants and parent phages.
# Pair superinfecting parent and superinfecting mutant data. Split averaged
# data into parent and mutant phage tables. Averaged, non-redundant data
# should only be for multi-titer assays, no low-confidence data, and only
# from lysogen data. There should be no repressor clones or single-titer assays.
mutant_data_lys <- subset(immunity_average,
                      immunity_average$strain_type == 'lysogen' & 
                        immunity_average$assay_type == 'multiple_titer')

# Rename all fields to indicate the data is generated using the mutant
# challenger. Note: important to remember that "mutant" prefix refers to the
# entire immunity assay data, not to any particular column.
names(mutant_data_lys) <- paste('mutant_',names(mutant_data_lys),sep="")


# Drop all data from mutant table not involving a mutant challenging phage. 
# It is possible that a lab-derived phage is not an escape mutant, so this
# a broad list of mutants.
mutant_data_lys <- subset(mutant_data_lys,
                      mutant_data_lys$mutant_challenging_source == 'lab')


# Create parent data to match.
parent_data_lys <- subset(immunity_average,
                      immunity_average$strain_type == 'lysogen' & 
                        immunity_average$assay_type == 'multiple_titer')

names(parent_data_lys) <- paste('parent_',names(parent_data_lys),sep="")


# Now match parent challenging data to mutant challenging data.
mutant_data_lys$parent_defending_challenging <- 
  paste(mutant_data_lys$mutant_defending_phage,
        "_",
        mutant_data_lys$mutant_challenging_parent
        ,sep="")

mutant_parent_lys <- merge(mutant_data_lys,
                         parent_data_lys,
                         by.x="parent_defending_challenging",
                         by.y="parent_defending_challenging")


# Compute difference in infection profiles. The direction of change in
# infection between the mutant and parent is informative, so don't use
# absolute value.
mutant_parent_lys$averaged_score_diff <-
  mutant_parent_lys$mutant_averaged_score - 
  mutant_parent_lys$parent_averaged_score


# Export all mutant data.
# This table was not used for the publication.
# write.table(mutant_parent_lys,
#             paste(DIR_OUTPUT,"mutant_parent_lysogen_data.csv",sep=""),
#             sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


# QC
plot_hist1(mutant_parent_lys,
           "averaged_score_diff",
           10,
           c(-1,4),
           c(0,80),
           "mutant_averaged_score_diff.pdf")


# Reduce to only escape mutants. Not all lab-derived phages are escape mutants.
# So at this step, explicitly select escape mutants,

# Mutant-parent phage sets
# 1
# pioneer
# phitm35 (escape mutant)

# 2
# eagleeye
# phitm36 (escape mutant)

# 3
# phitm44 (isolated mutant)
# phitm38 (escape mutant)

# 4
# et2brutus
# phitm39 (escape mutant)	
# phitm40 (escape mutant)

# 5
# l5
# phitm41 (escape mutant)

# 6
# trixie
# phitm42 (escape mutant)

# 8 (clone EM)
# davinci
# phitm46 (escape mutant)

# 9 (clone EM)
# gladiator
# phitm47 (escape mutant)

escape_parent_lys <- 
  subset(mutant_parent_lys,
         mutant_parent_lys$mutant_challenging_phage == "phitm35" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm36" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm38" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm39" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm40" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm41" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm42" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm46" | 
           mutant_parent_lys$mutant_challenging_phage == "phitm47")

escape_interclade <- subset_interclade(escape_parent_lys,
                                       "parent_gene_content_clade_compare")

escape_intraclade2 <- subset_intraclade(escape_parent_lys,
                                        "parent_gene_content_clade_compare",
                                        "clade2")

# Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
escape_intraclade2_empty <- 
  subset_intraclade(escape_parent_lys,
                    "parent_gene_content_clade_compare",
                    "empty")

escape_intraclade2_parent_same <- subset_same(escape_intraclade2,
                                              "parent_challenging_phage",
                                                "parent_defending_phage")

escape_intraclade2_parent_diff <- subset_diff(escape_intraclade2,
                                                "parent_challenging_phage",
                                                  "parent_defending_phage")

escape_intraclade2_mutant_same <- subset_same(escape_intraclade2,
                                              "mutant_challenging_phage",
                                                "mutant_defending_phage")

escape_intraclade2_mutant_diff <- subset_diff(escape_intraclade2,
                                                "mutant_challenging_phage",
                                                  "mutant_defending_phage")

# QC
nrow(escape_intraclade2)
nrow(escape_intraclade2_empty)
nrow(escape_intraclade2_parent_same)
nrow(escape_intraclade2_parent_diff)
nrow(escape_intraclade2_mutant_same)
nrow(escape_intraclade2_mutant_diff)
summary(escape_intraclade2_parent_same$parent_defending_phage)
summary(escape_intraclade2_mutant_same$mutant_defending_phage)


# Fig. 7c sub-panel 1
plot_tricolor_scatter1(escape_intraclade2,
                       escape_interclade,
                       escape_intraclade2_empty,
                       "parent_averaged_score",
                       "mutant_averaged_score",
                       c(0,6),
                       c(0,6),
                       "parent_vs_escape_infection_scores_lysogen.pdf")


# Fig. 7f sub-panel 1
plot_tricolor_scatter2(escape_intraclade2_parent_diff,
                       escape_interclade,
                       escape_intraclade2_parent_same,
                       "parent_stoperator_pwd_dist_euc",
                       "parent_averaged_score",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_parent.pdf")


# Fig. 7f sub-panel 2
plot_tricolor_scatter2(escape_intraclade2_mutant_diff,
                       escape_interclade,
                       escape_intraclade2_mutant_same,
                       "mutant_stoperator_pwd_dist_euc",
                       "mutant_averaged_score",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_escape.pdf")


# Compare escape mutants to parents on repressor clone strains. Use same
# pipeline as for lysogen strains.
mutant_data_clone <- subset(immunity_average,
                      immunity_average$strain_type == 'repressor_clone' & 
                        immunity_average$assay_type == 'multiple_titer')

names(mutant_data_clone) <- paste('mutant_',
                                names(mutant_data_clone),
                                sep="")

mutant_data_clone <- subset(mutant_data_clone,
                          mutant_data_clone$mutant_challenging_source == 'lab')

parent_data_clone <- subset(immunity_average,
                      immunity_average$strain_type == 'repressor_clone' & 
                        immunity_average$assay_type == 'multiple_titer')

names(parent_data_clone) <- paste('parent_',names(parent_data_clone),sep="")

mutant_data_clone$parent_defending_challenging <- 
  paste(mutant_data_clone$mutant_defending_phage,
        "_",
        mutant_data_clone$mutant_challenging_parent,
        sep="")

mutant_parent_clone <- merge(mutant_data_clone,
                             parent_data_clone,
                             by.x="parent_defending_challenging",
                             by.y="parent_defending_challenging")

mutant_parent_clone$averaged_score_diff <-
  mutant_parent_clone$mutant_averaged_score - 
  mutant_parent_clone$parent_averaged_score

escape_parent_clone <- 
  subset(mutant_parent_clone,
         mutant_parent_clone$mutant_challenging_phage == "phitm35" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm36" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm38" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm39" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm40" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm41" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm42" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm46" | 
           mutant_parent_clone$mutant_challenging_phage == "phitm47")

escape_clone_intraclade2 <- 
  subset_intraclade(escape_parent_clone,
                    "parent_gene_content_clade_compare",
                    "clade2")

escape_clone_interclade <- 
  subset_interclade(escape_parent_clone,
                    "parent_gene_content_clade_compare")

escape_clone_intraclade2_empty <- 
  subset_intraclade(escape_parent_clone,
                    "parent_gene_content_clade_compare",
                    "empty")

# QC: Should equal 0.
nrow(escape_clone_intraclade2_empty)


# Fig. 7c sub-panel 2
plot_tricolor_scatter1(escape_clone_intraclade2,
                       escape_clone_interclade,
                       escape_clone_intraclade2_empty,
                       "parent_averaged_score",
                       "mutant_averaged_score",
                       c(0,6),
                       c(0,6),
                       "parent_vs_escape_infection_scores_crs.pdf")


###
###
###
###
###
###
###
###
###
### 11. Compute immunity profile correlations.
# Import table of averaged reciprocal infection data used to create the
# heatmap in Fig. 5a.
# Data structure:
# Rows: infecting phage
# Columns: defending phage
ave_infection_matrix <- read.csv(RECIPROCAL_INFECTION_TABLE_FILENAME,
                                 sep=",",
                                 header=TRUE,
                                 row.names=1,
                                 na.strings="unspecified")

ave_infection_matrix_t <- as.data.frame(t(ave_infection_matrix))

# No data should be missing in this dataset.
defending_cor <- cor(ave_infection_matrix,
                     method="pearson",
                     use="complete.obs")
challenging_cor <- cor(ave_infection_matrix_t,
                       method="pearson",
                       use="complete.obs")

# Convert to 3-column data frame.
# Resulting table contains reciprocal data and self-comparisons.
defending_cor_df <-melt(defending_cor)
challenging_cor_df <-melt(challenging_cor)

names(defending_cor_df) <- c("phage1",
                             "phage2",
                             "defending_cor")
names(challenging_cor_df) <- c("phage1",
                               "phage2",
                               "challenging_cor")

defending_cor_df$phage1_phage2 <- paste(defending_cor_df$phage1,
                                        "_",
                                        defending_cor_df$phage2,
                                        sep="")

challenging_cor_df$phage1_phage2 <- paste(challenging_cor_df$phage1,
                                          "_",
                                          challenging_cor_df$phage2,
                                          sep="")

defending_cor_df$phage1_phage2 <- as.factor(defending_cor_df$phage1_phage2)
challenging_cor_df$phage1_phage2 <- as.factor(challenging_cor_df$phage1_phage2)

defending_cor_df <- subset(defending_cor_df,
                           select=c("phage1_phage2",
                                    "defending_cor"))
challenging_cor_df <- subset(challenging_cor_df,
                             select=c("phage1_phage2",
                                      "challenging_cor"))

immunity_correlation_data <- merge(defending_cor_df,
                                   challenging_cor_df,
                                   by.x="phage1_phage2",
                                   by.y="phage1_phage2",
                                   all.x = TRUE,
                                   all.y = TRUE)

# This merged dataset contains reciprocal data and self-comparisons,
# and can now be merged with other tables elsewhere in the analysis.
###
###
###
###
###
###
###
###
###
### 12. Analyze misc. phage genome metrics.
# Assess correlation of genomic distance, gcd, stoperator_pwd, etc.
# Aside from the immunity profile correlations, this analysis is 
# independent of any immunity assay data.

# Create list of all phage pairs.
actino1321_phages <- levels(phage_metadata$phageid)

distance_metrics <- expand.grid(phage1 = actino1321_phages,
                                phage2 = actino1321_phages)

distance_metrics$phage1_phage2 <- paste(distance_metrics$phage1,
                                        "_",
                                        distance_metrics$phage2,
                                        sep="")

# This dataset now contains reciprocal and self-comparison data.
# Remove all rows in which phage1 and phage2 are not alphabetically ordered.
# This does NOT retain self-comparisons (e.g. alma_alma).
distance_metrics$alpha_ordered <-
  as.character(distance_metrics$phage1) < as.character(distance_metrics$phage2)

# Now retain only the unique pairwise comparisons.
distance_metrics <- subset(distance_metrics,
                           distance_metrics$alpha_ordered == TRUE)

# Match the genome/gene distance data. Many comparisons will not be matched.
distance_metrics <- merge(distance_metrics,
                          genomic_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          repressor_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          cas4_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          endovii_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          dnapol_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          portal_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          stoperator_pwm_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          immunity_correlation_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)


# Match the phage metadata.
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('phage1',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")
distance_metrics <- merge(distance_metrics,
                          phage_metadata_to_match,
                          by.x="phage1",
                          by.y="phage1_phageid")

phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('phage2',
                                        '_',
                                        names(phage_metadata_to_match),
                                        sep="")
distance_metrics <- merge(distance_metrics,
                          phage_metadata_to_match,
                          by.x="phage2",
                          by.y="phage2_phageid")


# Compute comparisons. Since this dataset is not limited to immunity assays, 
# it contains inter-cluster comparisons.
distance_metrics$cluster_compare <- 
  ifelse(distance_metrics$phage1_cluster ==
           distance_metrics$phage2_cluster,
         as.character(distance_metrics$phage1_cluster),
         "different")

distance_metrics$subcluster_compare <- 
  ifelse(distance_metrics$phage1_subcluster ==
           distance_metrics$phage2_subcluster,
         as.character(distance_metrics$phage1_subcluster),
         "different")

distance_metrics$source_compare <- 
  ifelse(distance_metrics$phage1_source ==
           distance_metrics$phage2_source,
         as.character(distance_metrics$phage1_source),
         "different")

distance_metrics$temperate_empirical_compare <- 
  ifelse(distance_metrics$phage1_temperate_empirical ==
           distance_metrics$phage2_temperate_empirical,
         as.character(distance_metrics$phage1_temperate_empirical),
         "different")

distance_metrics$functional_repressor_compare <- 
  ifelse(distance_metrics$phage1_repressor_functional ==
           distance_metrics$phage2_repressor_functional,
         as.character(distance_metrics$phage1_repressor_functional),
         "different")

distance_metrics$lysogen_type_compare <- 
  ifelse(distance_metrics$phage1_lysogen_type == 
           distance_metrics$phage2_lysogen_type,
         as.character(distance_metrics$phage1_lysogen_type),
         "different")

distance_metrics$integrase_compare <- 
  ifelse(distance_metrics$phage1_pham_integrase == 
           distance_metrics$phage2_pham_integrase,
         as.character(distance_metrics$phage1_pham_integrase),
         "different")

distance_metrics$parb_compare <- 
  ifelse(distance_metrics$phage1_pham_parb ==
           distance_metrics$phage2_pham_parb,
         as.character(distance_metrics$phage1_pham_parb),
         "different")

distance_metrics$repressor_hth_compare <- 
  stringdist(as.character(distance_metrics$phage1_repressor_hth),
             as.character(distance_metrics$phage2_repressor_hth),
             method="hamming")

distance_metrics$repressor_length_full_compare <- 
  abs(distance_metrics$phage1_repressor_length_full - 
        distance_metrics$phage2_repressor_length_full)

distance_metrics$repressor_length_nterm_compare <- 
  abs(distance_metrics$phage1_repressor_length_nterm - 
        distance_metrics$phage2_repressor_length_nterm)

distance_metrics$repressor_length_cterm_compare <- 
  abs(distance_metrics$phage1_repressor_length_cterm - 
        distance_metrics$phage2_repressor_length_cterm)

distance_metrics$subcluster_compare2 <- 
  ifelse(distance_metrics$phage1_subcluster ==
           distance_metrics$phage2_subcluster,
         "same",
         "different")

distance_metrics$host_compare <- 
  ifelse(distance_metrics$phage1_host ==
           distance_metrics$phage2_host,
         as.character(distance_metrics$phage1_host),
         "different")

distance_metrics$gene_content_clade_compare <- 
  ifelse(distance_metrics$phage1_gene_content_clade ==
           distance_metrics$phage2_gene_content_clade,
         as.character(distance_metrics$phage1_gene_content_clade),
         "different")

distance_metrics$gene_content_clade_compare2 <- 
  ifelse(distance_metrics$phage1_gene_content_clade ==
           distance_metrics$phage2_gene_content_clade,
         "same",
         "different")

distance_metrics$cluster_compare <-
  as.factor(distance_metrics$cluster_compare)

distance_metrics$subcluster_compare <-
  as.factor(distance_metrics$subcluster_compare)

distance_metrics$source_compare <-
  as.factor(distance_metrics$source_compare)

distance_metrics$temperate_empirical_compare <-
  as.factor(distance_metrics$temperate_empirical_compare)

distance_metrics$functional_repressor_compare <-
  as.factor(distance_metrics$functional_repressor_compare)

distance_metrics$lysogen_type_compare <-
  as.factor(distance_metrics$lysogen_type_compare)

distance_metrics$integrase_compare <-
  as.factor(distance_metrics$integrase_compare)

distance_metrics$parb_compare <-
  as.factor(distance_metrics$parb_compare)

distance_metrics$host_compare <-
  as.factor(distance_metrics$host_compare)

distance_metrics$subcluster_compare2 <-
  factor(distance_metrics$subcluster_compare2, c("same", "different"))

distance_metrics$gene_content_clade_compare <-
  as.factor(distance_metrics$gene_content_clade_compare)

distance_metrics$gene_content_clade_compare2 <-
  factor(distance_metrics$gene_content_clade_compare2, c("same", "different"))


# Keep only Cluster A Mycobacterium phage data. Do not include Cluster A
# Gordonia phages, and do not include lab-derived mutants.
clusterA_data <- subset(
  distance_metrics,
  distance_metrics$cluster_compare == "A" &
    distance_metrics$source_compare == "environment" &
    distance_metrics$host_compare == "Mycobacterium"
)

clusterA_intraclade2 <-
  subset_intraclade(clusterA_data,
         "gene_content_clade_compare",
         "clade2")

clusterA_clade2_diff <-
  subset(
    clusterA_data,
    clusterA_data$gene_content_clade_compare == "different" &
      (
        clusterA_data$phage1_gene_content_clade == 'clade2' |
          clusterA_data$phage2_gene_content_clade == 'clade2'
      )
  )

# Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
clusterA_intraclade2_empty <-
  subset_intraclade(clusterA_data,
         "gene_content_clade_compare",
         "empty")

# QC: Should equal 0.
nrow(clusterA_intraclade2_empty)


# Fig. 2c
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "nuc_dist",
                      "gcd",
                      c(0,0.5),
                      c(0,1),
                      "nuc_vs_gcd.pdf")


# Fig. 2d
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "gcd",
                      "repressor_full_mafft_dist_uncorrected",
                      c(0,1),
                      c(0,70),
                      "gcd_vs_rep.pdf")


# Fig. 2f sub-panel 1
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "gcd",
                      "stoperator_pwd_dist_euc",
                      c(0,1),
                      c(0,5),
                      "gcd_vs_stop.pdf")


# Fig. 2f sub-panel 2
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "repressor_full_mafft_dist_uncorrected",
                      "stoperator_pwd_dist_euc",
                      c(0,70),
                      c(0,5),
                      "rep_vs_stop.pdf")


# Fig. 9a
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "repressor_nterm_mafft_dist_uncorrected",
                      "repressor_cterm_mafft_dist_uncorrected",
                      c(0,70),
                      c(0,70),
                      "rep_nterm_vs_rep_cterm.pdf",
                      "yes")

# Fig. 5c sub-panel 1
plot_tricolor_scatter2(clusterA_intraclade2,
                       clusterA_clade2_diff,
                       clusterA_intraclade2_empty,
                       "stoperator_pwd_dist_euc",
                       "challenging_cor",
                       c(0,5),
                       c(-1,1),
                       "stop_vs_cor_challenging.pdf")


# Fig. 5c sub-panel 2
plot_tricolor_scatter2(clusterA_intraclade2,
                       clusterA_clade2_diff,
                       clusterA_intraclade2_empty,
                       "stoperator_pwd_dist_euc",
                       "defending_cor",
                       c(0,5),
                       c(-1,1),
                       "stop_vs_cor_defending.pdf")


# Summary subsets and stats.
cluster_a <- 
  subset(phage_metadata,
         phage_metadata$cluster == "A")

cluster_a$subcluster <- factor(cluster_a$subcluster)

cluster_a_env <- 
  subset(cluster_a,
         cluster_a$source == "environment")

cluster_a_env_clade2 <- 
  subset(cluster_a_env,
         cluster_a_env$gene_content_clade == "clade2")

cluster_a_env_clade2$subcluster <- factor(cluster_a_env_clade2$subcluster)

cluster_a_env_clade2_myco <- 
  subset(cluster_a_env_clade2,
         cluster_a_env_clade2$host == "Mycobacterium")

cluster_a_env_clade2_myco_rep <- 
  subset(cluster_a_env_clade2_myco,
         cluster_a_env_clade2_myco$repressor_functional == "yes")

nrow(cluster_a)
nrow(cluster_a_env)
nrow(cluster_a_env_clade2)
nrow(cluster_a_env_clade2_myco)
nrow(cluster_a_env_clade2_myco_rep)
nlevels(cluster_a$subcluster)
nlevels(cluster_a_env_clade2$subcluster)


# Fig. S1a
par(mar=c(4,8,16,20))
boxplot(cluster_a_env_clade2_myco_rep$repressor_length_full,
        las=1,cex.axis=2,ann=FALSE,main=NULL,outline=FALSE,ylim=c(150,250),
        col="light grey")
par(new=TRUE)
stripchart(cluster_a_env_clade2_myco_rep$repressor_length_full,
           vertical=TRUE,las=1,cex.axis=2,pch=16,method="jitter",cex=1,
           ann=FALSE,main=NULL,ylim=c(150,250))
dev.copy(pdf,"rep_sizes.pdf")
dev.off()

###
###
###
###
###
###
###
###
###
### 13. Compare stoperator site data.
# Stoperator site prediction data, containing a list of predicted stoperator
# sites in all 327 Cluster A genomes (including escape mutants) from the
# Actino1321 database, for each of the 327 stoperator PWMs derived from all
# Cluster A genomes. An 88% relative score cutoff was used.
# Data structure:
# 1. tfbs88_stop_site_id
# 2. tfbs88_motif_target
# 3. tfbs88_stoperator_target
# 4. tfbs88_stoperator_motif
# 5. tfbs88_site_start (regardless of orientation)
# 6. tfbs88_site_end (regardless of orientation)
# 7. tfbs88_site_strand2
# 8. tfbs88_site_seq
# 9. tfbs88_site_abs_score
# 10. tfbs88_site_rel_score
# 11. tfbs88_site_self
stoperator_sites <-  read.csv(STOPERATOR_SITES_FILENAME,
                              sep = ",",
                              header = TRUE)

# In the Actino1321 database, there are expected to be 
# 327 PWMs * 327 target genomes = 106,929 combinations. This list only
# contains 58,690 levels. This likely reflects that for TFBS88 dataset, since
# low quality data has been removed (using the <88% cutoff score), not all
# target genomes contain predicted sites for all PWMs.

# Create a list of all possible pairwise combinations.
all_levels <-
  expand.grid(
    tfbs88_stoperator_target = as.character(levels(
      stoperator_sites$tfbs88_stoperator_target
    )),
    tfbs88_stoperator_motif = as.character(levels(
      stoperator_sites$tfbs88_stoperator_motif
    ))
  )

all_levels$tfbs88_motif_target <- paste(all_levels$tfbs88_stoperator_motif,
                                        "_",
                                        all_levels$tfbs88_stoperator_target,
                                        sep="")

all_levels$tfbs88_motif_target <- as.factor(all_levels$tfbs88_motif_target)

# Re-factor tfbs88_motif_target using all possible pairwise combinations.This
# will enable quantification for all combinations, even those that are missing.
stoperator_sites$tfbs88_motif_target <-
  factor(stoperator_sites$tfbs88_motif_target,
         levels(all_levels$tfbs88_motif_target))


# Match up phage metadata to stoperator site data.
metadata_to_match <- phage_metadata

names(metadata_to_match) <- paste("target_",names(metadata_to_match),sep="")

stoperator_sites <- merge(stoperator_sites,
                          metadata_to_match,
                          by.x = "tfbs88_stoperator_target",
                          by.y = "target_phageid",
                          all.x=TRUE)


# Compute the middle coordinate for each stoperator. The start is always the
# smaller of the two coordinates, so just add 6.
stoperator_sites$tfbs88_site_middle <- stoperator_sites$tfbs88_site_start + 6


# Now compute how far each site is from the alignment point.
stoperator_sites$site_pleft_dist <- stoperator_sites$tfbs88_site_middle - 
  stoperator_sites$target_coordinate_pleft

stoperator_sites$site_right_end_dist <- stoperator_sites$tfbs88_site_middle -
  stoperator_sites$target_size

stoperator_sites$site_rep_dist <- stoperator_sites$tfbs88_site_middle - 
  stoperator_sites$target_coordinate_repressor

stoperator_sites$site_center_dist <- stoperator_sites$tfbs88_site_middle - 
  stoperator_sites$target_coordinate_genome_center


# Analyze distribution of endogenous sites within each genome.
# Use only phages that infect Mycobacterium, are in Cluster A, and are
# isolated from the environment.
stops_endogenous <- subset(
  stoperator_sites,
  as.character(stoperator_sites$tfbs88_stoperator_target) ==
    as.character(stoperator_sites$tfbs88_stoperator_motif) &
    stoperator_sites$target_source == "environment" &
    stoperator_sites$target_gene_content_clade == "clade2" &
    stoperator_sites$target_host == "Mycobacterium"
)

stops_endogenous$tfbs88_motif_target <-
  factor(stops_endogenous$tfbs88_motif_target)

stops_endogenous$tfbs88_stoperator_target <-
  factor(stops_endogenous$tfbs88_stoperator_target)

stops_endogenous$tfbs88_stoperator_motif <-
  factor(stops_endogenous$tfbs88_stoperator_motif)

stops_endogenous$target_host <- 
  factor(stops_endogenous$target_host)


# Number of stoperators per genome.
stops_endo_freq <-
  as.data.frame(table(stops_endogenous$tfbs88_stoperator_target))

names(stops_endo_freq) <- c("phage","frequency")


# Fig. 3b - distance from repressor
plot_hist1(stops_endogenous,
           "site_rep_dist",
           2000,
           c(-4000,1000),
           c(0,50),
           "stoperators_near_repressor.pdf")


# Fig. 3d - distance from Pleft
plot_hist1(stops_endogenous,
           "site_pleft_dist",
           2500,
           c(-2000,500),
           c(0,100),
           "stoperators_near_pleft.pdf")


# Summary - 22% of all sites are positioned towards the right end of the genome.
nrow(subset(stops_endogenous, stops_endogenous$site_pleft_dist > -1000)) /
  nrow(stops_endogenous)


# Fig. S1e
par(mar=c(4,8,4,4))
hist(stops_endogenous$site_pleft_dist,
     main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
dev.copy(pdf,'stoperators_across_genome.pdf')
dev.off()


# Fig. S1b
plot_hist1(stops_endo_freq,
           "frequency",
           50,
           c(0,50),
           c(0,15),
           "stoperators_per_genome.pdf")


# Compute site orientation. Tally the number of stoperators on each strand
# of each side of genome center.
stops_left <- subset(stops_endogenous,
                     stops_endogenous$site_center_dist < 0)
stops_right <- subset(stops_endogenous,
                      stops_endogenous$site_center_dist >= 0)

stops_left_for <- subset(stops_left,
                         stops_left$tfbs88_site_strand2 == "forward")
stops_left_rev <- subset(stops_left,
                         stops_left$tfbs88_site_strand2 == "reverse")

stops_right_for <- subset(stops_right,
                          stops_right$tfbs88_site_strand2 == "forward")
stops_right_rev <- subset(stops_right,
                          stops_right$tfbs88_site_strand2 == "reverse")

# QC: Should equal 0.
nrow(stops_endogenous) - nrow(stops_left) - nrow(stops_right)
nrow(stops_left) - nrow(stops_left_for) - nrow(stops_left_rev)
nrow(stops_right) - nrow(stops_right_for) - nrow(stops_right_rev)


# Create frequency tables.
stops_left_freq <-
  as.data.frame(table(stops_left$tfbs88_stoperator_target))
names(stops_left_freq) <- c("phage","left_sites_freq")

stops_left_for_freq <-
  as.data.frame(table(stops_left_for$tfbs88_stoperator_target))
names(stops_left_for_freq) <- c("phage","left_sites_forward_freq")

stops_left_rev_freq <-
  as.data.frame(table(stops_left_rev$tfbs88_stoperator_target))
names(stops_left_rev_freq) <- c("phage","left_sites_reverse_freq")

stops_right_freq <-
  as.data.frame(table(stops_right$tfbs88_stoperator_target))
names(stops_right_freq) <- c("phage","right_sites_freq")

stops_right_for_freq <-
  as.data.frame(table(stops_right_for$tfbs88_stoperator_target))
names(stops_right_for_freq) <- c("phage","right_sites_forward_freq")

stops_right_rev_freq <-
  as.data.frame(table(stops_right_rev$tfbs88_stoperator_target))
names(stops_right_rev_freq) <- c("phage","right_sites_reverse_freq")

# QC: All should have the same number of rows.
nrow(stops_endo_freq)
nrow(stops_left_freq)
nrow(stops_left_for_freq)
nrow(stops_left_rev_freq)
nrow(stops_right_freq)
nrow(stops_right_for_freq)
nrow(stops_right_rev_freq)


# Combine site tally data.
stops_freq_summary <- merge(stops_endo_freq,
                            stops_left_freq,
                            by.x="phage",
                            by.y="phage")
stops_freq_summary <- merge(stops_freq_summary,
                            stops_left_for_freq,
                            by.x="phage",
                            by.y="phage")
stops_freq_summary <- merge(stops_freq_summary,
                            stops_left_rev_freq,
                            by.x="phage",
                            by.y="phage")
stops_freq_summary <- merge(stops_freq_summary,
                            stops_right_freq,
                            by.x="phage",
                            by.y="phage")
stops_freq_summary <- merge(stops_freq_summary,
                            stops_right_for_freq,
                            by.x="phage",
                            by.y="phage")
stops_freq_summary <- merge(stops_freq_summary,
                            stops_right_rev_freq,
                            by.x="phage",
                            by.y="phage")

# QC: all should equal 0.
stops_freq_summary$check_total_sites <-
  stops_freq_summary$frequency - 
  stops_freq_summary$right_sites_freq - 
  stops_freq_summary$left_sites_freq

stops_freq_summary$check_left_sites <-
  stops_freq_summary$left_sites_freq - 
  stops_freq_summary$left_sites_forward_freq - 
  stops_freq_summary$left_sites_reverse_freq

stops_freq_summary$check_right_sites <-
  stops_freq_summary$right_sites_freq - 
  stops_freq_summary$right_sites_forward_freq - 
  stops_freq_summary$right_sites_reverse_freq

summary(stops_freq_summary$check_total_sites)
summary(stops_freq_summary$check_left_sites)
summary(stops_freq_summary$check_right_sites)


# Compute percentage of total sites.
stops_freq_summary$left_sites_forward_percent <-
  stops_freq_summary$left_sites_forward_freq / 
  stops_freq_summary$left_sites_freq

stops_freq_summary$left_sites_reverse_percent <-
  stops_freq_summary$left_sites_reverse_freq / 
  stops_freq_summary$left_sites_freq

stops_freq_summary$right_sites_forward_percent <-
  stops_freq_summary$right_sites_forward_freq / 
  stops_freq_summary$right_sites_freq

stops_freq_summary$right_sites_reverse_percent <-
  stops_freq_summary$right_sites_reverse_freq / 
  stops_freq_summary$right_sites_freq

# QC: all should equal 1.
stops_freq_summary$check_left_percent <-
  stops_freq_summary$left_sites_forward_percent + 
  stops_freq_summary$left_sites_reverse_percent

stops_freq_summary$check_right_percent <-
  stops_freq_summary$right_sites_forward_percent + 
  stops_freq_summary$right_sites_reverse_percent

summary(stops_freq_summary$check_left_percent)
summary(stops_freq_summary$check_right_percent)


# Fig. S1c
par(mar=c(4,8,8,4))
plot(stops_freq_summary$right_sites_reverse_percent,
     stops_freq_summary$left_sites_forward_percent,
     xlim=c(0,1),ylim=c(0,1),
     cex.axis=2,ann=FALSE,main=NULL,las=1,
     col="black",pch=16,cex=2)
abline(0,1)
dev.copy(pdf,"stoperators_percent_syn_orientation.pdf")
dev.off()

###