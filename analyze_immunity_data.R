# R script to perform misc. data analyses for Mavrich & Hatfull, mBio, 2017.
# Travis Mavrich.
# Note: this code merges and analyzes various input data sets 
# prepared by other tools including Python and Excel.
# Distinct analyses and code blocks separated by "###".





### Install dependencies

# The melt function of reshape2 package is needed to convert matrix
# to unique-pair table.
#install.packages("reshape2")
library(reshape2)

# Stringdist is needed to compute hamming distance between strings
# install.packages("stringdist")
library(stringdist)




### Set working directory variables.
# The current code assumes a specific working directory structure.
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
        "repressor336_distance_data.csv",
        sep="")
CAS4_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "cas4311_distance_data.csv",
        sep="")
ENDOVII_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "endovii306_distance_data.csv",
        sep="")
DNAPOL_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "dnapol311_distance_data.csv",
        sep="")
PORTAL_DISTANCE_DATA_FILENAME = 
  paste(DIR_INPUT,
        "portal311_distance_data.csv",
        sep="")
STOPERATOR_PWM_DATA_FILENAME = 
  paste(DIR_INPUT,
        "stoperator_pwm_data.csv",
        sep="")
INFECTION_TABLE_REDUCED_FILENAME = 
  paste(DIR_INPUT,
        "multi_lys_ave_infection_table_reduced.csv",
        sep="")
STOPERATOR_SITES_FILENAME = 
  paste(DIR_INPUT,
        "stoperator_site_predictions.csv",
        sep="")


setwd(DIR_OUTPUT)


###Define functions

# Compute comparison fields
compute_comparisons <- function(table){

  table$subcluster_compare <-
    ifelse(
      table$defending_subcluster == table$challenging_subcluster,
      as.character(table$defending_subcluster),
      "different"
    )
  
  table$source_compare <-
    ifelse(
      table$defending_source == table$challenging_source,
      as.character(table$defending_source),
      "different"
    )
  
  table$temperate_empirical_compare <-
    ifelse(
      table$defending_temperate_empirical == 
        table$challenging_temperate_empirical,
      as.character(table$defending_temperate_empirical),
      "different"
    )
  
  table$functional_repressor_compare <-
    ifelse(
      table$defending_repressor_functional == 
        table$challenging_repressor_functional,
      as.character(table$defending_repressor_functional),
      "different"
    )
  
  table$lysogen_type_compare <-
    ifelse(
      table$defending_lysogen_type == 
        table$challenging_lysogen_type,
      as.character(table$defending_lysogen_type),
      "different"
    )
  
  table$integrase_compare <-
    ifelse(
      table$defending_pham_integrase == table$challenging_pham_integrase,
      as.character(table$defending_pham_integrase),
      "different"
    )
  
  table$parb_compare <-
    ifelse(
      table$defending_pham_parb == table$challenging_pham_parb,
      as.character(table$defending_pham_parb),
      "different"
    )
  
  table$repressor_hth_compare <-
    stringdist(
      as.character(table$defending_repressor_hth),
      as.character(table$challenging_repressor_hth),
      method = "hamming"
    )
  
  table$repressor_length_full_compare <-
    abs(table$defending_repressor_length_full - 
          table$challenging_repressor_length_full)

  table$repressor_length_nterm_compare <-
    abs(table$defending_repressor_length_nterm - 
          table$challenging_repressor_length_nterm)

  table$repressor_length_cterm_compare <-
    abs(
      table$defending_repressor_length_cterm - 
        table$challenging_repressor_length_cterm)
  
  table$gene_content_clade_compare <-
    ifelse(
      table$defending_gene_content_clade == 
        table$challenging_gene_content_clade,
      as.character(table$defending_gene_content_clade),
      "different"
    )
  
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
  
  clone_match_columns <- c('defending_challenging',
                           'defending_phage',
                           'challenging_phage',
                           'averaged_rank6',
                           'nuc_dist',
                           'gcd',
                           'repressor_cterm_mafft_dist_uncorrected',
                           'stoperator_pwd_dist_euc',
                           'gene_content_clade_compare')
  
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
  
  lys_clone_compare$lys_averaged_rank6 <- 
    as.numeric(as.character(lys_clone_compare$lys_averaged_rank6))
  
  lys_clone_compare$clone_averaged_rank6 <- 
    as.numeric(as.character(lys_clone_compare$clone_averaged_rank6))
  
  lys_clone_compare$rank6_diff <- 
    lys_clone_compare$clone_averaged_rank6 - 
    lys_clone_compare$lys_averaged_rank6
  
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

subset_homotypic <- function(table,phage1,phage2){
  table_subset <- subset(table,
                         as.character(table[,phage1]) ==
                           as.character(table[,phage2]))
  return(table_subset)
}

subset_heterotypic <- function(table,phage1,phage2){
  table_subset <- subset(table,
                         as.character(table[,phage1]) !=
                           as.character(table[,phage2]))
  return(table_subset)
}





compute_binned_freq1 <- function(table,bin_num){
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





#Plot to compare infection scores of two groups of phages
plot_tricolor_scatter1 <- function(table1,
                                  table2,
                                  table3,
                                  x_data,
                                  y_data,
                                  x_range,
                                  y_range,
                                  filename){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  table3 <- subset(table3,select = c(x_data,y_data))
  table3 <- na.omit(table3)
  
  #For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))
  print(paste("Number of data points in third table: ",nrow(table3)))

  #Plot data
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


#Plot to compare infection scores against a distance metric
plot_tricolor_scatter2 <- function(table1,
                                   table2,
                                   table3,
                                   x_data,
                                   y_data,
                                   x_range,
                                   y_range,
                                   filename){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  
  table3 <- subset(table3,select = c(x_data,y_data))
  table3 <- na.omit(table3)
  
  #For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))
  print(paste("Number of data points in third table: ",nrow(table3)))

  #Plot data
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
  
  dev.copy(pdf,filename)
  dev.off()
  
  lm_gcd <- lm(table1[,x_data] ~ table1[,y_data],data = table1)
  summary(lm_gcd)
}


#Plot to compare infection scores against a distance metric
#with h=0 line
plot_tricolor_scatter3 <- function(table1,
                                   table2,
                                   table3,
                                   x_data,
                                   y_data,
                                   x_range,
                                   y_range,
                                   filename){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  table3 <- subset(table3,select = c(x_data,y_data))
  table3 <- na.omit(table3)
  
  #For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))
  print(paste("Number of data points in third table: ",nrow(table3)))

  #Plot data
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
  abline(h=0)
  
  dev.copy(pdf,filename)
  dev.off()

  lm_gcd <- lm(table1[,x_data] ~ table1[,y_data],data = table1)
  summary(lm_gcd)
}


#Plot to compare genome metrics
plot_bicolor_scatter1 <- function(table1,
                                  table2,
                                  x_data,
                                  y_data,
                                  x_range,
                                  y_range,
                                  filename){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  #For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))

  
  #Plot data
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
  
  dev.copy(pdf,filename)
  dev.off()
}


#Plot to compare genome metrics with abline
plot_bicolor_scatter2 <- function(table1,
                                  table2,
                                  x_data,
                                  y_data,
                                  x_range,
                                  y_range,
                                  filename){
  
  #For each table, remove all rows with missing values.
  table1 <- subset(table1,select = c(x_data,y_data))
  table1 <- na.omit(table1)
  
  table2 <- subset(table2,select = c(x_data,y_data))
  table2 <- na.omit(table2)
  
  
  #For each table, report number of data points plotted.
  print(paste("Number of data points in first table: ",nrow(table1)))
  print(paste("Number of data points in second table: ",nrow(table2)))

  
  #Plot data
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
  abline(0,1)
  dev.copy(pdf,filename)
  dev.off()
}


#Plot bargraph assessing binned frequencies
plot_bargraph1 <- function(table1,value1,value2,y_range,filename){
  
  par(mar=c(4,8,4,4))
  barplot(table1[,value1],
          names.arg=table1[,value2],
          col='black',ylim=y_range)
  dev.copy(pdf,filename)
  dev.off()
}


#Plot bargraph
plot_bargraph2 <- function(table1,value1,y_range,filename){

  par(mar=c(4,8,8,4))
  barplot(summary(table1[,value1]),
          ylim=y_range,
          col="black",ann=FALSE,main=NULL,las=1)
  dev.copy(pdf,filename)
  dev.off()
}


#Plot histograms
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


###End of functions























### Import datasets

#setwd(DIR_INPUT)


# Expected structure of immunity data:
# "immunity_assay_id" (unique identifier)
# "immunity_set"
# "date"
# "notebook"
# "page"
# "strain" (systematic strain name)
# "prophage"
# "repressor_clone"
# "strain_type" (lysogen or repressor_clone)
# "defending_phage"
# "challenging_phage"
# "assay_type" (multiple_titer or single_titer)
# "lawn_notes"
# "lawn_reliability" (0 = unreliable; 3 = reliable)
# "tested_titer"
# "phage_reliability" (0 = unreliable; 3 = reliable)
# "observed_infection_strength"
# "observed_turbidity"
# "observed_plaque_size"
# "observed_plaques"
# "rank6" (systematically scored infection phenotype)

immunity_data <- read.csv(IMMUNITY_DATA_FILENAME,sep=",",header=TRUE)
immunity_data$immunity_assay_id <- as.factor(immunity_data$immunity_assay_id)
immunity_data$immunity_set <- as.factor(immunity_data$immunity_set)
immunity_data$notebook <- as.factor(immunity_data$notebook)
immunity_data$page <- as.factor(immunity_data$page)
immunity_data$lawn_reliability <- as.factor(immunity_data$lawn_reliability)
immunity_data$phage_reliability <- as.factor(immunity_data$phage_reliability)
immunity_data$rank6 <- as.factor(immunity_data$rank6)

# Several fields contain 'unspecified' which can be converted to NA and
# then need to be re-factored:
# "prophage"
# "repressor_clone"
# "tested_titer"
# "phage_reliability"
immunity_data[immunity_data == "Unspecified"] <- NA
immunity_data[immunity_data == "unspecified"] <- NA
immunity_data$prophage <- factor(immunity_data$prophage)
immunity_data$repressor_clone <- factor(immunity_data$repressor_clone)
immunity_data$phage_reliability <- factor(immunity_data$phage_reliability)
immunity_data$rank6 <- factor(immunity_data$rank6)

# Convert titer to numeric
immunity_data$tested_titer <-
  as.numeric(as.character(immunity_data$tested_titer))


# These fields contain NA's, but these are descriptive columns so no need to
# convert them to Unspecified:
# "observed_infection_strength"
# "observed_turbidity"
# "observed_plaque_size"
# "observed_plaques"


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




# TODO confirm how this data was re-generated.
# Mmash genomic distance data for all actino1321 phages, including all
# reciprocal data and self-comparison data.
# Data structure:
# "phage1_phage2"
# "modified_mash_distance" (whole genome nucleotide distance, nuc_dist)
# "pham_pham_dissimilarity" (gene content dissimilarity, gcd)
genomic_distance_data <- read.csv(GENOMIC_DISTANCE_DATA_FILENAME,
                                  sep=",",
                                  header=TRUE)


#Change column names for better readability
names(genomic_distance_data) <- c("phage1_phage2",
                                "nuc_dist",
                                "gcd")


#Use Actino1321 data, in which escape mutants have been added
#actino1321 phage metadata
# "phageid"
# "host"
# "cluster"
# "subcluster"
# "size"
# "lysogen_type" (extrachromosomal, integration, none)
# "pham_integrase"
# "pham_para" (imported as int)
# "source" (environment or lab)
# "parent"
# "repressor_functional" (yes, no, NA; For Cluster A phages, 
  #is the immunity repressor predicted to be functional?)
# "temperate_empirical" (no, unknown, yes, NA; For Cluster A phages,
  #can a lysogen be generated?)
# "repressor_hth"
# "repressor_length_full" (imported as factor)
# "repressor_length_nterm" (imported as factor)
# "repressor_length_cterm" (imported as factor)
# "pham_parb" (imported as int)
# "gene_content_clade" (clade2 = the "L5 clade")
# "coordinate_pleft" (coordinate for manual alignment)
# "coordinate_repressor" (coordinate for manual alignment)
# "coordinate_genome_center" (coordinate for manual alignment)

phage_metadata <- read.csv(PHAGE_METADATA_FILENAME,sep=",",header=TRUE)
phage_metadata$pham_para <- as.factor(phage_metadata$pham_para)
phage_metadata$pham_parb <- as.factor(phage_metadata$pham_parb)







#Several fields contain variants of 'unspecified'
#subcluster = Unspecified
#pham_integrase = ""
#repressor_functional = "not_applicable"
#temperate_empirical = "not_applicable"

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



# Aactino1321 Immunity Repressor protein distance data, including
# reciprocal data and self-comparison data for 336 homologs. There is no
# data for escape mutant or for parent phages that are natural mutants
# (e.g. d29, misswhite, jeffabunny) with no repressor annotated.
# Data structure:
# "phage1_phage2"
# "repressor_full_mafft_phyml_dist"
# "repressor_nterm_mafft_phyml_dist"
# "repressor_cterm_mafft_phyml_dist"
# "repressor_full_mafft_dist_uncorrected"
# "repressor_nterm_mafft_dist_uncorrected"
# "repressor_cterm_mafft_dist_uncorrected"
repressor336_distance_data <-
  read.csv(REPRESSOR_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Aactino1321 Cas4-family protein mafft distance data, including
# reciprocal data and self-comparison data for 311 homologs present 
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# "phage1_phage2"
# "cas4_mafft_dist_uncorrected"
cas4311_distance_data <-
  read.csv(CAS4_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 EndoVII protein mafft distance data, including
# reciprocal data and self-comparison data for 306 homologs present
# in Cluster A parent phages. There is not data escape mutants.
# Data structure:
# "phage1_phage2"
# "endovii_mafft_dist_uncorrected"
endovii306_distance_data <-
  read.csv(ENDOVII_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 DNA Polymerase protein mafft distance data, including 
# reciprocal data and self-comparison data for 311 homologs present
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# "phage1_phage2"
# "dnapol_mafft_dist_uncorrected"
dnapol311_distance_data <-
  read.csv(DNAPOL_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Actino1321 Portal protein mafft distance data, including
# reciprocal data and self-comparison data for 311 homologs present
# in Cluster A parent phages. There is no data for escape mutants.
# Data structure:
# "phage1_phage2"
# "portal_mafft_dist_uncorrected"
portal311_distance_data <-
  read.csv(PORTAL_DISTANCE_DATA_FILENAME,
           sep = ",",
           header = TRUE)


# Position weight matrix distance data, including reciprocal data
# and self-comparison data for all 327 Cluster A phage genomes 
# (including escape mutants) from Actino1321 database.
# Data structure:
# "phage1"
# "phage2"
# "dist_pearson" (Pairwise distance using Pearson metric)
# "dist_euc" (Pairwise distance using Euclidean metric)
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
# Retain only data for phages present in actino1321
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
                            repressor336_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            cas4311_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            endovii306_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            dnapol311_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,
                            portal311_distance_data,
                            by.x="defending_challenging",
                            by.y="phage1_phage2",
                            all.x=TRUE)


# Match PWM distance data. Contains PWM data for 327 Cluster A phages.
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



# Some immunity assay data was derived from single-titer assays. These will
# not be used for this analysis, so they can be removed.

#TODO add code to remove single-titer data here instead of later?



# QC Summary 
# Number of assays.
nrow(main_immunity_data)


# Histogram of titers to assess the range of titers used, which can be used to
# compute multiplicity of infection.

#setwd(DIR_OUTPUT)

par(mar=c(4,8,8,4))
hist(log(main_immunity_data$tested_titer,10),xlim=c(0,10),
     main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",
     breaks=20)
dev.copy(pdf,"tested_titers.pdf")
dev.off()


### At this point, all data in main_immunity_data is derived from phages 
# that are present in the Actino1321 database AND only confident data.





### Average data



# For most analyses, the averaged non-duplicated immunity data is needed. 
# Averages can be computed by unique experiment_id.
main_immunity_data$experiment_id <-
  factor(main_immunity_data$experiment_id)
main_immunity_data$rank6 <- as.numeric(as.character(main_immunity_data$rank6))






#Average infection scores for each unique experiment_id identifier

immunity_average <- aggregate(main_immunity_data[,'rank6'],
                              list(main_immunity_data$experiment_id),mean)
names(immunity_average) <- c('experiment_id',
                             'averaged_rank6') 
immunity_average$experiment_id <- factor(immunity_average$experiment_id)



# Compute the range of scores for each unique assay. First compute the minimum
# score and maximum score, then compute the range.
immunity_min <- aggregate(main_immunity_data[,'rank6'],
                          list(main_immunity_data$experiment_id),min)
names(immunity_min) <- c('experiment_id',
                         'min_rank6') 
immunity_min$experiment_id <- factor(immunity_min$experiment_id)



immunity_max <- aggregate(main_immunity_data[,'rank6'],
                          list(main_immunity_data$experiment_id),max)
names(immunity_max) <- c('experiment_id',
                         'max_rank6') 
immunity_max$experiment_id <- factor(immunity_max$experiment_id)


immunity_average <- merge(immunity_average,
                          immunity_min,
                          by.x="experiment_id",
                          by.y="experiment_id")

immunity_average <- merge(immunity_average,
                          immunity_max,
                          by.x="experiment_id",
                          by.y="experiment_id")


immunity_average$range_rank6 <- immunity_average$max_rank6 -
  immunity_average$min_rank6
immunity_average$range_rank6 <- as.factor(immunity_average$range_rank6)



# Create table of immunity data metadata that will be added back to the
# averaged experiment_id data.
# Keep the following immunity data columns after averaging:
# "experiment_id"
# "defending_challenging"
# "prophage"
# "repressor_clone"
# "strain_type"
# "defending_phage"
# "challenging_phage"
# "assay_type"
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


# Match the phage metadata
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
                          repressor336_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          cas4311_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          endovii306_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          dnapol311_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          portal311_distance_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)
immunity_average <- merge(immunity_average,
                          stoperator_pwm_data,
                          by.x="defending_challenging",
                          by.y="phage1_phage2",
                          all.x=TRUE)


# Compute comparison fields
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





# Export averaged data
#setwd(DIR_OUTPUT)


# All averaged data
write.table(immunity_average,
            paste(DIR_OUTPUT,"immunity_data_averaged.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)




# Table S1
output_fields <- c("strain_type",
                   "defending_phage",
                   "challenging_phage",
                   "nuc_dist",
                   "gcd",
                   "repressor_full_mafft_dist_uncorrected",
                   "stoperator_pwd_dist_euc",
                   "frequency",
                   "min_rank6",
                   "max_rank6",
                   "averaged_rank6")

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

write.table(immunity_average_reduced_for_output,
            paste(DIR_OUTPUT,"Table_S1_averaged_immunity_data.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


# Create average experiment_id data above











# Map averages back to main immunity table, subtract average from original
# values, then plot histogram of differences to show how reliable the
# dataset is.
immunity_ave_to_match <- subset(immunity_average,
                                select = c("experiment_id",
                                           "averaged_rank6",
                                           "min_rank6",
                                           "max_rank6",
                                           "range_rank6",
                                           "frequency"))


#TODO = confirm prefix change works okay
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


main_immunity_data$ave_data_averaged_rank6_diff <-
  as.numeric(as.character(main_immunity_data$rank6)) -
  as.numeric(as.character(main_immunity_data$ave_data_averaged_rank6))
# Note: if the diff is positive = assay's rank is higher than the average.



immunity_average$frequency <-
  as.numeric(as.character(immunity_average$frequency))

immunity_average$frequency2 <-
  ifelse(immunity_average$frequency > 10,
         11,
         immunity_average$frequency)

immunity_average$frequency <- as.factor(immunity_average$frequency)
immunity_average$frequency2 <- as.factor(immunity_average$frequency2)






# Plots
#setwd(DIR_OUTPUT)



#QC: assess how many replicates there are for each unique assay.
plot_bargraph2(immunity_average,
               "frequency",
               c(0,800),
               "immunity_assay_average_frequency.pdf")

plot_bargraph2(immunity_average,
               "frequency2",
               c(0,800),
               "immunity_assay_average_frequency.pdf_adjusted.pdf")



# QC: Number of unique assays with >1 replicate. Note:this includes all
# single-titer assays.
1 - 
  nrow(subset(immunity_average, immunity_average$frequency == "1")) / 
  nrow(immunity_average)
#48% of unique comparisons has > 1 replicate 

# QC: Assess variability in scoring between replicates for each unique assay.
plot_bargraph2(immunity_average,
               "range_rank6",
               c(0,1200),
               "immunity_assay_replicate_range.pdf")
# 92% of unique comparisons with > 1 replicate have a score range of < 2.




# Only look at multiple_titer assays.
immunity_ave_multi <- subset(immunity_average,
                             immunity_average$assay_type == "multiple_titer")



# Summary - number of unique multi-titer assays.
nrow(immunity_ave_multi)

# Summary - number of unique multi-titer assays with >1 replicate.
immunity_ave_multi_reps <- subset(immunity_ave_multi,
                                  immunity_ave_multi$frequency != "1")
nrow(immunity_ave_multi_reps)/nrow(immunity_ave_multi)



# Summary - number of unique multi-titer assays with > 1 replicate and
# score range < 2.
nrow(
  subset(
    immunity_ave_multi_reps,
    immunity_ave_multi_reps$range_rank6 == "0" |
      immunity_ave_multi_reps$range_rank6 == "1"
  )
) / nrow(immunity_ave_multi_reps)


###End of data average step











###Compute immunity profile correlation coefficients 
#setwd(DIR_INPUT)





# Import table of averaged data manipulated in Excel. This is a reduced dataset
# consisting of a complete reciprocal matrix.
infection_table_reduced <- read.csv(INFECTION_TABLE_REDUCED_FILENAME,
                                    sep=",",
                                    header=TRUE,
                                    row.names=1,
                                    na.strings="unspecified")


infection_table_reduced_t <- as.data.frame(t(infection_table_reduced))



# No data should be missing in this dataset.
defending_cor_reduced <- cor(infection_table_reduced,
                             method="pearson",
                             use="complete.obs")
challenging_cor_reduced <- cor(infection_table_reduced_t,
                               method="pearson",
                               use="complete.obs")



# Convert to 3-column data frame
# Resulting table contains reciprocal data and self-comparisons
defending_cor_reduced_df <-melt(defending_cor_reduced)
challenging_cor_reduced_df <-melt(challenging_cor_reduced)


names(defending_cor_reduced_df) <- c("phage1",
                                     "phage2",
                                     "defending_cor_reduced")
names(challenging_cor_reduced_df) <- c("phage1",
                                       "phage2",
                                       "challenging_cor_reduced")

defending_cor_reduced_df$phage1_phage2 <- 
  paste(defending_cor_reduced_df$phage1,
        "_",
        defending_cor_reduced_df$phage2,
        sep="")

challenging_cor_reduced_df$phage1_phage2 <- 
  paste(challenging_cor_reduced_df$phage1,
        "_",
        challenging_cor_reduced_df$phage2,
        sep="")


defending_cor_reduced_df$phage1_phage2 <-
  as.factor(defending_cor_reduced_df$phage1_phage2)

challenging_cor_reduced_df$phage1_phage2 <-
  as.factor(challenging_cor_reduced_df$phage1_phage2)

defending_cor_reduced_df <- subset(defending_cor_reduced_df,
                                   select=c("phage1_phage2",
                                            "defending_cor_reduced"))
challenging_cor_reduced_df <- subset(challenging_cor_reduced_df,
                                     select=c("phage1_phage2",
                                              "challenging_cor_reduced"))






# Merge datasets




# If not using cor_all_df data
immunity_correlation_data <- merge(defending_cor_reduced_df,
                                   challenging_cor_reduced_df,
                                   by.x="phage1_phage2",
                                   by.y="phage1_phage2",
                                   all.x = TRUE,
                                   all.y = TRUE)

# This merged dataset contains reciprocal data and self-comparisons,
# and can now be merged with other tables elsewhere in the analysis.




###Compute immunity profile correlation coefficients above




















###Plot immunity phenotypes by genome metrics to assess all data






# Below: averaged data
#setwd(DIR_OUTPUT)

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


# Split data by clade and by homotypic/heterotypic.
immunity_ave_subset1_interclade <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$gene_content_clade_compare == "different"
)

immunity_ave_subset1_intraclade2 <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$gene_content_clade_compare == "clade2")


immunity_ave_subset1_intraclade2_homotypic <-
  subset(
    immunity_ave_subset1_intraclade2,
    as.character(immunity_ave_subset1_intraclade2$defending_phage) ==
      as.character(immunity_ave_subset1_intraclade2$challenging_phage)
  )

immunity_ave_subset1_intraclade2_heterotypic <-
  subset(
    immunity_ave_subset1_intraclade2,
    as.character(immunity_ave_subset1_intraclade2$defending_phage) !=
      as.character(immunity_ave_subset1_intraclade2$challenging_phage)
  )












# Plots

# Fig. 5b sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_vs_infection_score.pdf")


# Fig. 5b sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_environment_stop_vs_infection_score.pdf")


# Fig. S5a sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "gcd",
                       "averaged_rank6",
                       c(0,1),
                       c(0,6),
                       "lysogen_environment_gcd_vs_infection_score.pdf")


# Fig. S5a sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "nuc_dist",
                       "averaged_rank6",
                       c(0,0.5),
                       c(0,6),
                       "lysogen_environment_nuc_vs_infection_score.pdf")


# Fig. S5d sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "portal_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_portal_vs_infection_score.pdf")


# Fig. S5d sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_pol_vs_infection_score.pdf")


# Fig. S5d sub-panel 3
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_endovii_vs_infection_score.pdf")


# Fig. S5d sub-panel 4
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_cas4_vs_infection_score.pdf")


# Fig. 9b sub-panel 1
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_nterm_vs_infection_score.pdf")


# Fig. 9b sub-panel 2
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "repressor_hth_compare",
                       "averaged_rank6",
                       c(0,10),
                       c(0,6),
                       "lysogen_environment_rep_hth_vs_infection_score.pdf")


# Fig. 9b sub-panel 3
plot_tricolor_scatter2(immunity_ave_subset1_intraclade2_heterotypic,
                       immunity_ave_subset1_interclade,
                       immunity_ave_subset1_intraclade2_homotypic,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "lysogen_environment_rep_cterm_vs_infection_score.pdf")








# Compare differences between prophage maintenance strategy
int_int <- subset(immunity_ave_subset1,
                  immunity_ave_subset1$lysogen_type_compare == 'integration')

extra_extra <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'extrachromosomal'
)

int_extra <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'different' &
    immunity_ave_subset1$defending_lysogen_type == 'integration' &
    immunity_ave_subset1$challenging_lysogen_type == 'extrachromosomal'
)

extra_int <- subset(
  immunity_ave_subset1,
  immunity_ave_subset1$lysogen_type_compare == 'different' &
    immunity_ave_subset1$defending_lysogen_type == 'extrachromosomal' &
    immunity_ave_subset1$challenging_lysogen_type == 'integration'
)


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


# Integrating and extrachromosomal
int_extra_interclade <- subset_interclade(int_extra,
                                          "gene_content_clade_compare")
int_extra_intraclade2 <- subset_intraclade(int_extra,
                                           "gene_content_clade_compare",
                                           "clade2")
int_extra_intraclade2_homotypic <- subset_homotypic(int_extra_intraclade2,
                                                    "defending_phage",
                                                    "challenging_phage")
int_extra_intraclade2_heterotypic <- subset_heterotypic(int_extra_intraclade2,
                                                        "defending_phage",
                                                        "challenging_phage")


# QC
nrow(int_extra)
nrow(int_extra_interclade)
nrow(int_extra_intraclade2)
nrow(int_extra_intraclade2_homotypic)
nrow(int_extra_intraclade2_heterotypic)




# Extrachromosomal and integrating
extra_int_interclade <- subset_interclade(extra_int,
                                          "gene_content_clade_compare")
extra_int_intraclade2 <- subset_intraclade(extra_int,
                                           "gene_content_clade_compare",
                                           "clade2")
extra_int_intraclade2_homotypic <- subset_homotypic(extra_int_intraclade2,
                                                    "defending_phage",
                                                    "challenging_phage")
extra_int_intraclade2_heterotypic <- subset_heterotypic(extra_int_intraclade2,
                                                        "defending_phage",
                                                        "challenging_phage")

# QC
nrow(extra_int)
nrow(extra_int_interclade)
nrow(extra_int_intraclade2)
nrow(extra_int_intraclade2_homotypic)
nrow(extra_int_intraclade2_heterotypic)
#





# Both integrating and same integrase pham.
int_int_same_interclade <- subset_interclade(int_int_same,
                                             "gene_content_clade_compare")
int_int_same_intraclade2 <- subset_intraclade(int_int_same,
                                              "gene_content_clade_compare",
                                              "clade2")

int_int_same_intraclade2_homotypic <-
  subset_homotypic(int_int_same_intraclade2,
                   "defending_phage",
                   "challenging_phage")

int_int_same_intraclade2_heterotypic <- 
  subset_heterotypic(
    int_int_same_intraclade2,
    "defending_phage",
    "challenging_phage")


# QC
nrow(int_int_same)
nrow(int_int_same_interclade)
nrow(int_int_same_intraclade2)
nrow(int_int_same_intraclade2_homotypic)
nrow(int_int_same_intraclade2_heterotypic)
#







# Both integrating but different integrase phams.
int_int_diff_interclade <- subset_interclade(int_int_diff,
                                             "gene_content_clade_compare")
int_int_diff_intraclade2 <- subset_intraclade(int_int_diff,
                                              "gene_content_clade_compare",
                                              "clade2")

int_int_diff_intraclade2_homotypic <-
  subset_homotypic(int_int_diff_intraclade2,
                   "defending_phage",
                   "challenging_phage")

int_int_diff_intraclade2_heterotypic <-
  subset_heterotypic(int_int_diff_intraclade2,
                     "defending_phage",
                     "challenging_phage")


# QC
nrow(int_int_diff)
nrow(int_int_diff_interclade)
nrow(int_int_diff_intraclade2)
nrow(int_int_diff_intraclade2_homotypic)
nrow(int_int_diff_intraclade2_heterotypic)

#







# Both extrachromosomal and same ParB pham.
extra_extra_same_interclade <- subset_interclade(extra_extra_same,
                                                 "gene_content_clade_compare")
extra_extra_same_intraclade2 <- subset_intraclade(extra_extra_same,
                                                  "gene_content_clade_compare",
                                                  "clade2")

extra_extra_same_intraclade2_homotypic <-
  subset_homotypic(extra_extra_same_intraclade2,
                   "defending_phage",
                   "challenging_phage")

extra_extra_same_intraclade2_heterotypic <-
  subset_heterotypic(extra_extra_same_intraclade2,
                     "defending_phage",
                     "challenging_phage")


# QC
nrow(extra_extra_same)
nrow(extra_extra_same_interclade)
nrow(extra_extra_same_intraclade2)
nrow(extra_extra_same_intraclade2_homotypic)
nrow(extra_extra_same_intraclade2_heterotypic)
#









# Both extrachromosomal but different ParB phams.
extra_extra_diff_interclade <- subset_interclade(extra_extra_diff,
                                                 "gene_content_clade_compare")
extra_extra_diff_intraclade2 <- subset_intraclade(extra_extra_diff,
                                                  "gene_content_clade_compare",
                                                  "clade2")

extra_extra_diff_intraclade2_homotypic <-
  subset_homotypic(extra_extra_diff_intraclade2,
                   "defending_phage",
                   "challenging_phage")


extra_extra_diff_intraclade2_heterotypic <-
  subset_heterotypic(extra_extra_diff_intraclade2,
                     "defending_phage",
                     "challenging_phage")



# QC
nrow(extra_extra_diff)
nrow(extra_extra_diff_interclade)
nrow(extra_extra_diff_intraclade2)
nrow(extra_extra_diff_intraclade2_homotypic)
nrow(extra_extra_diff_intraclade2_heterotypic)
#








# Plots

# Fig. S5b sub-panel 1 - IntInt_IntPhamSame
plot_tricolor_scatter2(int_int_same_intraclade2_heterotypic,
                       int_int_same_interclade,
                       int_int_same_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_same_stop_vs_infection_score.pdf")


# Fig. S5b sub-panel 2 - IntInt_IntPhamDiff
plot_tricolor_scatter2(int_int_diff_intraclade2_heterotypic,
                       int_int_diff_interclade,
                       int_int_diff_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_diff_stop_vs_infection_score.pdf")


# Fig. S5b sub-panel 3 - IntExtra
plot_tricolor_scatter2(int_extra_intraclade2_heterotypic,
                       int_extra_interclade,
                       int_extra_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_int_parb_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 1 - ExtraExtra_ParBPhamSame
plot_tricolor_scatter2(extra_extra_same_intraclade2_heterotypic,
                       extra_extra_same_interclade,
                       extra_extra_same_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_same_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 2 - ExtraExtra_ParBPhamDiff
plot_tricolor_scatter2(extra_extra_diff_intraclade2_heterotypic,
                       extra_extra_diff_interclade,
                       extra_extra_diff_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_diff_stop_vs_infection_score.pdf")


# Fig. S5c sub-panel 3 - ExtraInt
plot_tricolor_scatter2(extra_int_intraclade2_heterotypic,
                       extra_int_interclade,
                       extra_int_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "lysogen_env_parb_int_stop_vs_infection_score.pdf")








#Binned Data Stats

#Plot only Clade 2 OR (Clade 2 non-homotypic) immunity data stats


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



#TODO might be able to add some of these steps fo the compute_bin_freq function if it doesn't interfere with other usage of the function
clade2_env_bin0 <- subset(clade2_env,
                          clade2_env$averaged_rank6 <= 0.5)
clade2_env_bin1 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 0.5 &
                            clade2_env$averaged_rank6 <= 1.5)
clade2_env_bin2 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 1.5 &
                            clade2_env$averaged_rank6 <= 2.5)
clade2_env_bin3 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 2.5 &
                            clade2_env$averaged_rank6 <= 3.5)
clade2_env_bin4 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 3.5 &
                            clade2_env$averaged_rank6 <= 4.5)
clade2_env_bin5 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 4.5 &
                            clade2_env$averaged_rank6 <= 5.5)
clade2_env_bin6 <- subset(clade2_env,
                          clade2_env$averaged_rank6 > 5.5)


clade2_env_bin0$defending_phage <- factor(clade2_env_bin0$defending_phage)
clade2_env_bin1$defending_phage <- factor(clade2_env_bin1$defending_phage)
clade2_env_bin2$defending_phage <- factor(clade2_env_bin2$defending_phage)
clade2_env_bin3$defending_phage <- factor(clade2_env_bin3$defending_phage)
clade2_env_bin4$defending_phage <- factor(clade2_env_bin4$defending_phage)
clade2_env_bin5$defending_phage <- factor(clade2_env_bin5$defending_phage)
clade2_env_bin6$defending_phage <- factor(clade2_env_bin6$defending_phage)

clade2_env_bin0$challenging_phage <- factor(clade2_env_bin0$challenging_phage)
clade2_env_bin1$challenging_phage <- factor(clade2_env_bin1$challenging_phage)
clade2_env_bin2$challenging_phage <- factor(clade2_env_bin2$challenging_phage)
clade2_env_bin3$challenging_phage <- factor(clade2_env_bin3$challenging_phage)
clade2_env_bin4$challenging_phage <- factor(clade2_env_bin4$challenging_phage)
clade2_env_bin5$challenging_phage <- factor(clade2_env_bin5$challenging_phage)
clade2_env_bin6$challenging_phage <- factor(clade2_env_bin6$challenging_phage)








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








#Plots

#Fig. S4a sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "defending_percent",
               "bin",
               c(0,1),
               "infection_score_percent_defending.pdf")


#Fig. S4a sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "challenging_percent",
               "bin",
               c(0,1),
               "infection_score_percent_challenging.pdf")


#Fig. S4a sub-panel 3
plot_bargraph1(clade2_binned_frequency,
               "total_assays_percent",
               "bin",
               c(0,0.3),
               "infection_score_percent_total.pdf")


#Fig. S4b sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "inter_subcluster_percent",
               "bin",
               c(0,0.5),
               "infection_score_percent_intersubcluster.pdf")


#Fig. S4b sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "intra_subcluster_percent",
               "bin",
               c(0,0.5),
               "infection_score_percent_intrasubcluster.pdf")



###Above: averaged data



















###Reciprocal analysis
# Using averaged data, compute differences in reciprocal immunity data.
# Split immunity_average columns into groups.


# Phage-specific data = metadata that is not impacted by immunity vector:
# "prophage"
# "repressor_clone"
# "strain_type"
# "defending_host"
# "defending_cluster"
# "defending_subcluster"
# "defending_size"
# "defending_lysogen_type"
# "defending_pham_integrase"
# "defending_pham_para"
# "defending_source"
# "defending_parent"
# "defending_repressor_functional"
# "defending_temperate_empirical"
# "defending_repressor_hth"
# "defending_repressor_length_full"
# "defending_repressor_length_nterm"
# "defending_repressor_length_cterm"
# "defending_pham_parb"
# "challenging_host"
# "challenging_cluster"
# "challenging_subcluster"
# "challenging_size"
# "challenging_lysogen_type"
# "challenging_pham_integrase"
# "challenging_pham_para"
# "challenging_source"
# "challenging_parent"
# "challenging_repressor_functional"
# "challenging_temperate_empirical"
# "challenging_repressor_hth"
# "challenging_repressor_length_full"
# "challenging_repressor_length_nterm"
# "challenging_repressor_length_cterm"
# "challenging_pham_parb"


# Phage metadata comparisons.Data specific to both phages used in immunity 
# but not impacted by vector orientation:
# "nuc_dist"
# "gcd"
# "repressor_muscle_bionj_distances"
# "repressor_prank_phyml_distances"
# "portal_muscle_bionj_distances"
# "portal_prank_phyml_distances"
# "recb_muscle_bionj_distances"
# "recb_prank_phyml_distances"
# "repressor_full_mafft_phyml_dist"
# "repressor_nterm_mafft_phyml_dist"
# "repressor_cterm_mafft_phyml_dist"
# "repressor_full_mafft_dist_uncorrected"
# "repressor_nterm_mafft_dist_uncorrected"
# "repressor_cterm_mafft_dist_uncorrected"
# "stoperator_pwd_dist_pearson"
# "stoperator_pwd_dist_euc"
# "source_compare"
# "temperate_empirical_compare"
# "functional_repressor_compare"
# "lysogen_type_compare"
# "integrase_compare"
# "repressor_hth_compare"
# "repressor_length_full_compare"
# "repressor_length_nterm_compare"
# "repressor_length_cterm_compare"
# "parb_compare"
# "subcluster_compare"
# "gene_content_clade_compare"



# Vectored_data = data that is specific to immunity assay
# and depends on vector orientation:
vectored_column_names <- c("experiment_id",
                           "defending_challenging",
                           "defending_phage",
                           "challenging_phage",
                           "averaged_rank6",
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
reciprocal_data$averaged_rank6_diff <-
  abs(reciprocal_data$vector1_averaged_rank6 - 
        reciprocal_data$vector2_averaged_rank6)


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



# Output the reciprocal dataset.
#setwd(DIR_OUTPUT)
write.table(reciprocal_data_alpha_ordered,
            paste(DIR_OUTPUT,"reciprocal_immunity_data.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)









# Bin reciprocal data for histogram


# Reduce the reciprocal dataset to only environmental phages.
reciprocal_unique_envY <- 
  subset(reciprocal_data,
         reciprocal_data$vector1_alpha_ordered == TRUE &
           reciprocal_data$source_compare == 'environment')


# Bin data
reciprocal_bin0 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff <= 0.5)
reciprocal_bin1 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 0.5 &
                            reciprocal_unique_envY$averaged_rank6_diff <= 1.5)
reciprocal_bin2 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 1.5 &
                            reciprocal_unique_envY$averaged_rank6_diff <= 2.5)
reciprocal_bin3 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 2.5 &
                            reciprocal_unique_envY$averaged_rank6_diff <= 3.5)
reciprocal_bin4 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 3.5 &
                            reciprocal_unique_envY$averaged_rank6_diff <= 4.5)
reciprocal_bin5 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 4.5 &
                            reciprocal_unique_envY$averaged_rank6_diff <= 5.5)
reciprocal_bin6 <- subset(reciprocal_unique_envY,
                          reciprocal_unique_envY$averaged_rank6_diff > 5.5)





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


# QC: should equal 0.
sum(reciprocal_binned_freq$freq) - nrow(reciprocal_unique_envY)


# Plot data
#setwd(DIR_OUTPUT)




# Fig. S4c
plot_bargraph1(reciprocal_binned_freq,
               "freq_percent",
               "bin",
               c(0,0.5),
               "reciprocal_infection_scores_percent_total.pdf")








# How do various genome/gene distances correlate with differences 
# in reciprocal profiles?




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

reciprocal_intraclade2_homotypic <- 
  subset_homotypic(reciprocal_intraclade2,
                   "vector1_defending_phage",
                   "vector1_challenging_phage")

reciprocal_intraclade2_heterotypic <- 
  subset_heterotypic(reciprocal_intraclade2,
                     "vector1_defending_phage",
                     "vector1_challenging_phage")





#Plots

#Fig. 5b sub-panel 3
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_vs_infection_score.pdf")


#Fig. 5b sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6_diff",
                       c(0,5),
                       c(0,6),
                       "reciprocal_lys_env_stop_vs_infection_score.pdf")


#Fig. S5a sub-panel 3
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "gcd",
                       "averaged_rank6_diff",
                       c(0,1),
                       c(0,6),
                       "reciprocal_lys_env_gcd_vs_infection_score.pdf")


#Fig. S5a sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "nuc_dist",
                       "averaged_rank6_diff",
                       c(0,0.5),
                       c(0,6),
                       "reciprocal_lys_env_nuc_vs_infection_score.pdf")


#Fig. S5d sub-panel 5
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "portal_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_portal_vs_infection_score.pdf")


#Fig. S5d sub-panel 6
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_pol_vs_infection_score.pdf")


#Fig. S5d sub-panel 7
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_endovii_vs_infection_score.pdf")


#Fig. S5d sub-panel 8
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_cas4_vs_infection_score.pdf")


#Fig. 9b sub-panel 4
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_nterm_vs_infection_score.pdf")


#Fig. 9b sub-panel 5
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "repressor_hth_compare",
                       "averaged_rank6_diff",
                       c(0,10),
                       c(0,6),
                       "reciprocal_lys_env_rep_hth_vs_infection_score.pdf")


#Fig. 9b sub-panel 6
plot_tricolor_scatter2(reciprocal_intraclade2_heterotypic,
                       reciprocal_interclade,
                       reciprocal_intraclade2_homotypic,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_lys_env_rep_cterm_vs_infection_score.pdf")



###Reciprocal analysis above














###Compare lysogen vs repressor clone data

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

env_lys_clone_intraclade2_heterotypic <- 
  subset_homotypic(env_lys_clone_intraclade2,
                   "lys_defending_phage",
                   "lys_challenging_phage")

env_lys_clone_intraclade2_homotypic <- 
  subset_heterotypic(env_lys_clone_intraclade2,
                     "lys_defending_phage",
                     "lys_challenging_phage")






# Escape mutant data from averaged multiple-titer lysogen and clone data

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


escape_lys_clone_interclade <- 
  subset_interclade(escape_lys_clone_matched,
                    "lys_gene_content_clade_compare")

escape_lys_clone_intraclade2 <- 
  subset_intraclade(escape_lys_clone_matched,
                    "lys_gene_content_clade_compare",
                    "clade2")

escape_lys_clone_intraclade2_homotypic <- 
  subset_homotypic(escape_lys_clone_intraclade2,
                   "lys_defending_phage",
                   "lys_challenging_phage")

escape_lys_clone_intraclade2_heterotypic <- 
  subset_heterotypic(escape_lys_clone_intraclade2,
                     "lys_defending_phage",
                     "lys_challenging_phage")


# QC
nrow(env_lys_clone_matched)
nrow(escape_lys_clone_matched)

nrow(env_lys_clone_intraclade2_homotypic)
nrow(env_lys_clone_intraclade2_heterotypic)
nrow(env_lys_clone_interclade)

nrow(env_lys_clone_intraclade2_homotypic) +
  nrow(env_lys_clone_intraclade2_heterotypic) + 
  nrow(env_lys_clone_interclade)


nrow(escape_lys_clone_intraclade2_heterotypic)
nrow(escape_lys_clone_interclade)
nrow(escape_lys_clone_intraclade2_homotypic)

nrow(escape_lys_clone_intraclade2_heterotypic) + 
  nrow(escape_lys_clone_interclade) + 
  nrow(escape_lys_clone_intraclade2_homotypic)

#





#Plots
#setwd(DIR_OUTPUT)


#Fig. 6c
plot_tricolor_scatter1(env_lys_clone_intraclade2_homotypic,
                       env_lys_clone_interclade,
                       env_lys_clone_intraclade2_heterotypic,
                       "lys_averaged_rank6","clone_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "lys_vs_crs_env_infection_score.pdf")


#Fig. 7d
plot_tricolor_scatter1(escape_lys_clone_intraclade2_heterotypic,
                       escape_lys_clone_interclade,
                       escape_lys_clone_intraclade2_homotypic,
                       "lys_averaged_rank6","clone_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "lys_vs_crs_escape_infection_score.pdf")


#Fig. 6d
plot_tricolor_scatter3(env_lys_clone_intraclade2_homotypic,
                       env_lys_clone_interclade,
                       env_lys_clone_intraclade2_heterotypic,
                       "lys_stoperator_pwd_dist_euc",
                       "rank6_diff",
                       c(0,5),
                       c(-4,4),
                       "lys_vs_crs_env_infection_score_stop_vs_score_diff.pdf")


#Fig. 7e
plot_tricolor_scatter3(escape_lys_clone_intraclade2_heterotypic,
                       escape_lys_clone_interclade,
                       escape_lys_clone_intraclade2_homotypic,
                       "lys_stoperator_pwd_dist_euc",
                       "rank6_diff",
                       c(0,5),
                       c(-4,4),
                       "lys_vs_crs_escape_infection_score_stop_vs_score_diff.pdf")


###Lysogen-clone comparison above












###L5-mutant comparisons using all averaged data



# All averaged multiple-titer lysogen and clone data.
l5_columns <- c('defending_challenging',
                'defending_phage',
                'challenging_phage',
                'averaged_rank6',
                'nuc_dist',
                'gcd',
                'repressor_cterm_mafft_dist_uncorrected')

immunity_ave_subset3 <- 
  subset(immunity_average,
         immunity_average$assay_type == 'multiple_titer' &
           immunity_average$strain_type == 'lysogen')


# Compare L5 and L5-derivative superinfection patterns.
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


chal_l5_assays$l5_phitm41_rank6_diff <-
  chal_l5_assays$phitm41_averaged_rank6 - chal_l5_assays$l5_averaged_rank6

chal_l5_assays$l5_phitm1_rank6_diff <-
  chal_l5_assays$phitm1_averaged_rank6 - chal_l5_assays$l5_averaged_rank6

chal_l5_assays$l5_phitm4_rank6_diff <-
  chal_l5_assays$phitm4_averaged_rank6 - chal_l5_assays$l5_averaged_rank6

chal_l5_assays$l5_phitm6_rank6_diff <-
  chal_l5_assays$phitm6_averaged_rank6 - chal_l5_assays$l5_averaged_rank6

chal_l5_assays$phitm1_phitm6_rank6_diff <-
  chal_l5_assays$phitm6_averaged_rank6 - chal_l5_assays$phitm1_averaged_rank6

chal_l5_assays$phitm1_phitm4_rank6_diff <-
  chal_l5_assays$phitm4_averaged_rank6 - chal_l5_assays$phitm1_averaged_rank6


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

chal_l5_assays_intraclade2_phitm1_homotypic <- 
  subset_homotypic(chal_l5_assays_intraclade2,
                   "phitm1_challenging_phage",
                   "l5_defending_phage")

chal_l5_assays_intraclade2_phitm1_heterotypic <- 
  subset_heterotypic(chal_l5_assays_intraclade2,
                     "phitm1_challenging_phage",
                     "l5_defending_phage")


# QC
nrow(chal_l5_assays_intraclade2)
nrow(subset(chal_l5_assays_intraclade2,
            is.na(chal_l5_assays_intraclade2$phitm1_challenging_phage)))
nrow(chal_l5_assays_intraclade2_phitm1_homotypic)
nrow(chal_l5_assays_intraclade2_phitm1_heterotypic)
nrow(chal_l5_assays_intraclade2_empty)
#







#Plots
#setwd(DIR_OUTPUT)


#Fig. 10d sub-panel 1
plot_tricolor_scatter2(chal_l5_assays_intraclade2_phitm1_heterotypic,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_phitm1_homotypic,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm1_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_phitm1.pdf")


#Fig. 10d sub-panel 2
plot_tricolor_scatter2(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm4_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_phitm4.pdf")


#Fig. S8g sub-panel 1
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_rank6",
                       "phitm41_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm41.pdf")


#Fig. S8g sub-panel 2
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_rank6",
                       "phitm6_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm6.pdf")


#Fig. S8g sub-panel 3
plot_tricolor_scatter1(chal_l5_assays_intraclade2,
                       chal_l5_assays_interclade,
                       chal_l5_assays_intraclade2_empty,
                       "l5_averaged_rank6",
                       "phitm1_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "infection_scores_l5_vs_phitm1.pdf")



#defending correlation


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


def_l5_assays$l5_phitm41_rank6_diff <-
  def_l5_assays$phitm41_averaged_rank6 - def_l5_assays$l5_averaged_rank6

def_l5_assays$l5_phitm1_rank6_diff <-
  def_l5_assays$phitm1_averaged_rank6 - def_l5_assays$l5_averaged_rank6

def_l5_assays$l5_phitm6_rank6_diff <-
  def_l5_assays$phitm6_averaged_rank6 - def_l5_assays$l5_averaged_rank6

def_l5_assays$phitm1_phitm6_rank6_diff <-
  def_l5_assays$phitm6_averaged_rank6 - def_l5_assays$phitm1_averaged_rank6


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


# QC:
nrow(def_l5_assays)
nrow(def_l5_assays_intraclade2)
nrow(def_l5_assays_interclade)
nrow(def_l5_assays_clade2_empty)

summary(def_l5_assays$l5_challenging_phage)
summary(def_l5_assays_intraclade2$l5_challenging_phage)
summary(def_l5_assays_interclade$l5_challenging_phage)






#Plots

#Fig. S8h sub-panel 1
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm41_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm41.pdf")


#Fig. S8h sub-panel 2
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm6_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm6.pdf")


#Fig. S8h sub-panel 3
plot_tricolor_scatter1(def_l5_assays_intraclade2,
                       def_l5_assays_interclade,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm1_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "defense_scores_l5_vs_phitm1.pdf")


###L5-mutant comparisons using all averaged data above











### Mutant-parent superinfection profile comparison
# Pair superinfecting parent and superinfecting mutant data
# Split averaged data into parent and mutant phage tables

# Averaged, non-redundant data should only be for multi-titer assays,
# no low-confidence data, and only from lysogen data. There should be no
# repressor clones or single-titer assays.

mutant_data <- subset(immunity_average,
                      immunity_average$strain_type == 'lysogen' & 
                        immunity_average$assay_type == 'multiple_titer')

# Rename all fields to indicate the data is generated using the mutant
# challenger. Note: important to remember that "mutant" prefix refers to the
# entire immunity assay data, not to any particular column.
names(mutant_data) <- paste('mutant_',names(mutant_data),sep="")



# Drop all data from mutant table not involving a mutant challenging phage. 
# It is possible that a lab-derived phage is not an escape mutant, so this
# a broad list of mutants.
mutant_data <- subset(mutant_data,
                      mutant_data$mutant_challenging_source == 'lab')




# Create parent data to match
parent_data <- subset(immunity_average,
                      immunity_average$strain_type == 'lysogen' & 
                        immunity_average$assay_type == 'multiple_titer')

names(parent_data) <- paste('parent_',names(parent_data),sep="")


# Now match parent challenging data to mutant challenging data.
mutant_data$parent_defending_challenging <- 
  paste(mutant_data$mutant_defending_phage,
        "_",
        mutant_data$mutant_challenging_parent
        ,sep="")

mutant_analysis <- merge(mutant_data,
                         parent_data,
                         by.x="parent_defending_challenging",
                         by.y="parent_defending_challenging")





# Compute difference in infection profiles. The direction of change in
# infection between the mutant and parent is informative, so don't use
# absolute value.
mutant_analysis$averaged_rank6_diff <-
  mutant_analysis$mutant_averaged_rank6 - mutant_analysis$parent_averaged_rank6



#setwd(DIR_OUTPUT)
write.table(mutant_analysis,
            paste(DIR_OUTPUT,"mutant_analysis.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)








# QC
plot_hist1(mutant_analysis,
           "averaged_rank6_diff",
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


escape_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm35" | 
           mutant_analysis$mutant_challenging_phage == "phitm36" | 
           mutant_analysis$mutant_challenging_phage == "phitm38" | 
           mutant_analysis$mutant_challenging_phage == "phitm39" | 
           mutant_analysis$mutant_challenging_phage == "phitm40" | 
           mutant_analysis$mutant_challenging_phage == "phitm41" | 
           mutant_analysis$mutant_challenging_phage == "phitm42" | 
           mutant_analysis$mutant_challenging_phage == "phitm46" | 
           mutant_analysis$mutant_challenging_phage == "phitm47")


phitm35_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm35")

phitm36_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm36")

phitm38_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm38")

phitm39_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm39")

phitm40_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm40")

phitm41_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm41")

phitm42_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm42")

phitm46_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm46")

phitm47_mutant_analysis <- 
  subset(mutant_analysis,
         mutant_analysis$mutant_challenging_phage == "phitm47")


escape_interclade <- subset_interclade(escape_mutant_analysis,
                                       "parent_gene_content_clade_compare")

escape_intraclade2 <- subset_intraclade(escape_mutant_analysis,
                                        "parent_gene_content_clade_compare",
                                        "clade2")

#Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
escape_intraclade2_empty <- 
  subset_intraclade(escape_mutant_analysis,
                    "parent_gene_content_clade_compare",
                    "empty")

escape_intraclade2_parent_homotypic <- subset_homotypic(escape_intraclade2,
                                              "parent_challenging_phage",
                                                "parent_defending_phage")

escape_intraclade2_parent_heterotypic <- subset_heterotypic(escape_intraclade2,
                                                "parent_challenging_phage",
                                                  "parent_defending_phage")

escape_intraclade2_mutant_homotypic <- subset_homotypic(escape_intraclade2,
                                              "mutant_challenging_phage",
                                                "mutant_defending_phage")

escape_intraclade2_mutant_heterotypic <- subset_heterotypic(escape_intraclade2,
                                                "mutant_challenging_phage",
                                                  "mutant_defending_phage")

# QC
nrow(escape_intraclade2)
nrow(escape_intraclade2_empty)
nrow(escape_intraclade2_parent_homotypic)
nrow(escape_intraclade2_parent_heterotypic)
nrow(escape_intraclade2_mutant_homotypic)
nrow(escape_intraclade2_mutant_heterotypic)


summary(escape_intraclade2_parent_homotypic$parent_defending_phage)
summary(escape_intraclade2_mutant_homotypic$mutant_defending_phage)
#







#Plots


#Fig. 7c sub-panel 1
plot_tricolor_scatter1(escape_intraclade2,
                       escape_interclade,
                       escape_intraclade2_empty,
                       "parent_averaged_rank6",
                       "mutant_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "parent_vs_escape_infection_scores_lysogen.pdf")


#Fig. 7f sub-panel 1
plot_tricolor_scatter2(escape_intraclade2_parent_heterotypic,
                       escape_interclade,
                       escape_intraclade2_parent_homotypic,
                       "parent_stoperator_pwd_dist_euc",
                       "parent_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_parent.pdf")


#Fig. 7f sub-panel 2
plot_tricolor_scatter2(escape_intraclade2_mutant_heterotypic,
                       escape_interclade,
                       escape_intraclade2_mutant_homotypic,
                       "mutant_stoperator_pwd_dist_euc",
                       "mutant_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "stop_vs_infection_score_escape.pdf")


#TODO change all references to 'homotypic' to 'same'
#TODO change all references to 'heterotypic' to 'diff'



# Compare escape mutants to parents on repressor clone strains. Use same
# pipeline as for lysogen strains.

#TODO change _rep to _clone
mutant_data_rep <- subset(immunity_average,
                      immunity_average$strain_type == 'repressor_clone' & 
                        immunity_average$assay_type == 'multiple_titer')

names(mutant_data_rep) <- paste('mutant_',
                                names(mutant_data_rep),
                                sep="")

mutant_data_rep <- subset(mutant_data_rep,
                          mutant_data_rep$mutant_challenging_source == 'lab')

parent_data_rep <- subset(immunity_average,
                      immunity_average$strain_type == 'repressor_clone' & 
                        immunity_average$assay_type == 'multiple_titer')

names(parent_data_rep) <- paste('parent_',names(parent_data_rep),sep="")



mutant_data_rep$parent_defending_challenging <- 
  paste(mutant_data_rep$mutant_defending_phage,
        "_",
        mutant_data_rep$mutant_challenging_parent,
        sep="")

mutant_rep_analysis <- merge(mutant_data_rep,
                             parent_data_rep,
                             by.x="parent_defending_challenging",
                             by.y="parent_defending_challenging")

mutant_rep_analysis$averaged_rank6_diff <-
  mutant_rep_analysis$mutant_averaged_rank6 - 
  mutant_rep_analysis$parent_averaged_rank6

escape_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm35" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm36" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm38" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm39" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm40" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm41" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm42" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm46" | 
           mutant_rep_analysis$mutant_challenging_phage == "phitm47")


phitm35_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm35")
phitm36_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm36")
phitm38_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm38")
phitm39_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm39")
phitm40_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm40")
phitm41_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm41")
phitm42_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm42")
phitm46_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm46")
phitm47_mutant_rep_analysis <- 
  subset(mutant_rep_analysis,
         mutant_rep_analysis$mutant_challenging_phage == "phitm47")


escape_clone_intraclade2 <- 
  subset_intraclade(escape_mutant_rep_analysis,
                    "parent_gene_content_clade_compare",
                    "clade2")

escape_clone_interclade <- 
  subset_interclade(escape_mutant_rep_analysis,
                    "parent_gene_content_clade_compare")

escape_clone_intraclade2_empty <- 
  subset_intraclade(escape_mutant_rep_analysis,
                    "parent_gene_content_clade_compare",
                    "empty")


# QC
nrow(escape_clone_intraclade2_empty)





#Plots


#Fig. 7c sub-panel 2
plot_tricolor_scatter1(escape_clone_intraclade2,
                       escape_clone_interclade,
                       escape_clone_intraclade2_empty,
                       "parent_averaged_rank6",
                       "mutant_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "parent_vs_escape_infection_scores_crs.pdf")


###Compare known empirical temperate to isolated mutant to escape mutant


















###Whole genome metrics
# Assess correlation of genomic distance, gcd, stoperator_pwd, etc.
# This should be independent of any immunity assay data





# Create list of all phage pairs  
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



#Match the genome/gene distance data. Many comparisons will not be matched
distance_metrics <- merge(distance_metrics,
                          genomic_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          repressor336_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          cas4311_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          endovii306_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          dnapol311_distance_data,
                          by.x="phage1_phage2",
                          by.y="phage1_phage2",
                          all.x=TRUE)
distance_metrics <- merge(distance_metrics,
                          portal311_distance_data,
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



#Match the phage metadata
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



# Compute comparisons


# Since this dataset is not limited to immunity assays, 
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

distance_metrics$subcluster_compare2 <-
  factor(distance_metrics$subcluster_compare2, c("same", "different"))

distance_metrics$gene_content_clade_compare <-
  as.factor(distance_metrics$gene_content_clade_compare)

distance_metrics$gene_content_clade_compare2 <-
  factor(distance_metrics$gene_content_clade_compare2, c("same", "different"))








# Plot relationship between repressors, stoperators, and genomes







# Keep only Mycobacteriophage data and drop Gordonia phages.
# Do not include escape mutants.
clusterA_data <- subset(
  distance_metrics,
  distance_metrics$cluster_compare == "A" &
    distance_metrics$source_compare == "environment" &
    distance_metrics$phage1_host == "Mycobacterium" &
    distance_metrics$phage2_host == "Mycobacterium"
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

#QC
nrow(clusterA_intraclade2_empty)








#Plots
#setwd(DIR_OUTPUT)


#Fig. 2c
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "nuc_dist",
                      "gcd",
                      c(0,0.5),
                      c(0,1),
                      "nuc_vs_gcd.pdf")


#Fig. 2d
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "gcd",
                      "repressor_full_mafft_dist_uncorrected",
                      c(0,1),
                      c(0,70),
                      "gcd_vs_rep.pdf")


#Fig. 2f sub-panel 1
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "gcd",
                      "stoperator_pwd_dist_euc",
                      c(0,1),
                      c(0,5),
                      "gcd_vs_stop.pdf")


#Fig. 2f sub-panel 2
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "repressor_full_mafft_dist_uncorrected",
                      "stoperator_pwd_dist_euc",
                      c(0,70),
                      c(0,5),
                      "rep_vs_stop.pdf")


#Fig. 9a
plot_bicolor_scatter2(clusterA_clade2_diff,
                      clusterA_intraclade2,
                      "repressor_nterm_mafft_dist_uncorrected",
                      "repressor_cterm_mafft_dist_uncorrected",
                      c(0,70),
                      c(0,70),
                      "rep_nterm_vs_rep_cterm.pdf")


#Fig. 5c sub-panel 1
plot_tricolor_scatter2(clusterA_intraclade2,
                       clusterA_clade2_diff,
                       clusterA_intraclade2_empty,
                       "stoperator_pwd_dist_euc",
                       "challenging_cor_reduced",
                       c(0,5),
                       c(-1,1),
                       "stop_vs_cor_challenging.pdf")


#Fig. 5c sub-panel 2
plot_tricolor_scatter2(clusterA_intraclade2,
                       clusterA_clade2_diff,
                       clusterA_intraclade2_empty,
                       "stoperator_pwd_dist_euc",
                       "defending_cor_reduced",
                       c(0,5),
                       c(-1,1),
                       "stop_vs_cor_defending.pdf")






#Repressor size
clusterA_subset <- subset(
  phage_metadata,
  phage_metadata$cluster == 'A' &
    phage_metadata$source == 'environment' &
    phage_metadata$host == 'Mycobacterium' &
    phage_metadata$repressor_functional == 'yes' &
    phage_metadata$gene_content_clade == 'clade2'
)


#Fig. S1a
par(mar=c(4,8,16,20))
boxplot(clusterA_subset$repressor_length_full,
        las=1,cex.axis=2,ann=FALSE,main=NULL,outline=FALSE,ylim=c(150,250),
        col="light grey")
par(new=TRUE)
stripchart(clusterA_subset$repressor_length_full,
           vertical=TRUE,las=1,cex.axis=2,pch=16,method="jitter",cex=1,
           ann=FALSE,main=NULL,ylim=c(150,250))
dev.copy(pdf,"rep_sizes.pdf")
dev.off()



### Whole genome metrics (regardless of immunity data) above












# Stoperator site prediction data, containing a list of predicted stoperator
# sites in all 327 Cluster A genomes (including escape mutants) from the
# Actino1321 database, for each of the 327 stoperator PWMs derived from all
# Cluster A genomes. An 88% relative score cutoff was used.
# Data structure:
# "tfbs88_stop_site_id"
# "tfbs88_motif_target"
# "tfbs88_stoperator_target"
# "tfbs88_stoperator_motif"
# "tfbs88_site_start" (regardless of orientation)
# "tfbs88_site_end" (regardless of orientation)
# "tfbs88_site_strand2"
# "tfbs88_site_seq"
# "tfbs88_site_abs_score"
# "tfbs88_site_rel_score"
# "tfbs88_site_self"

#setwd(DIR_INPUT)


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



# Number of stoperators per genome
stops_endo_freq <-
  as.data.frame(table(stops_endogenous$tfbs88_stoperator_target))

names(stops_endo_freq) <- c("phage","frequency")








# Distance from Pleft
#setwd(DIR_OUTPUT)



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


# QC: should equal 0
nrow(stops_endogenous) - nrow(stops_left) - nrow(stops_right)
nrow(stops_left) - nrow(stops_left_for) - nrow(stops_left_rev)
nrow(stops_right) - nrow(stops_right_for) - nrow(stops_right_rev)






#Create frequency tables
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


#QC: all shoudl have the same number of rows.
nrow(stops_endo_freq)
nrow(stops_left_freq)
nrow(stops_left_for_freq)
nrow(stops_left_rev_freq)
nrow(stops_right_freq)
nrow(stops_right_for_freq)
nrow(stops_right_rev_freq)



# Combine site tally data
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







# Plots
#setwd(DIR_OUTPUT)




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

###Above stoperator site analysis






