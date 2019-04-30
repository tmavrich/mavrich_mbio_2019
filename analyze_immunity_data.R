# R script to perform misc. data analyses for Mavrich & Hatfull, mBio, 2017.
# Travis Mavrich.
# Note: this code merges and analyzes various datasets generated from 
# input data sets prepared by other tools including Python and Excel.
# Distinct analyses and code blocks separated by "###".





### Install dependencies

# TODO check if reshape2 is still needed
# The melt function of reshape2 package is needed to convert matrix
# to unique-pair table.
# install.packages("reshape2")
# TODO: update R and RStudio to most appropriate version for reshape2
library(reshape2)

# Stringdist is needed to compute hamming distance between strings
# install.packages("stringdist")
library(stringdist)



### Set working directory
# The current code assumes a specific working directory structure.

setwd("~/scratch/immunity_analysis/input/")









###Define functions

# Compute comparison fields
compute_comparisons <- function(table){

  table$subcluster_compare <- ifelse(table$defending_subcluster == table$challenging_subcluster,
                                     as.character(table$defending_subcluster),"different")
  
  table$source_compare <- ifelse(table$defending_source==table$challenging_source,
                                        as.character(table$defending_source),
                                        "different")
  
  table$temperate_empirical_compare <- ifelse(table$defending_cluster_a_temperate_empirical==table$challenging_cluster_a_temperate_empirical,
                                                     as.character(table$defending_cluster_a_temperate_empirical),
                                                     "different")
  
  table$functional_repressor_compare <- ifelse(table$defending_cluster_a_functional_repressor_predicted==table$challenging_cluster_a_functional_repressor_predicted,
                                                     as.character(table$defending_cluster_a_functional_repressor_predicted),
                                                     "different")
  
  table$lysogen_type_compare <- ifelse(table$defending_lysogen_type==table$challenging_lysogen_type,
                                              as.character(table$defending_lysogen_type),
                                              "different")
  
  table$integrase_compare <- ifelse(table$defending_pham_integrase==table$challenging_pham_integrase,
                                           as.character(table$defending_pham_integrase),
                                           "different")

  table$parb_compare <- ifelse(table$defending_pham_parb==table$challenging_pham_parb,
                                    as.character(table$defending_pham_parb),
                                    "different")
  
  table$repressor_hth_compare <- stringdist(as.character(table$defending_repressor_hth_domain_sequence),
                                            as.character(table$challenging_repressor_hth_domain_sequence),
                                            method="hamming")
  
  table$repressor_length_full_compare <- abs(table$defending_repressor_length_full - table$challenging_repressor_length_full)
  table$repressor_length_nterm_compare <- abs(table$defending_repressor_length_nterm - table$challenging_repressor_length_nterm)
  table$repressor_length_cterm_compare <- abs(table$defending_repressor_length_cterm - table$challenging_repressor_length_cterm)
  
  table$gene_content_clade_compare <- ifelse(table$defending_gene_content_clade==table$challenging_gene_content_clade,
                               as.character(table$defending_gene_content_clade),
                               "different")
  
  table$subcluster_compare <- as.factor(table$subcluster_compare)
  table$source_compare <- as.factor(table$source_compare)
  table$temperate_empirical_compare <- as.factor(table$temperate_empirical_compare)
  table$functional_repressor_compare <- as.factor(table$functional_repressor_compare)
  table$lysogen_type_compare <- as.factor(table$lysogen_type_compare)
  table$integrase_compare <- as.factor(table$integrase_compare)
  table$parb_compare <- as.factor(table$parb_compare)  
  table$gene_content_clade_compare <- as.factor(table$gene_content_clade_compare)  

  return(table)
}


match_average_lysogen_and_clone_data <- function(lysogen_data,clone_data){

  lysY_reduced <- subset(lysogen_data,select = c(clone_match_columns))
  names(lysY_reduced) <- paste('lys_',names(lysY_reduced),sep="")
  
  cloneY_reduced <- subset(clone_data,select = c(clone_match_columns))
  names(cloneY_reduced) <- paste('clone_',names(cloneY_reduced),sep="")

  lys_clone_compare <- merge(lysY_reduced,cloneY_reduced,by.x='lys_defending_challenging',by.y='clone_defending_challenging')
  lys_clone_compare$lys_averaged_rank6 <- as.numeric(as.character(lys_clone_compare$lys_averaged_rank6))
  lys_clone_compare$clone_averaged_rank6 <- as.numeric(as.character(lys_clone_compare$clone_averaged_rank6))
  lys_clone_compare$rank6_diff <- lys_clone_compare$clone_averaged_rank6 - lys_clone_compare$lys_averaged_rank6

  return(lys_clone_compare)
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



# Expected structure of immunity data:
# "immunity_assay_id"
# "immunity_set"
# "date"
# "notebook"
# "page"
# "strain"
# "prophage"
# "repressor_clone"
# "strain_type"
# "defending_phage"
# "challenging_phage"
# "assay_type"
# "lawn_notes"
# "lawn_reliability"
# "tested_titer"
# "phage_reliability"
# "observed_infection_strength"
# "observed_turbidity"
# "observed_plaque_size"
# "observed_plaques"
# "rank6"

immunity_data <- read.csv("immunity_data.csv",sep=",",header=TRUE)
immunity_data$immunity_assay_id <- as.factor(immunity_data$immunity_assay_id)
immunity_data$immunity_set <- as.factor(immunity_data$immunity_set)
immunity_data$notebook <- as.factor(immunity_data$notebook)
immunity_data$page <- as.factor(immunity_data$page)
immunity_data$lawn_reliability <- as.factor(immunity_data$lawn_reliability)
immunity_data$phage_reliability <- as.factor(immunity_data$phage_reliability)
immunity_data$rank6 <- as.factor(immunity_data$rank6)

# Several fields contain 'unspecified' which can be converted to NA
# Afterwards, re-factor
# prophage
# repressor_clone
# tested_titer
# phage_reliability
immunity_data[immunity_data == "Unspecified"] <- NA
immunity_data[immunity_data == "unspecified"] <- NA
immunity_data$prophage <- factor(immunity_data$prophage)
immunity_data$repressor_clone <- factor(immunity_data$repressor_clone)
immunity_data$phage_reliability <- factor(immunity_data$phage_reliability)
immunity_data$rank6 <- factor(immunity_data$rank6)

#Convert titer to numeric
immunity_data$tested_titer <- as.numeric(as.character(immunity_data$tested_titer))


#These fields contain NA's, but these are descriptive columns so no need to conver them to Unspecified
#observed_infection_strength
#observed_turbidity
#observed_plaque_size
#observed_plaques


#Create identifiers
immunity_data$defending_challenging <- paste(immunity_data$defending_phage,"_",immunity_data$challenging_phage,sep="")
immunity_data$defending_challenging <- as.factor(immunity_data$defending_challenging)
immunity_data$assay_strain_defending_challenging <- paste(immunity_data$assay_type,
                                                          "_",
                                                          immunity_data$strain_type,
                                                          "_",
                                                          immunity_data$defending_phage,
                                                          "_",
                                                          immunity_data$challenging_phage,
                                                          sep="")


immunity_data$assay_strain_defending_challenging <- as.factor(immunity_data$assay_strain_defending_challenging)


immunity_data$strain_defending_challenging <- paste(immunity_data$strain_type,
                                                    "_",
                                                    immunity_data$defending_phage,
                                                    "_",
                                                    immunity_data$challenging_phage,
                                                    sep="")

immunity_data$strain_defending_challenging <- as.factor(immunity_data$strain_defending_challenging)



# TODO confirm how this data was re-generated.
#mash genomic distance data
#all actino1321 phages 
#reciprocal data and self-comparison data
# "phage1_phage2"
# "modified_mash_distance"
# "pham_pham_dissimilarity"
genomic_distance_data <- read.csv("genomic_distance_data.csv",sep=",",header=TRUE)



#Use Actino1321 data, in which escape mutants have been added
#actino1321 phage metadata
# "phageid"
# "host"
# "cluster"
# "subcluster"
# "size"
# "lysogen_type"
# "pham_integrase"
# "pham_para" (imported as int)
# "source"
# "parent"
# "cluster_a_functional_repressor_predicted"
# "cluster_a_temperate_empirical"
# "repressor_hth_domain_sequence"
# "repressor_length_full" (imported as factor)
# "repressor_length_nterm" (imported as factor)
# "repressor_length_cterm" (imported as factor)
# "pham_parb" (imported as int)
# "gene_content_clade"
# "pleft_alignment_reference"
# "immunity_repressor_alignment_reference"
# "genome_center_alignment_reference"

phage_metadata <- read.csv("phage_metadata.csv",sep=",",header=TRUE)
phage_metadata$pham_para <- as.factor(phage_metadata$pham_para)
phage_metadata$pham_parb <- as.factor(phage_metadata$pham_parb)



#Several fields contain variants of 'unspecified'
#subcluster = Unspecified
#pham_integrase = ""
#cluster_a_functional_repressor_predicted = "not_applicable"
#cluster_a_temperate_empirical = "not_applicable"

phage_metadata[phage_metadata == "Unspecified"] <- NA
phage_metadata[phage_metadata == "unspecified"] <- NA
phage_metadata[phage_metadata == ""] <- NA
phage_metadata[phage_metadata == "not_applicable"] <- NA
phage_metadata$subcluster <- factor(phage_metadata$subcluster)
phage_metadata$pham_integrase <- factor(phage_metadata$pham_integrase)
phage_metadata$cluster_a_functional_repressor_predicted <- factor(phage_metadata$cluster_a_functional_repressor_predicted)
phage_metadata$cluster_a_temperate_empirical <- factor(phage_metadata$cluster_a_temperate_empirical)
phage_metadata$gene_content_clade <- factor(phage_metadata$gene_content_clade)
phage_metadata$repressor_length_full <- as.numeric(as.character(phage_metadata$repressor_length_full))
phage_metadata$repressor_length_nterm <- as.numeric(as.character(phage_metadata$repressor_length_nterm))
phage_metadata$repressor_length_cterm <- as.numeric(as.character(phage_metadata$repressor_length_cterm))
phage_metadata$pleft_alignment_reference <- as.numeric(as.character(phage_metadata$pleft_alignment_reference))
phage_metadata$immunity_repressor_alignment_reference <- as.numeric(as.character(phage_metadata$immunity_repressor_alignment_reference))
phage_metadata$genome_center_alignment_reference <- as.numeric(as.character(phage_metadata$genome_center_alignment_reference))



#actino1321 cluster a repressor protein distance data
#reciprocal data and self-comparison data
#336 repressors
#no data for escape mutants
#no data for parent phages that are natural mutants (e.g. d29, misswhite, jeffabunny) with no repressor annotated
# "phage1_phage2"
# "repressor_full_mafft_phyml_dist"
# "repressor_nterm_mafft_phyml_dist"
# "repressor_cterm_mafft_phyml_dist"
# "repressor_full_mafft_dist_uncorrected"
# "repressor_nterm_mafft_dist_uncorrected"
# "repressor_cterm_mafft_dist_uncorrected"
repressor336_distance_data <- read.csv("repressor336_distance_data.csv",sep=",",header=TRUE)


#actino1321 cas4 protein mafft distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "cas4_mafft_dist_uncorrected"
cas4311_distance_data <- read.csv("cas4311_distance_data.csv",sep=",",header=TRUE)


#actino1321 endoVII protein mafft distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "endovii_mafft_dist_uncorrected"
endovii306_distance_data <- read.csv("endovii306_distance_data.csv",sep=",",header=TRUE)


#actino1321 DNA Polymerase protein mafft distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "dnapol_mafft_dist_uncorrected"
dnapol311_distance_data <- read.csv("dnapol311_distance_data.csv",sep=",",header=TRUE)


#actino1321 portal protein mafft distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "portal_mafft_dist_uncorrected"
portal311_distance_data <- read.csv("portal311_distance_data.csv",sep=",",header=TRUE)


#stoperator327 position weight matrix data
#contains reciprocal data and self-comparison data
#stoperator248 = does NOT contain data for escape mutants
#stoperator264 = DOES contain data for escape mutants
#stoperator327 = contains stoperator data for all 327 Cluster A phage genomes (including escape mutants) from actino1321 database
#phage1
#phage2
#dist_pearson
#dist_euc
stoperator_pwm_data <- read.csv("stoperator_pwm_data.csv",sep=",",header=TRUE)
stoperator_pwm_data$phage1_phage2 <- paste(stoperator_pwm_data$phage1,
                                           "_",
                                           stoperator_pwm_data$phage2,
                                           sep="")
stoperator_pwm_data$phage1_phage2 <- as.factor(stoperator_pwm_data$phage1_phage2)
names(stoperator_pwm_data) <- c("phage1",
                                "phage2",
                                "stoperator_pwd_dist_pearson",
                                "stoperator_pwd_dist_euc",
                                "phage1_phage2")

stoperator_pwm_data <- subset(stoperator_pwm_data,
                              select=c("phage1_phage2",
                                       "stoperator_pwd_dist_pearson",
                                       "stoperator_pwd_dist_euc"))


#stoperator327 site prediction
#contains list of prediction stoperator sites in all 327 Cluster A genomes from actino1321, for each of the 327 stoperaotr PWMs
#An 88% relative score cutoff was used
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
stoperator_site_predictions88 <- read.csv("stoperator_site_predictions.csv",sep=",",header=TRUE)



#actino1319 = 248 PWMs * 325 target genomes = 80,600 combinations. This list only contains 80,165 levels.
#This likely reflects the fact that target genomes may not contain any predicted sites.

#actino1321 = 264 PWMs * 327 target genomes = 86,238 combinations. 
#actino1321 = 327 PWMs * 327 target genomes = 106,929 combinations. This list only contains 58690 levels.
#This likely reflects that for TFBS88 dataset, since I have removed a lot of low quality data (<88% score),
#not all target genomes contain predicted sites for all PWMs.


all_levels <- expand.grid(tfbs88_stoperator_target = as.character(levels(stoperator_site_predictions88$tfbs88_stoperator_target)),
                          tfbs88_stoperator_motif = as.character(levels(stoperator_site_predictions88$tfbs88_stoperator_motif)))

all_levels$tfbs88_motif_target <- paste(all_levels$tfbs88_stoperator_motif,
                                 "_",
                                 all_levels$tfbs88_stoperator_target,
                                 sep="")

all_levels$tfbs88_motif_target <- as.factor(all_levels$tfbs88_motif_target)


#Re-factor tfbs88_motif_target using all possible pairwise combinations.
#This will enable quantification for all combinations, even those that are missing.
stoperator_site_predictions88$tfbs88_motif_target <- factor(stoperator_site_predictions88$tfbs88_motif_target,levels(all_levels$tfbs88_motif_target))







#Match up phage metadata to stoperator site data
metadata_to_match <- phage_metadata

names(metadata_to_match) <- paste("target_",names(metadata_to_match),sep="")

stoperator_site_predictions88 <- merge(stoperator_site_predictions88,
                                       metadata_to_match,
                                       by.x = "tfbs88_stoperator_target",
                                       by.y = "target_phageid",
                                       all.x=TRUE)



#Compute the middle coordinate for each stoperator.
#The start is always the smaller of the two coordinates, so just add 6.
stoperator_site_predictions88$tfbs88_site_middle <- stoperator_site_predictions88$tfbs88_site_start + 6


#Now compute how far each site is from the alignment point
stoperator_site_predictions88$site_pleft_dist <- stoperator_site_predictions88$tfbs88_site_middle - stoperator_site_predictions88$target_pleft_alignment_reference
stoperator_site_predictions88$site_right_end_dist <- stoperator_site_predictions88$tfbs88_site_middle - stoperator_site_predictions88$target_size
stoperator_site_predictions88$site_rep_dist <- stoperator_site_predictions88$tfbs88_site_middle - stoperator_site_predictions88$target_immunity_repressor_alignment_reference
stoperator_site_predictions88$site_center_dist <- stoperator_site_predictions88$tfbs88_site_middle - stoperator_site_predictions88$target_genome_center_alignment_reference


#Analyze distribution of sites

stop88_clade2_env_self <- subset(stoperator_site_predictions88,
                                 as.character(stoperator_site_predictions88$tfbs88_stoperator_target) == as.character(stoperator_site_predictions88$tfbs88_stoperator_motif) &
                                   stoperator_site_predictions88$target_source == "environment" &
                                   stoperator_site_predictions88$target_gene_content_clade == "clade2" &
                                   stoperator_site_predictions88$target_host == "Mycobacterium")


stop88_clade2_env_self$tfbs88_motif_target <- factor(stop88_clade2_env_self$tfbs88_motif_target)
stop88_clade2_env_self$tfbs88_stoperator_target <- factor(stop88_clade2_env_self$tfbs88_stoperator_target)
stop88_clade2_env_self$tfbs88_stoperator_motif <- factor(stop88_clade2_env_self$tfbs88_stoperator_motif)
stop88_clade2_env_self$target_host <- factor(stop88_clade2_env_self$target_host)

















#Distance from Pleft
setwd("~/scratch/immunity_analysis/output/")


#Fig. S1e
par(mar=c(4,8,4,4))
hist(stop88_clade2_env_self$site_pleft_dist,
     main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
dev.copy(pdf,'clade2_env_self_sites_near_pleft2.pdf')
dev.off()


#Fig. 3d
plot_hist1(stop88_clade2_env_self,
           "site_pleft_dist",
           2500,
           c(-2000,500),
           c(0,100),
           "clade2_env_self_sites_near_pleft.pdf")

#Summary
nrow(subset(stop88_clade2_env_self,stop88_clade2_env_self$site_pleft_dist > -1000))/nrow(stop88_clade2_env_self)






#Distance from repressor


#Fig. 3b
plot_hist1(stop88_clade2_env_self,
           "site_rep_dist",
           2000,
           c(-4000,1000),
           c(0,50),
           "clade2_env_self_sites_near_repressor2.pdf")


#Examine # of stoperators per genome
stop88_clade2_env_self_freq <- as.data.frame(table(stop88_clade2_env_self$tfbs88_stoperator_target))
names(stop88_clade2_env_self_freq) <- c("phage","frequency")


#Fig. S1b
plot_hist1(stop88_clade2_env_self_freq,
           "frequency",
           50,
           c(0,50),
           c(0,15),
           "clade2_env_self_stops_per_genome.pdf")












#Site orientation
#Tally # of stoperators on each strand of each side of genome center
stop88_clade2_env_self_left_sites <- subset(stop88_clade2_env_self,
                                            stop88_clade2_env_self$site_center_dist < 0)

stop88_clade2_env_self_right_sites <- subset(stop88_clade2_env_self,
                                             stop88_clade2_env_self$site_center_dist >= 0)




stop88_clade2_env_self_left_sites_forward <- subset(stop88_clade2_env_self_left_sites,
                                                    stop88_clade2_env_self_left_sites$tfbs88_site_strand2 == "forward")
stop88_clade2_env_self_left_sites_reverse <- subset(stop88_clade2_env_self_left_sites,
                                                    stop88_clade2_env_self_left_sites$tfbs88_site_strand2 == "reverse")

stop88_clade2_env_self_right_sites_forward <- subset(stop88_clade2_env_self_right_sites,
                                                     stop88_clade2_env_self_right_sites$tfbs88_site_strand2 == "forward")
stop88_clade2_env_self_right_sites_reverse <- subset(stop88_clade2_env_self_right_sites,
                                                     stop88_clade2_env_self_right_sites$tfbs88_site_strand2 == "reverse")


#QC: should equal 0
nrow(stop88_clade2_env_self) - nrow(stop88_clade2_env_self_left_sites) - nrow(stop88_clade2_env_self_right_sites)
nrow(stop88_clade2_env_self_left_sites) - nrow(stop88_clade2_env_self_left_sites_forward) - nrow(stop88_clade2_env_self_left_sites_reverse)
nrow(stop88_clade2_env_self_right_sites) - nrow(stop88_clade2_env_self_right_sites_forward) - nrow(stop88_clade2_env_self_right_sites_reverse)






#Create frequency tables
stop88_clade2_env_self_left_sites_freq <- as.data.frame(table(stop88_clade2_env_self_left_sites$tfbs88_stoperator_target))
names(stop88_clade2_env_self_left_sites_freq) <- c("phage","left_sites_freq")

stop88_clade2_env_self_left_sites_forward_freq <- as.data.frame(table(stop88_clade2_env_self_left_sites_forward$tfbs88_stoperator_target))
names(stop88_clade2_env_self_left_sites_forward_freq) <- c("phage","left_sites_forward_freq")

stop88_clade2_env_self_left_sites_reverse_freq <- as.data.frame(table(stop88_clade2_env_self_left_sites_reverse$tfbs88_stoperator_target))
names(stop88_clade2_env_self_left_sites_reverse_freq) <- c("phage","left_sites_reverse_freq")



stop88_clade2_env_self_right_sites_freq <- as.data.frame(table(stop88_clade2_env_self_right_sites$tfbs88_stoperator_target))
names(stop88_clade2_env_self_right_sites_freq) <- c("phage","right_sites_freq")

stop88_clade2_env_self_right_sites_forward_freq <- as.data.frame(table(stop88_clade2_env_self_right_sites_forward$tfbs88_stoperator_target))
names(stop88_clade2_env_self_right_sites_forward_freq) <- c("phage","right_sites_forward_freq")

stop88_clade2_env_self_right_sites_reverse_freq <- as.data.frame(table(stop88_clade2_env_self_right_sites_reverse$tfbs88_stoperator_target))
names(stop88_clade2_env_self_right_sites_reverse_freq) <- c("phage","right_sites_reverse_freq")


nrow(stop88_clade2_env_self_freq)
nrow(stop88_clade2_env_self_left_sites_freq)
nrow(stop88_clade2_env_self_left_sites_forward_freq)
nrow(stop88_clade2_env_self_left_sites_reverse_freq)
nrow(stop88_clade2_env_self_right_sites_freq)
nrow(stop88_clade2_env_self_right_sites_forward_freq)
nrow(stop88_clade2_env_self_right_sites_reverse_freq)



#Combine site tally data
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq,
                                             stop88_clade2_env_self_left_sites_freq,
                                             by.x="phage",
                                             by.y="phage")
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq_summary,
                                             stop88_clade2_env_self_left_sites_forward_freq,
                                             by.x="phage",
                                             by.y="phage")
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq_summary,
                                             stop88_clade2_env_self_left_sites_reverse_freq,
                                             by.x="phage",
                                             by.y="phage")
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq_summary,
                                             stop88_clade2_env_self_right_sites_freq,
                                             by.x="phage",
                                             by.y="phage")
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq_summary,
                                             stop88_clade2_env_self_right_sites_forward_freq,
                                             by.x="phage",
                                             by.y="phage")
stop88_clade2_env_self_freq_summary <- merge(stop88_clade2_env_self_freq_summary,
                                             stop88_clade2_env_self_right_sites_reverse_freq,
                                             by.x="phage",
                                             by.y="phage")


#QC should equal 0
stop88_clade2_env_self_freq_summary$check_total_sites <- stop88_clade2_env_self_freq_summary$frequency - stop88_clade2_env_self_freq_summary$right_sites_freq - stop88_clade2_env_self_freq_summary$left_sites_freq
stop88_clade2_env_self_freq_summary$check_left_sites <- stop88_clade2_env_self_freq_summary$left_sites_freq - stop88_clade2_env_self_freq_summary$left_sites_forward_freq - stop88_clade2_env_self_freq_summary$left_sites_reverse_freq
stop88_clade2_env_self_freq_summary$check_right_sites <- stop88_clade2_env_self_freq_summary$right_sites_freq - stop88_clade2_env_self_freq_summary$right_sites_forward_freq - stop88_clade2_env_self_freq_summary$right_sites_reverse_freq
summary(stop88_clade2_env_self_freq_summary$check_total_sites)
summary(stop88_clade2_env_self_freq_summary$check_left_sites)
summary(stop88_clade2_env_self_freq_summary$check_right_sites)






stop88_clade2_env_self_freq_summary$left_sites_percent <- stop88_clade2_env_self_freq_summary$left_sites_freq/stop88_clade2_env_self_freq_summary$frequency
stop88_clade2_env_self_freq_summary$right_sites_percent <- stop88_clade2_env_self_freq_summary$right_sites_freq/stop88_clade2_env_self_freq_summary$frequency


stop88_clade2_env_self_freq_summary$left_sites_forward_percent <- stop88_clade2_env_self_freq_summary$left_sites_forward_freq/stop88_clade2_env_self_freq_summary$left_sites_freq
stop88_clade2_env_self_freq_summary$left_sites_reverse_percent <- stop88_clade2_env_self_freq_summary$left_sites_reverse_freq/stop88_clade2_env_self_freq_summary$left_sites_freq


stop88_clade2_env_self_freq_summary$right_sites_forward_percent <- stop88_clade2_env_self_freq_summary$right_sites_forward_freq/stop88_clade2_env_self_freq_summary$right_sites_freq
stop88_clade2_env_self_freq_summary$right_sites_reverse_percent <- stop88_clade2_env_self_freq_summary$right_sites_reverse_freq/stop88_clade2_env_self_freq_summary$right_sites_freq


#QC
stop88_clade2_env_self_freq_summary$check_left_percent <- stop88_clade2_env_self_freq_summary$left_sites_forward_percent + stop88_clade2_env_self_freq_summary$left_sites_reverse_percent
stop88_clade2_env_self_freq_summary$check_right_percent <- stop88_clade2_env_self_freq_summary$right_sites_forward_percent + stop88_clade2_env_self_freq_summary$right_sites_reverse_percent
summary(stop88_clade2_env_self_freq_summary$check_left_percent)
summary(stop88_clade2_env_self_freq_summary$check_right_percent)





setwd("~/scratch/immunity_analysis/output/")

#Fig. S1c
par(mar=c(4,8,8,4))
plot(stop88_clade2_env_self_freq_summary$right_sites_reverse_percent,
     stop88_clade2_env_self_freq_summary$left_sites_forward_percent,
     xlim=c(0,1),ylim=c(0,1),
     cex.axis=2,ann=FALSE,main=NULL,las=1,
     col="black",pch=16,cex=2)
abline(0,1)
dev.copy(pdf,'clade2_env_self_percent_txn_oriented_sites_per_genome.pdf')
dev.off()






#Match the phage metadata
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('defending','_',names(phage_metadata_to_match),sep="")



#There are 2 defending phages not in actino1321 = 'l5_gp71-flag' and 'redrock_bxb1'. 
#Retain only data for phages present in actino1321
#main_immunity_data <- merge(immunity_data,phage_metadata_to_match,by.x="defending_phage",by.y="defending_phageid",all.x=TRUE)
main_immunity_data <- merge(immunity_data,phage_metadata_to_match,by.x="defending_phage",by.y="defending_phageid")

phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('challenging','_',names(phage_metadata_to_match),sep="")


#There are several challenging phages not in actino1321 = 'd29_mutant','l5_flag1-1','l5_ha1','l5_ha2','l5_ha3','l5_ha3_hazy'.
#Retain only data for phages present in actino1321
#main_immunity_data <- merge(main_immunity_data,phage_metadata_to_match,by.x="challenging_phage",by.y="challenging_phageid",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,phage_metadata_to_match,by.x="challenging_phage",by.y="challenging_phageid")






###At this point, all data in main_immunity_data is derived from phages present in actino1321 database###



#Match the genomic distance data. Defending phages and challenging phages not in actino1321 will not be matched
main_immunity_data <- merge(main_immunity_data,genomic_distance_data,by.x="defending_challenging",by.y="phage1_phage2")





#Match the gene distance data. Many comparisons will not be matched
main_immunity_data <- merge(main_immunity_data,repressor336_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,cas4311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,endovii306_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,dnapol311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,portal311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)

#TODO remove
# main_immunity_data <- merge(main_immunity_data,repressor_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
# main_immunity_data <- merge(main_immunity_data,portal_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
# main_immunity_data <- merge(main_immunity_data,recb_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)



#Match PWM distance data. Contains PWM data for 264 phages.
main_immunity_data <- merge(main_immunity_data,stoperator_pwm_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)



#Compute comparison fields
main_immunity_data <- compute_comparisons(main_immunity_data)





#Retain only confident data, and discard questionable data
main_immunity_data_unreduced <- main_immunity_data

main_immunity_data <- subset(main_immunity_data,
                          main_immunity_data$lawn_reliability != 1 &
                            main_immunity_data$phage_reliability != 1)




# QC Summary - histogram of titers to assess the range of titers used
# Useful to compute multiplicity of infection
par(mar=c(4,8,8,4))
hist(log(main_immunity_data$tested_titer,10),xlim=c(0,10),
     main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",
     breaks=20)
dev.copy(pdf,"tested_titers.pdf")
dev.off()


###At this point, all data in main_immunity_data is derived from phages present in actino1321 database AND only confident data###












### Average data



#For certain analyses, I need averaged non-duplicated immunity comparisons. 
#Averages can be computed by unique assay_strain_defending_challenging identifier, 
#or by unique strain_defending_challenging identifier (which merges multiple-titer and single-titer data)


conf_to_average <- main_immunity_data
conf_to_average$assay_strain_defending_challenging <- factor(conf_to_average$assay_strain_defending_challenging)
conf_to_average$strain_defending_challenging <- factor(conf_to_average$strain_defending_challenging)
conf_to_average$rank6 <- as.numeric(as.character(conf_to_average$rank6))



#Average infection scores for each unique assay_strain_defending_challenging identifier

conf_assay_strain_def_chal_average <- aggregate(conf_to_average[,'rank6'],
                                                list(conf_to_average$assay_strain_defending_challenging),mean)


names(conf_assay_strain_def_chal_average) <- c('assay_strain_defending_challenging',
                                               'averaged_rank6') 

conf_assay_strain_def_chal_average$assay_strain_defending_challenging <- factor(conf_assay_strain_def_chal_average$assay_strain_defending_challenging)



#Compute the range of scores for each unique assay
#First compute the minimum score and maximum score
#Then compute the range

conf_assay_strain_def_chal_min <- aggregate(conf_to_average[,'rank6'],
                                            list(conf_to_average$assay_strain_defending_challenging),min)


names(conf_assay_strain_def_chal_min) <- c('assay_strain_defending_challenging',
                                               'min_rank6') 

conf_assay_strain_def_chal_min$assay_strain_defending_challenging <- factor(conf_assay_strain_def_chal_min$assay_strain_defending_challenging)


conf_assay_strain_def_chal_max <- aggregate(conf_to_average[,'rank6'],
                                            list(conf_to_average$assay_strain_defending_challenging),max)

names(conf_assay_strain_def_chal_max) <- c('assay_strain_defending_challenging',
                                           'max_rank6') 

conf_assay_strain_def_chal_max$assay_strain_defending_challenging <- factor(conf_assay_strain_def_chal_max$assay_strain_defending_challenging)


conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            conf_assay_strain_def_chal_min,
                                            by.x="assay_strain_defending_challenging",
                                            by.y="assay_strain_defending_challenging")

conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            conf_assay_strain_def_chal_max,
                                            by.x="assay_strain_defending_challenging",
                                            by.y="assay_strain_defending_challenging")

conf_assay_strain_def_chal_average$range_rank6 <- conf_assay_strain_def_chal_average$max_rank6 - 
  conf_assay_strain_def_chal_average$min_rank6

conf_assay_strain_def_chal_average$range_rank6 <- as.factor(conf_assay_strain_def_chal_average$range_rank6)

#Create table of immunity data metadata that will be added back to the averaged assay_strain_defending_challenging data
#Immunity data columns to keep after averaging
# "assay_strain_defending_challenging"
# "strain_defending_challenging"
# "defending_challenging"
# "prophage"
# "repressor_clone"
# "strain_type"
# "defending_phage"
# "challenging_phage"
# "assay_type"

reduced_immunity_metadata1 <- subset(conf_to_average,select=c('assay_strain_defending_challenging',
                                                                  'strain_defending_challenging',
                                                                  'defending_challenging',
                                                                  'prophage',
                                                                  'repressor_clone',
                                                                  'strain_type',
                                                                  'defending_phage',
                                                                  'challenging_phage',
                                                                  'assay_type'))


reduced_immunity_metadata1 <- reduced_immunity_metadata1[!duplicated(reduced_immunity_metadata1),]
reduced_immunity_metadata1$assay_strain_defending_challenging <- factor(reduced_immunity_metadata1$assay_strain_defending_challenging)
reduced_immunity_metadata1$strain_defending_challenging <- factor(reduced_immunity_metadata1$strain_defending_challenging)
reduced_immunity_metadata1$defending_challenging <- factor(reduced_immunity_metadata1$defending_challenging)
reduced_immunity_metadata1$prophage <- factor(reduced_immunity_metadata1$prophage)
reduced_immunity_metadata1$repressor_clone <- factor(reduced_immunity_metadata1$repressor_clone)
reduced_immunity_metadata1$strain_type <- factor(reduced_immunity_metadata1$strain_type)
reduced_immunity_metadata1$defending_phage <- factor(reduced_immunity_metadata1$defending_phage)
reduced_immunity_metadata1$challenging_phage <- factor(reduced_immunity_metadata1$challenging_phage)
reduced_immunity_metadata1$assay_type <- factor(reduced_immunity_metadata1$assay_type)



conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            reduced_immunity_metadata1,
                                            by.x='assay_strain_defending_challenging',
                                            by.y='assay_strain_defending_challenging')


#Match the phage metadata
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('defending','_',names(phage_metadata_to_match),sep="")

#There are 2 defending phages not in actino1321 = 'l5_gp71-flag' and 'redrock_bxb1'. So for right now keep all rows.
#conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,phage_metadata_to_match,by.x="defending_phage",by.y="defending_phageid",all.x=TRUE)

conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            phage_metadata_to_match,
                                            by.x="defending_phage",
                                            by.y="defending_phageid")




phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('challenging','_',names(phage_metadata_to_match),sep="")

#There are several challenging phages not in actino1321 = 'd29_mutant','l5_flag1-1','l5_ha1','l5_ha2','l5_ha3','l5_ha3_hazy'.
#So for right now keep all rows.
#conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,phage_metadata_to_match,by.x="challenging_phage",by.y="challenging_phageid",all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            phage_metadata_to_match,
                                            by.x="challenging_phage",
                                            by.y="challenging_phageid")



#Match the genomic distance data. Defending phages and challenging phages not in actino1319 will not be matched
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            genomic_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2")

#Match the gene distance data. Many comparisons will not be matched
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            repressor336_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            cas4311_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            endovii306_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            dnapol311_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            portal311_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            stoperator_pwm_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)


#Compute comparison fields
conf_assay_strain_def_chal_average <- compute_comparisons(conf_assay_strain_def_chal_average)



#Compute frequency of each unique assay_strain identifier and append to averaged datasets
conf_to_average_assay_strain_def_chal_freq <- as.data.frame(table(conf_to_average$assay_strain_defending_challenging))
names(conf_to_average_assay_strain_def_chal_freq) <- c('assay_strain_defending_challenging','frequency')
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            conf_to_average_assay_strain_def_chal_freq,
                                            by.x='assay_strain_defending_challenging',
                                            by.y='assay_strain_defending_challenging')

conf_assay_strain_def_chal_average$frequency <- as.factor(conf_assay_strain_def_chal_average$frequency) 






# TODO rename to 'Table S2'? and remove various fields?
#Export all averaged data
setwd("~/scratch/immunity_analysis/output/")
write.table(conf_assay_strain_def_chal_average,
            "conf_assay_strain_def_chal_average.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)






#Create average assay_strain_defending_challenging data above





















#Map averages back to main immunity table, subtract average from original values, then plot histogram of differences to show how reliable the dataset is

conf_assay_strain_def_chal_average_reduced <- subset(conf_assay_strain_def_chal_average,
                                                     select = c("assay_strain_defending_challenging",
                                                                "averaged_rank6",
                                                                "min_rank6",
                                                                "max_rank6",
                                                                "range_rank6",
                                                                "frequency"))


names(conf_assay_strain_def_chal_average_reduced) <- paste('assay_strain',
                                                           '_',
                                                           names(conf_assay_strain_def_chal_average_reduced),
                                                           sep="")







#main_immunity_data at this point retains only phages that are present in actino1321 and only confident data.
main_immunity_data <- merge(main_immunity_data,conf_assay_strain_def_chal_average_reduced,
                            by.x='assay_strain_defending_challenging',
                            by.y='assay_strain_assay_strain_defending_challenging',
                            all.x=TRUE)

main_immunity_data$assay_strain_averaged_rank6_diff <- as.numeric(as.character(main_immunity_data$rank6)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_rank6))
#Note: if the diff is positive = assay's rank is higher than the average








#Assess reproducibility for all comparisons with frequency > 1
main_immunity_data_assay_strain_reduced <- subset(main_immunity_data,
                                                  as.numeric(as.character(main_immunity_data$assay_strain_frequency)) > 1)

# main_immunity_data_strain_reduced <- subset(main_immunity_data,
#                                             main_immunity_data$strain_frequency > 1)



conf_assay_strain_def_chal_average$frequency <- as.numeric(as.character(conf_assay_strain_def_chal_average$frequency))
conf_assay_strain_def_chal_average$frequency2 <- ifelse(conf_assay_strain_def_chal_average$frequency > 10,11,conf_assay_strain_def_chal_average$frequency)
conf_assay_strain_def_chal_average$frequency <- as.factor(conf_assay_strain_def_chal_average$frequency)
conf_assay_strain_def_chal_average$frequency2 <- as.factor(conf_assay_strain_def_chal_average$frequency2)













#Plots
setwd("~/scratch/immunity_analysis/output/")

plot_bargraph2(conf_assay_strain_def_chal_average,
               "frequency",
               c(0,800),
               "conf_assay_strain_def_chal_ave_frequency.pdf")

plot_bargraph2(conf_assay_strain_def_chal_average,
               "frequency2",
               c(0,800),
               "conf_assay_strain_def_chal_ave_frequency_adjusted.pdf")



#Number of unique assays with >1 replicate. This includes all single-titer assays
1 - nrow(subset(conf_assay_strain_def_chal_average,conf_assay_strain_def_chal_average$frequency == "1"))/nrow(conf_assay_strain_def_chal_average)
#48% of unique comparisons has > 1 replicate 


plot_bargraph2(conf_assay_strain_def_chal_average,
               "range_rank6",
               c(0,1200),
               "conf_assay_strain_def_chal_ave_rank6_range.pdf")
#92% of unique comparisons with > 1 replicate have a score range of < 2.










#Only look at multiple_titer assays
conf_assay_strain_def_chal_average_multi <- subset(conf_assay_strain_def_chal_average,conf_assay_strain_def_chal_average$assay_type == "multiple_titer")



# Summary - number of unique multi-titer assays
nrow(conf_assay_strain_def_chal_average_multi)

#Summary - number of unique multi-titer assays with >1 replicate.
conf_assay_strain_def_chal_average_multi_reps <- subset(conf_assay_strain_def_chal_average_multi,conf_assay_strain_def_chal_average_multi$frequency != "1")
nrow(conf_assay_strain_def_chal_average_multi_reps)/nrow(conf_assay_strain_def_chal_average_multi)



#Summary - number of unique multi-titer assays with > 1 replicate and score range < 2
nrow(subset(conf_assay_strain_def_chal_average_multi_reps,conf_assay_strain_def_chal_average_multi_reps$range_rank6 == "0" | conf_assay_strain_def_chal_average_multi_reps$range_rank6 == "1"))/nrow(conf_assay_strain_def_chal_average_multi_reps)


###End of data average step
























###Compute immunity profile correlation coefficients 
setwd("~/scratch/immunity_analysis/input/")





#Import table of averaged data manipulated in Excel: reduced dataset, complete reciprocal matrix
infection_table_reduced <- read.csv("multi_lys_ave_infection_table_reduced.csv",
                                    sep=",",
                                    header=TRUE,
                                    row.names=1,
                                    na.strings="unspecified")


infection_table_reduced_t <- as.data.frame(t(infection_table_reduced))



#No data should be missing in this dataset
defending_cor_reduced <- cor(infection_table_reduced,method="pearson",use="complete.obs")
challenging_cor_reduced <- cor(infection_table_reduced_t,method="pearson",use="complete.obs")



#Convert to 3-column data frame
#Resulting table contains reciprocal data and self-comparisons
defending_cor_reduced_df <-melt(defending_cor_reduced)
challenging_cor_reduced_df <-melt(challenging_cor_reduced)


names(defending_cor_reduced_df) <- c("phage1","phage2","defending_cor_reduced")
names(challenging_cor_reduced_df) <- c("phage1","phage2","challenging_cor_reduced")

defending_cor_reduced_df$phage1_phage2 <- paste(defending_cor_reduced_df$phage1,
                                            "_",
                                            defending_cor_reduced_df$phage2,
                                            sep="")

challenging_cor_reduced_df$phage1_phage2 <- paste(challenging_cor_reduced_df$phage1,
                                              "_",
                                              challenging_cor_reduced_df$phage2,
                                              sep="")


defending_cor_reduced_df$phage1_phage2 <- as.factor(defending_cor_reduced_df$phage1_phage2)
challenging_cor_reduced_df$phage1_phage2 <- as.factor(challenging_cor_reduced_df$phage1_phage2)

defending_cor_reduced_df <- subset(defending_cor_reduced_df,select=c("phage1_phage2","defending_cor_reduced"))
challenging_cor_reduced_df <- subset(challenging_cor_reduced_df,select=c("phage1_phage2","challenging_cor_reduced"))


#






#Merge datasets




#If not using cor_all_df data
immunity_correlation_data <- merge(defending_cor_reduced_df,
                                   challenging_cor_reduced_df,
                                   by.x="phage1_phage2",
                                   by.y="phage1_phage2",
                                   all.x = TRUE,
                                   all.y = TRUE)

#This merged dataset contains reciprocal data and self-comparisons,
#and can now be merged with other tables elsewhere in the analysis.




###Compute immunity profile correlation coefficients above in progress




















###Plot immunity phenotypes by genome metrics to assess all data






#Below: averaged data
setwd("~/scratch/immunity_analysis/output/")


conf_assay_strain_ave_lys_multi_env_temp_rep <- subset(conf_assay_strain_def_chal_average,
                                                       conf_assay_strain_def_chal_average$strain_type == 'lysogen' &
                                                         conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                                         conf_assay_strain_def_chal_average$source_compare == 'environment' &
                                                         conf_assay_strain_def_chal_average$temperate_empirical_compare == 'yes' &
                                                         conf_assay_strain_def_chal_average$functional_repressor_compare == 'yes')




conf_assay_strain_ave_lys_multi_env_temp_rep$defending_challenging <- factor(conf_assay_strain_ave_lys_multi_env_temp_rep$defending_challenging)
conf_assay_strain_ave_lys_multi_env_temp_rep$defending_phage <- factor(conf_assay_strain_ave_lys_multi_env_temp_rep$defending_phage)
conf_assay_strain_ave_lys_multi_env_temp_rep$challenging_phage <- factor(conf_assay_strain_ave_lys_multi_env_temp_rep$challenging_phage)



conf_assay_strain_ave_lys_multi_TempDiff_RepDiff <- subset(conf_assay_strain_def_chal_average,
                                                           conf_assay_strain_def_chal_average$strain_type == 'lysogen' &
                                                             conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                                             conf_assay_strain_def_chal_average$defending_source == 'environment' &
                                                             (conf_assay_strain_def_chal_average$temperate_empirical_compare == 'different' |
                                                                conf_assay_strain_def_chal_average$functional_repressor_compare == 'different'))



conf_assay_strain_ave_lys_multi_env_temp_rep_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$gene_content_clade_compare == "different")

conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                                   conf_assay_strain_ave_lys_multi_env_temp_rep$gene_content_clade_compare == "clade2")


conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2,
                                                                             as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2$defending_phage) == 
                                                                               as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2$challenging_phage))

conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2,
                                                                               as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2$defending_phage) != 
                                                                                 as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2$challenging_phage))





#split by clade and homo/heterotypic









#Fig. 5b sub-panel 1
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repFull.pdf")


#Fig. 5b sub-panel 2
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_stopEuc.pdf")


#Fig. S5a sub-panel 1
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "pham_pham_dissimilarity",
                       "averaged_rank6",
                       c(0,1),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_gcd.pdf")


#Fig. S5a sub-panel 2
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "modified_mash_distance",
                       "averaged_rank6",
                       c(0,0.5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_mash.pdf")


#Fig. S5d sub-panel 1
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "portal_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_portal.pdf")


#Fig. S5d sub-panel 2
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_dnapol.pdf")


#Fig. S5d sub-panel 3
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_endovii.pdf")


#Fig. S5d sub-panel 4
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_cas4.pdf")


#Fig. 9b sub-panel 1
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repNterm.pdf")


#Fig. 9b sub-panel 2
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "repressor_hth_compare",
                       "averaged_rank6",
                       c(0,10),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repHTH.pdf")


#Fig. 9b sub-panel 3
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_rank6",
                       c(0,70),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repCterm.pdf")








#Compare differences between prophage maintenance strategy
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                              conf_assay_strain_ave_lys_multi_env_temp_rep$lysogen_type_compare == 'integration')

conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$lysogen_type_compare == 'extrachromosomal')

conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                                conf_assay_strain_ave_lys_multi_env_temp_rep$lysogen_type_compare == 'different' &
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$defending_lysogen_type == 'integration' &
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$challenging_lysogen_type == 'extrachromosomal')

conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                                conf_assay_strain_ave_lys_multi_env_temp_rep$lysogen_type_compare == 'different' &
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$defending_lysogen_type == 'extrachromosomal' &
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep$challenging_lysogen_type == 'integration')


conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$defending_pham_integrase <- factor(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$defending_pham_integrase)
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$defending_pham_parb <- factor(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$defending_pham_parb)

conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt,
                                                                          conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$integrase_compare != "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt,
                                                                          conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$integrase_compare == "different")

conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra,
                                                                           conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$parb_compare != "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra,
                                                                           conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$parb_compare == "different")







#IntExtra
conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra,
                                                                           conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra,
                                                                            conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$gene_content_clade_compare == "clade2")

conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2,
                                                                                                 as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$defending_phage) ==
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2,
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$defending_phage) !=
                                                                                                     as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_heterotypic)
#








#ExtraInt
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt,
                                                                           conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt,
                                                                            conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$gene_content_clade_compare == "clade2")

conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2,
                                                                                                 as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$defending_phage) ==
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2,
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$defending_phage) !=
                                                                                                     as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_heterotypic)
#






#IntInt Same Pham
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame,
                                                                                     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame,
                                                                                      conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$gene_content_clade_compare == "clade2")
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2,
                                                                                                as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2$defending_phage) ==
                                                                                                  as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2,
                                                                                                as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2$defending_phage) !=
                                                                                                  as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic)
#




#IntInt Different Pham
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff,
                                                                                     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff,
                                                                                      conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$gene_content_clade_compare == "clade2")

conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2,
                                                                                                 as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$defending_phage) ==
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2,
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$defending_phage) !=
                                                                                                     as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_heterotypic)

#






#ExtraExtra Same ParB
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame,
                                                                                      conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame,
                                                                                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$gene_content_clade_compare == "clade2")
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2,
                                                                                                as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2$defending_phage) ==
                                                                                                  as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2,
                                                                                                  as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2$defending_phage) !=
                                                                                                    as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic)
#








#ExtraExtra Different ParB
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff,
                                                                                      conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$gene_content_clade_compare == "different")
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff,
                                                                                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$gene_content_clade_compare == "clade2")

conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_homotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2,
                                                                                                 as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$defending_phage) ==
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$challenging_phage))
conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_heterotypic <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2,
                                                                                                   as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$defending_phage) !=
                                                                                                     as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$challenging_phage))

#QC
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_heterotypic)
#




#Fig. S5b sub-panel 1 - IntInt_IntPhamSame
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntSame_stopEuc_subtypes.pdf")


#Fig. S5b sub-panel 2 - IntInt_IntPhamDiff
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntDiff_stopEuc_subtypes.pdf")


#Fig. S5b sub-panel 3 - IntExtra
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_stopEuc_subtypes.pdf")


#Fig. S5c sub-panel 1 - ExtraExtra_ParBPhamSame
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraSame_stopEuc_subtypes.pdf")


#Fig. S5c sub-panel 2 - ExtraExtra_ParBPhamDiff
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraDiff_stopEuc_subtypes.pdf")


#Fig. S5c sub-panel 3 - ExtraInt
plot_tricolor_scatter2(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_heterotypic,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade,
                       conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_stopEuc_subtypes.pdf")




#Binned Data Stats

#Plot only Clade 2 OR (Clade 2 non-homotypic) immunity data stats

conf_assay_strain_ave_lys_multi <- subset(conf_assay_strain_def_chal_average,
                                          conf_assay_strain_def_chal_average$strain_type == 'lysogen' &
                                            conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')


conf_assay_strain_ave_lys_multi$defending_lysogen_type <- factor(conf_assay_strain_ave_lys_multi$defending_lysogen_type)
conf_assay_strain_ave_lys_multi$challenging_lysogen_type <- factor(conf_assay_strain_ave_lys_multi$challenging_lysogen_type)

conf_assay_strain_ave_lys_multi$challenging_phage <- factor(conf_assay_strain_ave_lys_multi$challenging_phage)
conf_assay_strain_ave_lys_multi$defending_phage <- factor(conf_assay_strain_ave_lys_multi$defending_phage)


clade2_env <- subset(conf_assay_strain_ave_lys_multi,
                     conf_assay_strain_ave_lys_multi$gene_content_clade_compare == "clade2" &
                     conf_assay_strain_ave_lys_multi$source_compare == "environment")

clade2_env_temp <- clade2_env
clade2_env <- subset(clade2_env,as.character(clade2_env$defending_phage) != as.character(clade2_env$challenging_phage))



clade2_env$defending_phage <- factor(clade2_env$defending_phage)
clade2_env$challenging_phage <- factor(clade2_env$challenging_phage)




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



clade2_bin0_freq <- data.frame("bin0",
                        nlevels(clade2_env_bin0$defending_phage),
                        nlevels(clade2_env_bin0$challenging_phage),
                        nrow(clade2_env_bin0),
                        nrow(subset(clade2_env_bin0,
                                    clade2_env_bin0$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin0,
                                    clade2_env_bin0$subcluster_compare == "different")))

clade2_bin1_freq <- data.frame("bin1",
                        nlevels(clade2_env_bin1$defending_phage),
                        nlevels(clade2_env_bin1$challenging_phage),
                        nrow(clade2_env_bin1),
                        nrow(subset(clade2_env_bin1,
                                    clade2_env_bin1$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin1,
                                    clade2_env_bin1$subcluster_compare == "different")))

clade2_bin2_freq <- data.frame("bin2",
                        nlevels(clade2_env_bin2$defending_phage),
                        nlevels(clade2_env_bin2$challenging_phage),
                        nrow(clade2_env_bin2),
                        nrow(subset(clade2_env_bin2,
                                    clade2_env_bin2$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin2,
                                    clade2_env_bin2$subcluster_compare == "different")))

clade2_bin3_freq <- data.frame("bin3",
                        nlevels(clade2_env_bin3$defending_phage),
                        nlevels(clade2_env_bin3$challenging_phage),
                        nrow(clade2_env_bin3),
                        nrow(subset(clade2_env_bin3,
                                    clade2_env_bin3$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin3,
                                    clade2_env_bin3$subcluster_compare == "different")))

clade2_bin4_freq <- data.frame("bin4",
                        nlevels(clade2_env_bin4$defending_phage),
                        nlevels(clade2_env_bin4$challenging_phage),
                        nrow(clade2_env_bin4),
                        nrow(subset(clade2_env_bin4,
                                    clade2_env_bin4$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin4,
                                    clade2_env_bin4$subcluster_compare == "different")))

clade2_bin5_freq <- data.frame("bin5",
                        nlevels(clade2_env_bin5$defending_phage),
                        nlevels(clade2_env_bin5$challenging_phage),
                        nrow(clade2_env_bin5),
                        nrow(subset(clade2_env_bin5,
                                    clade2_env_bin5$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin5,
                                    clade2_env_bin5$subcluster_compare == "different")))

clade2_bin6_freq <- data.frame("bin6",
                        nlevels(clade2_env_bin6$defending_phage),
                        nlevels(clade2_env_bin6$challenging_phage),
                        nrow(clade2_env_bin6),
                        nrow(subset(clade2_env_bin6,
                                    clade2_env_bin6$subcluster_compare != "different")),
                        nrow(subset(clade2_env_bin6,
                                    clade2_env_bin6$subcluster_compare == "different")))

names(clade2_bin0_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin1_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin2_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin3_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin4_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin5_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")
names(clade2_bin6_freq) <- c("bin","defending_freq","challenging_freq","total_num_assays","num_intra_subcluster_assays","num_inter_subcluster_assays")

clade2_bin0_freq$defending_percent <- clade2_bin0_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin1_freq$defending_percent <- clade2_bin1_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin2_freq$defending_percent <- clade2_bin2_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin3_freq$defending_percent <- clade2_bin3_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin4_freq$defending_percent <- clade2_bin4_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin5_freq$defending_percent <- clade2_bin5_freq$defending_freq/nlevels(clade2_env$defending_phage)
clade2_bin6_freq$defending_percent <- clade2_bin6_freq$defending_freq/nlevels(clade2_env$defending_phage)

clade2_bin0_freq$challenging_percent <- clade2_bin0_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin1_freq$challenging_percent <- clade2_bin1_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin2_freq$challenging_percent <- clade2_bin2_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin3_freq$challenging_percent <- clade2_bin3_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin4_freq$challenging_percent <- clade2_bin4_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin5_freq$challenging_percent <- clade2_bin5_freq$challenging_freq/nlevels(clade2_env$challenging_phage)
clade2_bin6_freq$challenging_percent <- clade2_bin6_freq$challenging_freq/nlevels(clade2_env$challenging_phage)

clade2_bin0_freq$total_assays_percent <- clade2_bin0_freq$total_num_assays/nrow(clade2_env)
clade2_bin1_freq$total_assays_percent <- clade2_bin1_freq$total_num_assays/nrow(clade2_env)
clade2_bin2_freq$total_assays_percent <- clade2_bin2_freq$total_num_assays/nrow(clade2_env)
clade2_bin3_freq$total_assays_percent <- clade2_bin3_freq$total_num_assays/nrow(clade2_env)
clade2_bin4_freq$total_assays_percent <- clade2_bin4_freq$total_num_assays/nrow(clade2_env)
clade2_bin5_freq$total_assays_percent <- clade2_bin5_freq$total_num_assays/nrow(clade2_env)
clade2_bin6_freq$total_assays_percent <- clade2_bin6_freq$total_num_assays/nrow(clade2_env)

clade2_bin0_freq$intra_subcluster_percent <- clade2_bin0_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin1_freq$intra_subcluster_percent <- clade2_bin1_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin2_freq$intra_subcluster_percent <- clade2_bin2_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin3_freq$intra_subcluster_percent <- clade2_bin3_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin4_freq$intra_subcluster_percent <- clade2_bin4_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin5_freq$intra_subcluster_percent <- clade2_bin5_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))
clade2_bin6_freq$intra_subcluster_percent <- clade2_bin6_freq$num_intra_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare != "different"))

clade2_bin0_freq$inter_subcluster_percent <- clade2_bin0_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin1_freq$inter_subcluster_percent <- clade2_bin1_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin2_freq$inter_subcluster_percent <- clade2_bin2_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin3_freq$inter_subcluster_percent <- clade2_bin3_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin4_freq$inter_subcluster_percent <- clade2_bin4_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin5_freq$inter_subcluster_percent <- clade2_bin5_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))
clade2_bin6_freq$inter_subcluster_percent <- clade2_bin6_freq$num_inter_subcluster_assays/nrow(subset(clade2_env,clade2_env$subcluster_compare == "different"))






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







#Fig. S4a sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "defending_percent",
               "bin",
               c(0,1),
               "conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_defending.pdf")


#Fig. S4a sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "challenging_percent",
               "bin",
               c(0,1),
               "conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_challenging.pdf")


#Fig. S4a sub-panel 3
plot_bargraph1(clade2_binned_frequency,
               "total_assays_percent",
               "bin",
               c(0,0.3),
               "conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_total.pdf")


#Fig. S4b sub-panel 1
plot_bargraph1(clade2_binned_frequency,
               "inter_subcluster_percent",
               "bin",
               c(0,0.5),
               "conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_intersubcluster.pdf")


#Fig. S4b sub-panel 2
plot_bargraph1(clade2_binned_frequency,
               "intra_subcluster_percent",
               "bin",
               c(0,0.5),
               "conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_intrasubcluster.pdf")



###Above: averaged data

























###Reciprocal analysis
#Using averaged data, compute differences in reciprocal immunity data
#Split conf_assay_strain_def_chal_average columns into groups


#Phage-specific data = metadata that is not impacted by immunity vector
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
# "defending_cluster_a_functional_repressor_predicted"
# "defending_cluster_a_temperate_empirical"
# "defending_repressor_hth_domain_sequence"
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
# "challenging_cluster_a_functional_repressor_predicted"
# "challenging_cluster_a_temperate_empirical"
# "challenging_repressor_hth_domain_sequence"
# "challenging_repressor_length_full"
# "challenging_repressor_length_nterm"
# "challenging_repressor_length_cterm"
# "challenging_pham_parb"


#Phage metadata comparisons
# data specific to both phages used in immunity but not impacted by vector orientation
# "modified_mash_distance"
# "pham_pham_dissimilarity"
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


#Vectored_data = data that is specific to immunity assay and depends on vector orientation
# "assay_strain_defending_challenging"
# "defending_challenging"
# "challenging_phage"
# "defending_phage"
# "averaged_infection_strength"
# "averaged_turbidity"
# "averaged_plaque_size"
# "averaged_plaques"
# "averaged_four_factors"
# "averaged_rank6"
# "strain_defending_challenging"
# "assay_type"
# "frequency"

vectored_column_names <- c("assay_strain_defending_challenging",
                           "defending_challenging",
                           "defending_phage",
                           "challenging_phage",
                           "averaged_rank6",
                           "strain_defending_challenging",
                           "assay_type",
                           "frequency"
)

reciprocal_data <- subset(conf_assay_strain_def_chal_average,
                          conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                            conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')


vector2_data <- subset(reciprocal_data,select = vectored_column_names)
names(vector2_data) <- paste('vector2_',names(vector2_data),sep="")




for (column_name in vectored_column_names){
  
  new_column_name <- paste('vector1_',column_name,sep="")
  
  print(new_column_name)
  
  names(reciprocal_data)[names(reciprocal_data) == column_name] <- paste('vector1_',column_name,sep="")
  
}






reciprocal_data$vector2_defending_challenging <- paste(reciprocal_data$vector1_challenging_phage,"_",reciprocal_data$vector1_defending_phage,sep="")
reciprocal_data <- merge(reciprocal_data,vector2_data,by.x='vector2_defending_challenging',by.y='vector2_defending_challenging')


#Compute difference in infection profiles
reciprocal_data$averaged_rank6_diff <- abs(reciprocal_data$vector1_averaged_rank6 - reciprocal_data$vector2_averaged_rank6)


#Now that all data has been matched, the data is duplicated. (vector1 will contain l5_alma and alma_l5, so both vector2's will be present as well)
#To solve this, remove all vector1 rows in which the defending_phage and challenging_phage are not alphabetically ordered.
#This does NOT retain self-comparisons (e.g. alma_alma)
reciprocal_data$vector1_alpha_ordered <- as.character(reciprocal_data$vector1_defending_phage) < as.character(reciprocal_data$vector1_challenging_phage)



#Now retain only the unique pairwise comparisons
reciprocal_data_alpha_ordered <- subset(reciprocal_data,reciprocal_data$vector1_alpha_ordered == TRUE)



#Output the reciprocal dataset
setwd("~/scratch/immunity_analysis/output/")
write.table(reciprocal_data_alpha_ordered,
            "reciprocal_data_alpha_ordered.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)





#Reduce the reciprocal dataset to only environmental phages 
reciprocal_unique_envY_intY <- subset(reciprocal_data,
                                      reciprocal_data$vector1_alpha_ordered == TRUE &
                                        reciprocal_data$source_compare == 'environment' &
                                        reciprocal_data$lysogen_type_compare == 'integration')

reciprocal_unique_envY_extraY <- subset(reciprocal_data,
                                        reciprocal_data$vector1_alpha_ordered == TRUE &
                                          reciprocal_data$source_compare == 'environment' &
                                          reciprocal_data$lysogen_type_compare == 'extrachromosomal')

reciprocal_unique_envY <- subset(reciprocal_data,
                                 reciprocal_data$vector1_alpha_ordered == TRUE &
                                   reciprocal_data$source_compare == 'environment')



#Binning

#Bin reciprocal data for histogram

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




reciprocal_bin0_freq <- data.frame("bin0",nrow(reciprocal_bin0))
reciprocal_bin1_freq <- data.frame("bin1",nrow(reciprocal_bin1))
reciprocal_bin2_freq <- data.frame("bin2",nrow(reciprocal_bin2))
reciprocal_bin3_freq <- data.frame("bin3",nrow(reciprocal_bin3))
reciprocal_bin4_freq <- data.frame("bin4",nrow(reciprocal_bin4))
reciprocal_bin5_freq <- data.frame("bin5",nrow(reciprocal_bin5))
reciprocal_bin6_freq <- data.frame("bin6",nrow(reciprocal_bin6))


names(reciprocal_bin0_freq) <- c("bin","freq")
names(reciprocal_bin1_freq) <- c("bin","freq")
names(reciprocal_bin2_freq) <- c("bin","freq")
names(reciprocal_bin3_freq) <- c("bin","freq")
names(reciprocal_bin4_freq) <- c("bin","freq")
names(reciprocal_bin5_freq) <- c("bin","freq")
names(reciprocal_bin6_freq) <- c("bin","freq")


reciprocal_bin0_freq$freq_percent <- reciprocal_bin0_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin1_freq$freq_percent <- reciprocal_bin1_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin2_freq$freq_percent <- reciprocal_bin2_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin3_freq$freq_percent <- reciprocal_bin3_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin4_freq$freq_percent <- reciprocal_bin4_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin5_freq$freq_percent <- reciprocal_bin5_freq$freq/nrow(reciprocal_unique_envY)
reciprocal_bin6_freq$freq_percent <- reciprocal_bin6_freq$freq/nrow(reciprocal_unique_envY)


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


#QC
sum(reciprocal_binned_freq$freq) - nrow(reciprocal_unique_envY)






#Plot data
setwd("~/scratch/immunity_analysis/output/")




#Fig. S4c
plot_bargraph1(reciprocal_binned_freq,
               "freq_percent",
               "bin",
               c(0,0.5),
               "reciprocal_unique_env_binned_freq_percent.pdf")










#How do various genome/gene distances correlate with differences in reciprocal profiles?



#rank6

#reciprocal_unique_envY: subsetted from conf, assay_strain_ave, lys, multi, env data. Since it is reciprocal, it is by definition tempY and repY data.
reciprocal_unique_envY_interclade <- subset(reciprocal_unique_envY,
                                            reciprocal_unique_envY$gene_content_clade_compare == "different")

reciprocal_unique_envY_intraclade2 <- subset(reciprocal_unique_envY,
                                             reciprocal_unique_envY$gene_content_clade_compare == "clade2")

reciprocal_unique_envY_intraclade2_homotypic <- subset(reciprocal_unique_envY_intraclade2,
                                                         as.character(reciprocal_unique_envY_intraclade2$vector1_defending_phage) == 
                                                           as.character(reciprocal_unique_envY_intraclade2$vector1_challenging_phage))


reciprocal_unique_envY_intraclade2_heterotypic <- subset(reciprocal_unique_envY_intraclade2,
                                                         as.character(reciprocal_unique_envY_intraclade2$vector1_defending_phage) != 
                                                           as.character(reciprocal_unique_envY_intraclade2$vector1_challenging_phage))



#Fig. 5b sub-panel 3
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "repressor_full_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_repFull.pdf")


#Fig. 5b sub-panel 4
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "stoperator_pwd_dist_euc",
                       "averaged_rank6_diff",
                       c(0,5),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_stopEuc.pdf")


#Fig. S5a sub-panel 3
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "pham_pham_dissimilarity",
                       "averaged_rank6_diff",
                       c(0,1),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_gcd.pdf")


#Fig. S5a sub-panel 4
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "modified_mash_distance",
                       "averaged_rank6_diff",
                       c(0,0.5),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_mash.pdf")


#Fig. S5d sub-panel 5
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "portal_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_portal.pdf")


#Fig. S5d sub-panel 6
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "dnapol_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_dnapol.pdf")


#Fig. S5d sub-panel 7
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "endovii_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_endovii.pdf")


#Fig. S5d sub-panel 8
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "cas4_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_cas4.pdf")


#Fig. 9b sub-panel 4
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "repressor_nterm_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_repNterm.pdf")


#Fig. 9b sub-panel 5
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "repressor_hth_compare",
                       "averaged_rank6_diff",
                       c(0,10),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_repHTH.pdf")


#Fig. 9b sub-panel 6
plot_tricolor_scatter2(reciprocal_unique_envY_intraclade2_heterotypic,
                       reciprocal_unique_envY_interclade,
                       reciprocal_unique_envY_intraclade2_homotypic,
                       "repressor_cterm_mafft_dist_uncorrected",
                       "averaged_rank6_diff",
                       c(0,70),
                       c(0,6),
                       "reciprocal_unique_envY_subtypes_by_repCterm.pdf")



###Reciprocal analysis above
















###Compare lysogen vs repressor clone data




#Lysogen-Repressor analyses using averaged data
clone_match_columns <- c('defending_challenging',
                         'defending_phage',
                         'challenging_phage',
                         'averaged_rank6',
                         'modified_mash_distance',
                         'pham_pham_dissimilarity',
                         'repressor_cterm_mafft_dist_uncorrected',
                         'stoperator_pwd_dist_euc',
                         'gene_content_clade_compare')


#All averaged multiple-titer lysogen and clone data
ave_multiY_confY_lysY <- subset(conf_assay_strain_def_chal_average,
                                  conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                  conf_assay_strain_def_chal_average$strain_type == 'lysogen')

ave_multiY_confY_cloneY <- subset(conf_assay_strain_def_chal_average,
                                    conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                    conf_assay_strain_def_chal_average$strain_type == 'repressor_clone')

ave_multi_conf_lys_clone_matched <- match_average_lysogen_and_clone_data(ave_multiY_confY_lysY,
                                                                         ave_multiY_confY_cloneY)



#Looks like I lose some comparisons when i subset the matched data into different types of phages.
nrow(ave_multiY_confY_lysY)
nrow(ave_multiY_confY_cloneY)
nrow(ave_multi_conf_lys_clone_matched)
lost_clones1 <- merge(ave_multiY_confY_cloneY,ave_multi_conf_lys_clone_matched,by.x="defending_challenging",by.y="lys_defending_challenging",all.x=TRUE)
lost_clones2 <- subset(lost_clones1,is.na(lost_clones1$clone_challenging_phage))
lost_clones2$defending_challenging <- factor(lost_clones2$defending_challenging)
#There are 11 repressor clone assays with no matching lysogen. 9 involve the Bxb1 repressor. One involves DaVinci-phiTM1, and one involves Trixie-phiTM46.





#Environment, func. repressor, and empirically temperate data from averaged multiple-titer lysogen and clone data
ave_multiY_confY_envY_repY_empY_lysY <- subset(conf_assay_strain_def_chal_average,
                                                 conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                                 conf_assay_strain_def_chal_average$source_compare == 'environment' &
                                                 conf_assay_strain_def_chal_average$functional_repressor_compare == 'yes' &
                                                 conf_assay_strain_def_chal_average$temperate_empirical_compare == 'yes' &
                                                 conf_assay_strain_def_chal_average$strain_type == 'lysogen')

ave_multiY_confY_envY_repY_empY_cloneY <- subset(conf_assay_strain_def_chal_average,
                                                   conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                                   conf_assay_strain_def_chal_average$source_compare == 'environment' &
                                                   conf_assay_strain_def_chal_average$functional_repressor_compare == 'yes' &
                                                   conf_assay_strain_def_chal_average$temperate_empirical_compare == 'yes' &
                                                   conf_assay_strain_def_chal_average$strain_type == 'repressor_clone')

ave_multi_conf_env_rep_emp_lys_clone_matched <- match_average_lysogen_and_clone_data(ave_multiY_confY_envY_repY_empY_lysY,
                                                                                     ave_multiY_confY_envY_repY_empY_cloneY)







#
ave_multi_conf_env_rep_emp_lys_clone_matched_interclade <- subset(ave_multi_conf_env_rep_emp_lys_clone_matched,
                                                                  ave_multi_conf_env_rep_emp_lys_clone_matched$lys_gene_content_clade_compare == "different")

ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2 <- subset(ave_multi_conf_env_rep_emp_lys_clone_matched,
                                                                   ave_multi_conf_env_rep_emp_lys_clone_matched$lys_gene_content_clade_compare == "clade2")

ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic <- subset(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2,
                                                                               as.character(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2$lys_defending_phage) == 
                                                                                 as.character(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2$lys_challenging_phage))

ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic <- subset(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2,
                                                                               as.character(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2$lys_defending_phage) != 
                                                                                 as.character(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2$lys_challenging_phage))




#









#Escape mutant data from averaged multiple-titer lysogen and clone data
ave_multiY_confY_envN_lysY <- subset(conf_assay_strain_def_chal_average,
                                     conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                       conf_assay_strain_def_chal_average$source_compare != 'environment' &
                                       conf_assay_strain_def_chal_average$strain_type == 'lysogen')

ave_multiY_confY_envN_cloneY <- subset(conf_assay_strain_def_chal_average,
                                       conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                         conf_assay_strain_def_chal_average$source_compare != 'environment' &
                                         conf_assay_strain_def_chal_average$strain_type == 'repressor_clone')

ave_multi_conf_envN_lys_clone_matched <- match_average_lysogen_and_clone_data(ave_multiY_confY_envN_lysY,
                                                                              ave_multiY_confY_envN_cloneY)



ave_multi_conf_envN_lys_clone_matched_interclade <- subset(ave_multi_conf_envN_lys_clone_matched,
                                                           ave_multi_conf_envN_lys_clone_matched$lys_gene_content_clade_compare == "different")

ave_multi_conf_envN_lys_clone_matched_intraclade2 <- subset(ave_multi_conf_envN_lys_clone_matched,
                                                            ave_multi_conf_envN_lys_clone_matched$lys_gene_content_clade_compare == "clade2")

ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic <- subset(ave_multi_conf_envN_lys_clone_matched_intraclade2,
                                                                      as.character(ave_multi_conf_envN_lys_clone_matched_intraclade2$lys_defending_phage) == 
                                                                        as.character(ave_multi_conf_envN_lys_clone_matched_intraclade2$lys_challenging_phage))

ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic <- subset(ave_multi_conf_envN_lys_clone_matched_intraclade2,
                                                                        as.character(ave_multi_conf_envN_lys_clone_matched_intraclade2$lys_defending_phage) != 
                                                                          as.character(ave_multi_conf_envN_lys_clone_matched_intraclade2$lys_challenging_phage))


#




#QC
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched)
nrow(ave_multi_conf_envN_lys_clone_matched)

nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade)

nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic) +
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic) + 
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade)


nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_envN_lys_clone_matched_interclade)
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)

nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic) + 
  nrow(ave_multi_conf_envN_lys_clone_matched_interclade) + 
  nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)

#





#Plots
setwd("~/scratch/immunity_analysis/output/")


#Fig. 6c
plot_tricolor_scatter1(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic,
                       ave_multi_conf_env_rep_emp_lys_clone_matched_interclade,
                       ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic,
                       "lys_averaged_rank6","clone_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "ave_multi_conf_env_rep_emp_lys_clone_matched_rank6compare_subtypes.pdf")


#Fig. 7d
plot_tricolor_scatter1(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic,
                       ave_multi_conf_envN_lys_clone_matched_interclade,
                       ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic,
                       "lys_averaged_rank6","clone_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "ave_multi_conf_envN_lys_clone_matched_rank6compare_subtypes.pdf")


#Fig. 6d
plot_tricolor_scatter3(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic,
                       ave_multi_conf_env_rep_emp_lys_clone_matched_interclade,
                       ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic,
                       "lys_stoperator_pwd_dist_euc",
                       "rank6_diff",
                       c(0,5),
                       c(-4,4),
                       "ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_stopEuc_subtypes.pdf")


#Fig. 7e
plot_tricolor_scatter3(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic,
                       ave_multi_conf_envN_lys_clone_matched_interclade,
                       ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic,
                       "lys_stoperator_pwd_dist_euc",
                       "rank6_diff",
                       c(0,5),
                       c(-4,4),
                       "ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_stopEuc_subtypes.pdf")


###Lysogen-clone comparison above














###L5-mutant comparisons using all averaged data



#All averaged multiple-titer lysogen and clone data
l5_columns <- c('defending_challenging',
                'defending_phage',
                'challenging_phage',
                'averaged_rank6',
                'modified_mash_distance',
                'pham_pham_dissimilarity',
                'repressor_cterm_mafft_dist_uncorrected')


ave_multiY_confY_lysY <- subset(conf_assay_strain_def_chal_average,
                                conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                  conf_assay_strain_def_chal_average$strain_type == 'lysogen')


#challenging correlation
chal_l5 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$challenging_phage == 'l5')


chal_phitm41 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$challenging_phage == 'phitm41',select = c(l5_columns))
chal_phitm1 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$challenging_phage == 'phitm1',select = c(l5_columns))
chal_phitm4 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$challenging_phage == 'phitm4',select = c(l5_columns))
chal_phitm6 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$challenging_phage == 'phitm6',select = c(l5_columns))


names(chal_l5) <- paste('l5_',names(chal_l5),sep="")
names(chal_phitm1) <- paste('phitm1_',names(chal_phitm1),sep="")
names(chal_phitm41) <- paste('phitm41_',names(chal_phitm41),sep="")
names(chal_phitm4) <- paste('phitm4_',names(chal_phitm4),sep="")
names(chal_phitm6) <- paste('phitm6_',names(chal_phitm6),sep="")



chal_l5_assays <- merge(chal_l5,chal_phitm41,by.x='l5_defending_phage',by.y='phitm41_defending_phage',all.x=TRUE,all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,chal_phitm1,by.x='l5_defending_phage',by.y='phitm1_defending_phage',all.x=TRUE,all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,chal_phitm4,by.x='l5_defending_phage',by.y='phitm4_defending_phage',all.x=TRUE,all.y=TRUE)
chal_l5_assays <- merge(chal_l5_assays,chal_phitm6,by.x='l5_defending_phage',by.y='phitm6_defending_phage',all.x=TRUE,all.y=TRUE)


chal_l5_assays$l5_phitm41_rank6_diff <- chal_l5_assays$phitm41_averaged_rank6 - chal_l5_assays$l5_averaged_rank6
chal_l5_assays$l5_phitm1_rank6_diff <- chal_l5_assays$phitm1_averaged_rank6 - chal_l5_assays$l5_averaged_rank6
chal_l5_assays$l5_phitm4_rank6_diff <- chal_l5_assays$phitm4_averaged_rank6 - chal_l5_assays$l5_averaged_rank6
chal_l5_assays$l5_phitm6_rank6_diff <- chal_l5_assays$phitm6_averaged_rank6 - chal_l5_assays$l5_averaged_rank6

chal_l5_assays$phitm1_phitm6_rank6_diff <- chal_l5_assays$phitm6_averaged_rank6 - chal_l5_assays$phitm1_averaged_rank6
chal_l5_assays$phitm1_phitm4_rank6_diff <- chal_l5_assays$phitm4_averaged_rank6 - chal_l5_assays$phitm1_averaged_rank6



chal_l5_assays_nonclade2 <- subset(chal_l5_assays,chal_l5_assays$l5_defending_gene_content_clade != "clade2")
chal_l5_assays_clade2 <- subset(chal_l5_assays,chal_l5_assays$l5_defending_gene_content_clade == "clade2")





#Right now I only need to compute homotypic comparisons for phiTM1
chal_l5_assays_clade2_phitm1homotypic <- subset(chal_l5_assays_clade2,
                                                as.character(chal_l5_assays_clade2$phitm1_challenging_phage) == 
                                                  as.character(chal_l5_assays_clade2$l5_defending_phage))
chal_l5_assays_clade2_phitm1heterotypic <- subset(chal_l5_assays_clade2,
                                                as.character(chal_l5_assays_clade2$phitm1_challenging_phage) != 
                                                  as.character(chal_l5_assays_clade2$l5_defending_phage))

#Used as a dummy table for plotting. phiTM4 is lytic, so unable to form a lysogen.
chal_l5_assays_clade2_empty <- subset(chal_l5_assays_clade2,chal_l5_assays_clade2$l5_defending_phage == 'phiTM4')

#QC
nrow(chal_l5_assays_clade2)
nrow(subset(chal_l5_assays_clade2,is.na(chal_l5_assays_clade2$phitm1_challenging_phage)))
nrow(chal_l5_assays_clade2_phitm1homotypic)
nrow(chal_l5_assays_clade2_phitm1heterotypic)
nrow(chal_l5_assays_clade2_empty)
#





#Plots
setwd("~/scratch/immunity_analysis/output/")


#Fig. 10d sub-panel 1
plot_tricolor_scatter2(chal_l5_assays_clade2_phitm1heterotypic,
                       chal_l5_assays_nonclade2,
                       chal_l5_assays_clade2_phitm1homotypic,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm1_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "chal_phitm1_rank6_vs_stopEuc_by_subtype.pdf")


#Fig. 10d sub-panel 2
plot_tricolor_scatter2(chal_l5_assays_clade2,
                       chal_l5_assays_nonclade2,
                       chal_l5_assays_clade2_empty,
                       "l5_stoperator_pwd_dist_euc",
                       "phitm4_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "chal_phitm4_rank6_vs_stopEuc_by_subtype.pdf")


#Fig. S8g sub-panel 1
plot_tricolor_scatter1(chal_l5_assays_clade2,
                       chal_l5_assays_nonclade2,
                       chal_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm41_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "chal_l5_phitm41_rank6_by_subtype.pdf")


#Fig. S8g sub-panel 2
plot_tricolor_scatter1(chal_l5_assays_clade2,
                       chal_l5_assays_nonclade2,
                       chal_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm6_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "chal_l5_phitm6_rank6_by_subtype.pdf")


#Fig. S8g sub-panel 3
plot_tricolor_scatter1(chal_l5_assays_clade2,
                       chal_l5_assays_nonclade2,
                       chal_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm1_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "chal_l5_phitm1_rank6_by_subtype.pdf")



#defending correlation


def_l5 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$defending_phage == 'l5')

def_phitm41 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$defending_phage == 'phitm41',select = c(l5_columns))
def_phitm1 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$defending_phage == 'phitm1',select = c(l5_columns))
def_phitm6 <- subset(ave_multiY_confY_lysY,ave_multiY_confY_lysY$defending_phage == 'phitm6',select = c(l5_columns))


names(def_l5) <- paste('l5_',names(def_l5),sep="")
names(def_phitm41) <- paste('phitm41_',names(def_phitm41),sep="")
names(def_phitm1) <- paste('phitm1_',names(def_phitm1),sep="")
names(def_phitm6) <- paste('phitm6_',names(def_phitm6),sep="")



def_l5_assays <- merge(def_l5,def_phitm41,by.x='l5_challenging_phage',by.y='phitm41_challenging_phage',all.x=TRUE,all.y=TRUE)
def_l5_assays <- merge(def_l5_assays,def_phitm1,by.x='l5_challenging_phage',by.y='phitm1_challenging_phage',all.x=TRUE,all.y=TRUE)
def_l5_assays <- merge(def_l5_assays,def_phitm6,by.x='l5_challenging_phage',by.y='phitm6_challenging_phage',all.x=TRUE,all.y=TRUE)


def_l5_assays$l5_phitm41_rank6_diff <- def_l5_assays$phitm41_averaged_rank6 - def_l5_assays$l5_averaged_rank6
def_l5_assays$l5_phitm1_rank6_diff <- def_l5_assays$phitm1_averaged_rank6 - def_l5_assays$l5_averaged_rank6
def_l5_assays$l5_phitm6_rank6_diff <- def_l5_assays$phitm6_averaged_rank6 - def_l5_assays$l5_averaged_rank6


def_l5_assays$phitm1_phitm6_rank6_diff <- def_l5_assays$phitm6_averaged_rank6 - def_l5_assays$phitm1_averaged_rank6

def_l5_assays_clade2 <- subset(def_l5_assays,def_l5_assays$l5_challenging_gene_content_clade == "clade2")
def_l5_assays_nonclade2 <- subset(def_l5_assays,def_l5_assays$l5_challenging_gene_content_clade != "clade2")



#Used as a dummy table for plotting. phiTM4 is lytic, so unable to form a lysogen.
def_l5_assays_clade2_empty <- subset(def_l5_assays,def_l5_assays$l5_defending_phage == "phiTM4")


nrow(def_l5_assays)
nrow(def_l5_assays_clade2)
nrow(def_l5_assays_nonclade2)
nrow(def_l5_assays_clade2_empty)

summary(def_l5_assays$l5_challenging_phage)
summary(def_l5_assays_clade2$l5_challenging_phage)
summary(def_l5_assays_nonclade2$l5_challenging_phage)






#Plots



#Fig. S8h sub-panel 1
plot_tricolor_scatter1(def_l5_assays_clade2,
                       def_l5_assays_nonclade2,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm41_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "def_l5_phitm41_rank6_by_subtype.pdf")


#Fig. S8h sub-panel 2
plot_tricolor_scatter1(def_l5_assays_clade2,
                       def_l5_assays_nonclade2,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm6_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "def_l5_phitm6_rank6_by_subtype.pdf")


#Fig. S8h sub-panel 3
plot_tricolor_scatter1(def_l5_assays_clade2,
                       def_l5_assays_nonclade2,
                       def_l5_assays_clade2_empty,
                       "l5_averaged_rank6",
                       "phitm1_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "def_l5_phitm1_rank6_by_subtype.pdf")


###L5-mutant comparisons using all averaged data above











###Mutant-parent challenging profile comparison
#Pair challenging parent and challenging mutant data
#Split reduced, averaged data into parent and mutant phage tables

#Averaged, non-redundant data should only be for multi-titer assays, no low-confident data, and only from lysogen data.
#There should be no repressor clones or single-titer assays.

mutant_data <- subset(conf_assay_strain_def_chal_average,
                      conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                        conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

#Rename all fields to indicate the data is generated using the mutant challenger
#Note: important to remember that "mutant" prefix refers to the entire immunity assay data, not to any particular column
names(mutant_data) <- paste('mutant_',names(mutant_data),sep="")



#Drop all data from mutant table not involving a mutant challenging phage
mutant_data <- subset(mutant_data,mutant_data$mutant_challenging_source == 'lab')




#Create parent data to match
parent_data <- subset(conf_assay_strain_def_chal_average,
                      conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                        conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

names(parent_data) <- paste('parent_',names(parent_data),sep="")


#Now match parent challenging data to mutant challenging data
mutant_data$parent_defending_challenging <- paste(mutant_data$mutant_defending_phage,"_",mutant_data$mutant_challenging_parent,sep="")
mutant_analysis <- merge(mutant_data,parent_data,by.x="parent_defending_challenging", by.y="parent_defending_challenging")





#Compute difference in infection profiles
#Mutant phenotype should be stronger than parent, so don't use absolute value
mutant_analysis$averaged_rank6_diff <- mutant_analysis$mutant_averaged_rank6 - mutant_analysis$parent_averaged_rank6







setwd("~/scratch/immunity_analysis/output/")
write.table(mutant_analysis,
            "mutant_analysis.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)








#QC
plot_hist1(mutant_analysis,
           "averaged_rank6_diff",
           10,
           c(-1,4),
           c(0,80),
           "mutant_averaged_rank6_diff.pdf")





#Reduce to only escape mutants

#Mutant-parent phage sets
#1
#pioneer
#phitm35 (escape mutant)

#2
#eagleeye
#phitm36 (escape mutant)

#3
#phitm44 (isolated mutant)
#phitm38 (escape mutant)

#4
#et2brutus
#phitm39 (escape mutant)	
#phitm40 (escape mutant)

#5
#l5
#phitm41 (escape mutant)

#6
#trixie
#phitm42 (escape mutant)

#8 (clone EM)
#davinci
#phitm46 (escape mutant)

#9 (clone EM)
#gladiator
#phitm47 (escape mutant)


escape_mutant_analysis <- subset(mutant_analysis,
                                 mutant_analysis$mutant_challenging_phage == "phitm35" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm36" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm38" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm39" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm40" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm41" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm42" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm46" | 
                                   mutant_analysis$mutant_challenging_phage == "phitm47")


phitm35_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm35")
phitm36_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm36")
phitm38_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm38")
phitm39_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm39")
phitm40_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm40")
phitm41_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm41")
phitm42_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm42")
phitm46_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm46")
phitm47_mutant_analysis <- subset(mutant_analysis,mutant_analysis$mutant_challenging_phage == "phitm47")



escape_mutant_analysis_interclade <- subset(escape_mutant_analysis,
                                            escape_mutant_analysis$parent_gene_content_clade_compare == "different")
escape_mutant_analysis_intraclade2 <- subset(escape_mutant_analysis,
                                             escape_mutant_analysis$parent_gene_content_clade_compare == "clade2")

#Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
escape_mutant_analysis_intraclade2_empty <- subset(escape_mutant_analysis,
                                             escape_mutant_analysis$parent_gene_content_clade_compare == "empty")


escape_mutant_analysis_intraclade2_parenthomotypic <- subset(escape_mutant_analysis_intraclade2,
                                                       as.character(escape_mutant_analysis_intraclade2$parent_challenging_phage) ==
                                                         as.character(escape_mutant_analysis_intraclade2$parent_defending_phage))
escape_mutant_analysis_intraclade2_parentheterotypic <- subset(escape_mutant_analysis_intraclade2,
                                                             as.character(escape_mutant_analysis_intraclade2$parent_challenging_phage) !=
                                                               as.character(escape_mutant_analysis_intraclade2$parent_defending_phage))

escape_mutant_analysis_intraclade2_mutanthomotypic <- subset(escape_mutant_analysis_intraclade2,
                                                             as.character(escape_mutant_analysis_intraclade2$mutant_challenging_phage) ==
                                                               as.character(escape_mutant_analysis_intraclade2$mutant_defending_phage))
escape_mutant_analysis_intraclade2_mutantheterotypic <- subset(escape_mutant_analysis_intraclade2,
                                                               as.character(escape_mutant_analysis_intraclade2$mutant_challenging_phage) !=
                                                                 as.character(escape_mutant_analysis_intraclade2$mutant_defending_phage))

#QC
nrow(escape_mutant_analysis_intraclade2)
nrow(escape_mutant_analysis_intraclade2_empty)
nrow(escape_mutant_analysis_intraclade2_parenthomotypic)
nrow(escape_mutant_analysis_intraclade2_parentheterotypic)
nrow(escape_mutant_analysis_intraclade2_mutanthomotypic)
nrow(escape_mutant_analysis_intraclade2_mutantheterotypic)


summary(escape_mutant_analysis_intraclade2_parenthomotypic$parent_defending_phage)
summary(escape_mutant_analysis_intraclade2_mutanthomotypic$mutant_defending_phage)
#







#Plots


#Fig. 7c sub-panel 1
plot_tricolor_scatter1(escape_mutant_analysis_intraclade2,
                       escape_mutant_analysis_interclade,
                       escape_mutant_analysis_intraclade2_empty,
                       "parent_averaged_rank6",
                       "mutant_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "ave_rank6_parent_vs_mutant_subtypes.pdf")


#Fig. 7f sub-panel 1
plot_tricolor_scatter2(escape_mutant_analysis_intraclade2_parentheterotypic,
                       escape_mutant_analysis_interclade,
                       escape_mutant_analysis_intraclade2_parenthomotypic,
                       "parent_stoperator_pwd_dist_euc",
                       "parent_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "ave_rank6_parent_vs_parent_stopEuc.pdf")


#Fig. 7f sub-panel 2
plot_tricolor_scatter2(escape_mutant_analysis_intraclade2_mutantheterotypic,
                       escape_mutant_analysis_interclade,
                       escape_mutant_analysis_intraclade2_mutanthomotypic,
                       "mutant_stoperator_pwd_dist_euc",
                       "mutant_averaged_rank6",
                       c(0,5),
                       c(0,6),
                       "ave_rank6_mutant_vs_mutant_stopEuc.pdf")






#Compare escape mutants to parents on repressor clone strains

mutant_data_rep <- subset(conf_assay_strain_def_chal_average,
                      conf_assay_strain_def_chal_average$strain_type == 'repressor_clone' & 
                        conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')



#Rename all fields to indicate the data is generated using the mutant challenger
#Note: important to remember that "mutant" prefix refers to the entire immunity assay data, not to any particular column
names(mutant_data_rep) <- paste('mutant_',names(mutant_data_rep),sep="")



#Drop all data from mutant table not involving a mutant challenging phage
mutant_data_rep <- subset(mutant_data_rep,mutant_data_rep$mutant_challenging_source == 'lab')




#Create parent data to match
parent_data_rep <- subset(conf_assay_strain_def_chal_average,
                      conf_assay_strain_def_chal_average$strain_type == 'repressor_clone' & 
                        conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

names(parent_data_rep) <- paste('parent_',names(parent_data_rep),sep="")


#Now match parent challenging data to mutant challenging data
mutant_data_rep$parent_defending_challenging <- paste(mutant_data_rep$mutant_defending_phage,"_",mutant_data_rep$mutant_challenging_parent,sep="")
mutant_rep_analysis <- merge(mutant_data_rep,parent_data_rep,by.x="parent_defending_challenging", by.y="parent_defending_challenging")





#Compute difference in infection profiles
#Mutant phenotype should be stronger than parent, so don't use absolute value
mutant_rep_analysis$averaged_rank6_diff <- mutant_rep_analysis$mutant_averaged_rank6 - mutant_rep_analysis$parent_averaged_rank6

# mutant_rep_analysis$averaged_infection_strength_diff <- mutant_rep_analysis$mutant_averaged_infection_strength - mutant_rep_analysis$parent_averaged_infection_strength
# mutant_rep_analysis$averaged_four_factors_diff <- mutant_rep_analysis$mutant_averaged_four_factors - mutant_rep_analysis$parent_averaged_four_factors




#Reduce to only escape mutants

escape_mutant_rep_analysis <- subset(mutant_rep_analysis,
                                     mutant_rep_analysis$mutant_challenging_phage == "phitm35" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm36" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm38" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm39" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm40" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm41" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm42" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm46" | 
                                       mutant_rep_analysis$mutant_challenging_phage == "phitm47")


phitm35_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm35")
phitm36_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm36")
phitm38_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm38")
phitm39_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm39")
phitm40_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm40")
phitm41_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm41")
phitm42_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm42")
phitm46_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm46")
phitm47_mutant_rep_analysis <- subset(mutant_rep_analysis,mutant_rep_analysis$mutant_challenging_phage == "phitm47")



escape_mutant_rep_analysis_intraclade2 <- subset(escape_mutant_rep_analysis,
                                                 escape_mutant_rep_analysis$parent_gene_content_clade_compare == "clade2")
escape_mutant_rep_analysis_interclade <- subset(escape_mutant_rep_analysis,
                                                 escape_mutant_rep_analysis$parent_gene_content_clade_compare == "different")


#Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
escape_mutant_rep_analysis_intraclade2_empty <- subset(escape_mutant_rep_analysis,
                                                 escape_mutant_rep_analysis$parent_gene_content_clade_compare == "empty")


#QC
nrow(escape_mutant_rep_analysis_intraclade2_empty)




#Fig. 7c sub-panel 2
plot_tricolor_scatter1(escape_mutant_rep_analysis_intraclade2,
                       escape_mutant_rep_analysis_interclade,
                       escape_mutant_rep_analysis_intraclade2_empty,
                       "parent_averaged_rank6",
                       "mutant_averaged_rank6",
                       c(0,6),
                       c(0,6),
                       "ave_rank6_parent_vs_mutant_rep_subtypes.pdf")


###Compare known empirical temperate to isolated mutant to escape mutant


















###Whole genome metrics
#Assess correlation of genomic distance, gcd, stoperator_pwd, repressor phylogeny, etc.
#This should be independent of any immunity assay data





#Create list of all phage pairs  
#actino1319_phages <- levels(phage_metadata$phageid)
actino1321_phages <- levels(phage_metadata$phageid)

distance_metrics <- expand.grid(phage1 = actino1321_phages,
                                phage2 = actino1321_phages)

distance_metrics$phage1_phage2 <- paste(distance_metrics$phage1,
                                        "_",
                                        distance_metrics$phage2,
                                        sep="")

#This dataset now contains reciprocal and self-comparison data
#To solve this, remove all rows in which phage1 and phage2 are not alphabetically ordered.
#This does NOT retain self-comparisons (e.g. alma_alma)
distance_metrics$alpha_ordered <- as.character(distance_metrics$phage1) < as.character(distance_metrics$phage2)

#Now retain only the unique pairwise comparisons
distance_metrics <- subset(distance_metrics,distance_metrics$alpha_ordered == TRUE)



#Match the genome/gene distance data. Many comparisons will not be matched
distance_metrics <- merge(distance_metrics,genomic_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,repressor336_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,cas4311_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,endovii306_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,dnapol311_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,portal311_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,stoperator_pwm_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,immunity_correlation_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)



#Match the phage metadata
phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('phage1','_',names(phage_metadata_to_match),sep="")
distance_metrics <- merge(distance_metrics,phage_metadata_to_match,by.x="phage1",by.y="phage1_phageid")

phage_metadata_to_match <- phage_metadata
names(phage_metadata_to_match) <- paste('phage2','_',names(phage_metadata_to_match),sep="")
distance_metrics <- merge(distance_metrics,phage_metadata_to_match,by.x="phage2",by.y="phage2_phageid")



#Compute comparisons


#Since this dataset is not limited to immunity assays, it contains inter-cluster comparisons.
distance_metrics$cluster_compare <- ifelse(distance_metrics$phage1_cluster==distance_metrics$phage2_cluster,
                                           as.character(distance_metrics$phage1_cluster),
                                           "different")

distance_metrics$subcluster_compare <- ifelse(distance_metrics$phage1_subcluster==distance_metrics$phage2_subcluster,
                                              as.character(distance_metrics$phage1_subcluster),
                                              "different")

distance_metrics$source_compare <- ifelse(distance_metrics$phage1_source==distance_metrics$phage2_source,
                                          as.character(distance_metrics$phage1_source),
                                          "different")

distance_metrics$temperate_empirical_compare <- ifelse(distance_metrics$phage1_cluster_a_temperate_empirical==distance_metrics$phage2_cluster_a_temperate_empirical,
                                                       as.character(distance_metrics$phage1_cluster_a_temperate_empirical),
                                                       "different")

distance_metrics$functional_repressor_compare <- ifelse(distance_metrics$phage1_cluster_a_functional_repressor_predicted==distance_metrics$phage2_cluster_a_functional_repressor_predicted,
                                                        as.character(distance_metrics$phage1_cluster_a_functional_repressor_predicted),
                                                        "different")

distance_metrics$lysogen_type_compare <- ifelse(distance_metrics$phage1_lysogen_type==distance_metrics$phage2_lysogen_type,
                                                as.character(distance_metrics$phage1_lysogen_type),
                                                "different")

distance_metrics$integrase_compare <- ifelse(distance_metrics$phage1_pham_integrase==distance_metrics$phage2_pham_integrase,
                                             as.character(distance_metrics$phage1_pham_integrase),
                                             "different")

distance_metrics$parb_compare <- ifelse(distance_metrics$phage1_pham_parb==distance_metrics$phage2_pham_parb,
                                        as.character(distance_metrics$phage1_pham_parb),
                                        "different")

distance_metrics$repressor_hth_compare <- stringdist(as.character(distance_metrics$phage1_repressor_hth_domain_sequence),
                                                     as.character(distance_metrics$phage2_repressor_hth_domain_sequence),
                                                     method="hamming")

distance_metrics$repressor_length_full_compare <- abs(distance_metrics$phage1_repressor_length_full - distance_metrics$phage2_repressor_length_full)
distance_metrics$repressor_length_nterm_compare <- abs(distance_metrics$phage1_repressor_length_nterm - distance_metrics$phage2_repressor_length_nterm)
distance_metrics$repressor_length_cterm_compare <- abs(distance_metrics$phage1_repressor_length_cterm - distance_metrics$phage2_repressor_length_cterm)


distance_metrics$subcluster_compare2 <- ifelse(distance_metrics$phage1_subcluster==distance_metrics$phage2_subcluster,
                                            "same",
                                            "different")

distance_metrics$gene_content_clade_compare <- ifelse(distance_metrics$phage1_gene_content_clade==distance_metrics$phage2_gene_content_clade,
                                                      as.character(distance_metrics$phage1_gene_content_clade),
                                                      "different")

distance_metrics$gene_content_clade_compare2 <- ifelse(distance_metrics$phage1_gene_content_clade==distance_metrics$phage2_gene_content_clade,
                                               "same",
                                               "different")


distance_metrics$cluster_compare <- as.factor(distance_metrics$cluster_compare)
distance_metrics$subcluster_compare <- as.factor(distance_metrics$subcluster_compare)
distance_metrics$source_compare <- as.factor(distance_metrics$source_compare)
distance_metrics$temperate_empirical_compare <- as.factor(distance_metrics$temperate_empirical_compare)
distance_metrics$functional_repressor_compare <- as.factor(distance_metrics$functional_repressor_compare)
distance_metrics$lysogen_type_compare <- as.factor(distance_metrics$lysogen_type_compare)
distance_metrics$integrase_compare <- as.factor(distance_metrics$integrase_compare)
distance_metrics$parb_compare <- as.factor(distance_metrics$parb_compare)
distance_metrics$subcluster_compare2 <- factor(distance_metrics$subcluster_compare2,c("same","different"))
distance_metrics$gene_content_clade_compare <- as.factor(distance_metrics$gene_content_clade_compare)  
distance_metrics$gene_content_clade_compare2 <- factor(distance_metrics$gene_content_clade_compare2,c("same","different"))








#Plot relationship between repressors, stoperators, and genomes




#Keep only Mycobacteriophage data and drop Gordonia phages
#Do not include escape mutants.
clusterA_data <- subset(distance_metrics,
                        distance_metrics$cluster_compare == "A" &
                          distance_metrics$source_compare == "environment" &
                          distance_metrics$phage1_host == "Mycobacterium" &
                          distance_metrics$phage2_host == "Mycobacterium")


clusterA_clade2 <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "clade2")

clusterA_clade2_diff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different" &
                                 (clusterA_data$phage1_gene_content_clade == 'clade2' |
                                    clusterA_data$phage2_gene_content_clade == 'clade2'))

#Used as a dummy table for plotting. "empty" is not a valid clade_comparison.
clusterA_clade2_empty <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "empty")

#QC
nrow(clusterA_clade2_empty)








#Plots
setwd("~/scratch/immunity_analysis/output/")


#Fig. 2c
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_clade2,
                      "modified_mash_distance",
                      "pham_pham_dissimilarity",
                      c(0,0.5),
                      c(0,1),
                      "mash_vs_gcd_clade2.pdf")


#Fig. 2d
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_clade2,
                      "pham_pham_dissimilarity",
                      "repressor_full_mafft_dist_uncorrected",
                      c(0,1),
                      c(0,70),
                      "gcd_vs_repFullMafUn_clade2.pdf")


#Fig. 2f sub-panel 1
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_clade2,
                      "pham_pham_dissimilarity",
                      "stoperator_pwd_dist_euc",
                      c(0,1),
                      c(0,5),
                      "gcd_vs_stopEuc_clade2.pdf")


#Fig. 2f sub-panel 2
plot_bicolor_scatter1(clusterA_clade2_diff,
                      clusterA_clade2,
                      "repressor_full_mafft_dist_uncorrected",
                      "stoperator_pwd_dist_euc",
                      c(0,70),
                      c(0,5),
                      "repMafftUn_vs_stopEuc_clade2.pdf")


#Fig. 9a
plot_bicolor_scatter2(clusterA_clade2_diff,
                      clusterA_clade2,
                      "repressor_nterm_mafft_dist_uncorrected",
                      "repressor_cterm_mafft_dist_uncorrected",
                      c(0,70),
                      c(0,70),
                      "repNtermMafftUn_vs_repCtermMaffUn_clade2.pdf")


#Fig. 5c sub-panel 1
plot_tricolor_scatter2(clusterA_clade2,
                       clusterA_clade2_diff,
                       clusterA_clade2_empty,
                       "stoperator_pwd_dist_euc",
                       "challenging_cor_reduced",
                       c(0,5),
                       c(-1,1),
                       "clusterA_clade2_challCor_vs_stopEuc_subtypes.pdf")


#Fig. 5c sub-panel 2
plot_tricolor_scatter2(clusterA_clade2,
                       clusterA_clade2_diff,
                       clusterA_clade2_empty,
                       "stoperator_pwd_dist_euc",
                       "defending_cor_reduced",
                       c(0,5),
                       c(-1,1),
                       "clusterA_clade2_defCor_vs_stopEuc_subtypes.pdf")



#Repressor size
clusterA_subset <- subset(phage_metadata,
                          phage_metadata$cluster == 'A' &
                            phage_metadata$source == 'environment' &
                            phage_metadata$host == 'Mycobacterium' &
                            phage_metadata$cluster_a_functional_repressor_predicted == 'yes' &
                            phage_metadata$gene_content_clade == 'clade2')


#Fig. S1a
par(mar=c(4,8,16,20))
boxplot(clusterA_subset$repressor_length_full,
        las=1,cex.axis=2,ann=FALSE,main=NULL,outline=FALSE,ylim=c(150,250),
        col="light grey")
par(new=TRUE)
stripchart(clusterA_subset$repressor_length_full,
           vertical=TRUE,las=1,cex.axis=2,pch=16,method="jitter",cex=1,
           ann=FALSE,main=NULL,ylim=c(150,250))
dev.copy(pdf,"clusterA_subset_repFull_sizes.pdf")
dev.off()



### Whole genome metrics (regardless of immunity data) above





