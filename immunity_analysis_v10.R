#Analyze pairwise distances for immunity assay analysis
#Use Actino1321 data, in which escape mutants have been added and Che12 sequence has been corrected.

#setwd("~/Desktop/Project_stoperator/immunity_bfx/")




#The melt function of reshape2 package is needed to convert matrix to unique-pair table
#install.packages("reshape2")
#TODO: update R and RStudio to most appropriate version for reshape2
library(reshape2)

#Stringdist is needed to compute hamming distance between strings
#install.packages("stringdist")
library(stringdist)


#Pysch is needed to produce 2d error bar plots
#install.packages("psych")
library(psych)

setwd("~/scratch/immunity_analysis/input/")



#Functions

#Compute comparison fields
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




plot_scored_infection_strength <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,
       as.numeric(as.character(table$scored_infection_strength)),
       xlim=c(0,0.4),ylim=c(0,3),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}



plot_scored_four_factors_mash <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,
       as.numeric(as.character(table$scored_four_factors)),
       xlim=c(0,0.4),ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}

plot_scored_four_factors_gcd <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$pham_pham_dissimilarity,
       as.numeric(as.character(table$scored_four_factors)),
       ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}

plot_scored_four_factors_repressor <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$repressor_muscle_bionj_distances,
       as.numeric(as.character(table$scored_four_factors)),
       ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}

plot_scored_four_factors_repressor_phyml <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$repressor_prank_phyml_distances,
       as.numeric(as.character(table$scored_four_factors)),
       ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}


plot_scored_four_factors_recb <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$recb_muscle_bionj_distances,
       as.numeric(as.character(table$scored_four_factors)),
       ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}

plot_scored_four_factors_portal <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$portal_muscle_bionj_distances,
       as.numeric(as.character(table$scored_four_factors)),
       ylim=c(0,9),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
}


plot_reduce_data1 <- function(table,source1,source2,lys1,lys2){
  
  table <- subset(table,
                  table$defending_source == source1
                  & table$challenging_source == source2
                  & table$defending_lysogen_type == lys1
                  & table$challenging_lysogen_type == lys2
  )
  
  plot_scored_four_factors(table)
  plot_name <- paste(source1,"_",source2,"_",lys1,"_",lys2,'.pdf',sep='')
  dev.copy(pdf,plot_name)
  dev.off()
  
}



###End of functions























###Import datasets

#immunity data
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
# "scored_infection_strength"
# "scored_turbidity"
# "scored_plaque_size"
# "scored_plaques"
# "scored_four_factors"
# "rank6"        


immunity_data <- read.csv("immunity_data.csv",sep=",",header=TRUE)
immunity_data$immunity_assay_id <- as.factor(immunity_data$immunity_assay_id)
immunity_data$immunity_set <- as.factor(immunity_data$immunity_set)
immunity_data$notebook <- as.factor(immunity_data$notebook)
immunity_data$page <- as.factor(immunity_data$page)
immunity_data$lawn_reliability <- as.factor(immunity_data$lawn_reliability)
immunity_data$phage_reliability <- as.factor(immunity_data$phage_reliability)
immunity_data$scored_infection_strength <- as.factor(immunity_data$scored_infection_strength)
immunity_data$scored_turbidity <- as.factor(immunity_data$scored_turbidity)
immunity_data$scored_plaque_size <- as.factor(immunity_data$scored_plaque_size)
immunity_data$scored_plaques <- as.factor(immunity_data$scored_plaques)
immunity_data$scored_four_factors <- as.factor(immunity_data$scored_four_factors)
immunity_data$rank6 <- as.factor(immunity_data$rank6)

#Several fields contain 'unspecified' which can be converted to NA
#Afterwards, re-factor
#prophage
#repressor_clone
#tested_titer
#phage_reliability
#scored_infection_strength
immunity_data[immunity_data == "Unspecified"] <- NA
immunity_data[immunity_data == "unspecified"] <- NA
immunity_data$prophage <- factor(immunity_data$prophage)
immunity_data$repressor_clone <- factor(immunity_data$repressor_clone)
immunity_data$phage_reliability <- factor(immunity_data$phage_reliability)
immunity_data$scored_infection_strength <- factor(immunity_data$scored_infection_strength)
immunity_data$rank6 <- factor(immunity_data$rank6)

#Convert titer to numeric
immunity_data$tested_titer <- as.numeric(as.character(immunity_data$tested_titer))


#These fields contain NA's, but these are descriptive columns so no need to conver them to Unspecified
#observed_infection_strength = NA
#observed_turbidity = NA
#observed_plaque_size = NA
#observed_plaques = NA


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



#mash genomic distance data
#all actino1321 phages 
#reciprocal data and self-comparison data
# "phage1_phage2"
# "modified_mash_distance"
# "pham_pham_dissimilarity"
genomic_distance_data <- read.csv("genomic_distance_data.csv",sep=",",header=TRUE)




#actino1321 phage metadata
# "phageid"
# "host"
# "cluster"
# "subcluster"
# "size"
# "status"
# "author"
# "datelastmodified"
# "mode_approx_80_percent"
# "network005_interaction_tally"
# "network005_group"
# "network005_group_tally"
# "lysogen_type"
# "pham_integrase"
# "pham_para" (imported as int)
# "source"
# "parent"
# "cluster_a_functional_repressor_predicted"
# "cluster_a_temperate_empirical"
# "maxgcdgap_all"
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
#datelastmodified = Unspecified
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
phage_metadata$datelastmodified <- factor(phage_metadata$datelastmodified)
phage_metadata$gene_content_clade <- factor(phage_metadata$gene_content_clade)
phage_metadata$repressor_length_full <- as.numeric(as.character(phage_metadata$repressor_length_full))
phage_metadata$repressor_length_nterm <- as.numeric(as.character(phage_metadata$repressor_length_nterm))
phage_metadata$repressor_length_cterm <- as.numeric(as.character(phage_metadata$repressor_length_cterm))
phage_metadata$pleft_alignment_reference <- as.numeric(as.character(phage_metadata$pleft_alignment_reference))
phage_metadata$immunity_repressor_alignment_reference <- as.numeric(as.character(phage_metadata$immunity_repressor_alignment_reference))
phage_metadata$genome_center_alignment_reference <- as.numeric(as.character(phage_metadata$genome_center_alignment_reference))





#actino1321 repressor protein distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
#no data for parent phages that are natural mutants (e.g. d29, misswhite, jeffabunny) with no repressor annotated
# "phage1_phage2"
# "repressor_muscle_bionj_distances"
# "repressor_prank_phyml_distances"
repressor_distance_data <- read.csv("repressor_distance_data.csv",sep=",",header=TRUE)











#actino1321 portal protein distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "portal_muscle_bionj_distances"
# "portal_prank_phyml_distances"
portal_distance_data <- read.csv("portal_distance_data.csv",sep=",",header=TRUE)


#actino1321 recb protein distance data
#reciprocal data and self-comparison data
#only Cluster A parent phages used in immunity assays
#no data for escape mutants
# "phage1_phage2"
# "recb_muscle_bionj_distances"
# "recb_prank_phyml_distances"
recb_distance_data <- read.csv("recb_distance_data.csv",sep=",",header=TRUE)






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
#contains list of prediction stopeator sites in all 327 Cluster A genomes from actino1321, for each of the 327 stoperaotr PWMs
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



par(mar=c(4,8,4,4))
hist(stop88_clade2_env_self$site_pleft_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
dev.copy(pdf,'clade2_env_self_sites_near_pleft2.pdf')
dev.off()




# par(mar=c(4,8,4,4))
# hist(stop88_clade2_env_self$site_pleft_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=1000,xlim=c(-10000,500))
# abline(v=-1500)
# 

#FinalFig
par(mar=c(4,8,8,4))
hist(stop88_clade2_env_self$site_pleft_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=2500,xlim=c(-2000,500),ylim=c(0,100))
dev.copy(pdf,'clade2_env_self_sites_near_pleft.pdf')
dev.off()


nrow(subset(stop88_clade2_env_self,stop88_clade2_env_self$site_pleft_dist > -1000))/nrow(stop88_clade2_env_self)







#Distance from right end of genome
hist(stop88_clade2_env_self$site_right_end_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
hist(stop88_clade2_env_self$site_right_end_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=1000,xlim=c(-2500,0))
abline(v=-2000)


#FinalFig
par(mar=c(4,8,4,4))
hist(stop88_clade2_env_self$site_right_end_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
dev.copy(pdf,'clade2_env_self_sites_right_termini.pdf')
dev.off()



nrow(subset(stop88_clade2_env_self,stop88_clade2_env_self$site_right_end_dist > -2000))/nrow(stop88_clade2_env_self)
#22% of all sites are within 2kb of the right end of the genome






#Distance from repressor
# par(mar=c(4,8,8,4))
# hist(stop88_clade2_env_self$site_rep_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
# 
# par(mar=c(4,8,8,4))
# hist(stop88_clade2_env_self$site_rep_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=2000,xlim=c(-1000,1000),ylim=c(0,50))
# dev.copy(pdf,'clade2_env_self_sites_near_repressor.pdf')
# dev.off()

par(mar=c(4,8,8,4))
hist(stop88_clade2_env_self$site_rep_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=2000,xlim=c(-4000,1000),ylim=c(0,50))
dev.copy(pdf,'clade2_env_self_sites_near_repressor2.pdf')
dev.off()







#Distance from center of genome

hist(stop88_clade2_env_self$site_center_dist,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)

nrow(subset(stop88_clade2_env_self,stop88_clade2_env_self$site_center_dist > 0))/nrow(stop88_clade2_env_self)
#74% of all sites are to the right of the genome center







#Distance from left end of genome
par(mar=c(4,8,8,4))
hist(stop88_clade2_env_self$tfbs88_site_middle,main=NULL,ann=FALSE,las=1,cex.axis=2,col="black",breaks=100)
















#Examine # of stoperators per genome
stop88_clade2_env_self_freq <- as.data.frame(table(stop88_clade2_env_self$tfbs88_stoperator_target))
names(stop88_clade2_env_self_freq) <- c("phage","frequency")
  







setwd("~/scratch/immunity_analysis/output/")

#FinalFig
par(mar=c(4,8,8,4))
hist(stop88_clade2_env_self_freq$frequency,col="black",breaks=50,xlim=c(0,50),ylim=c(0,15),main=NULL,ann=FALSE,las=1,cex.axis=2)
dev.copy(pdf,'clade2_env_self_stops_per_genome.pdf')
dev.off()













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

#FinalFig
par(mar=c(4,8,8,4))
plot(stop88_clade2_env_self_freq_summary$right_sites_reverse_percent,
     stop88_clade2_env_self_freq_summary$left_sites_forward_percent,
     xlim=c(0,1),ylim=c(0,1),
     cex.axis=2,ann=FALSE,main=NULL,las=1,
     col="black",pch=16,cex=2)
abline(0,1)
dev.copy(pdf,'clade2_env_self_percent_txn_oriented_sites_per_genome.pdf')
dev.off()











#Table creates a two column table = the levels and the frequency (int)


#Compute # of sites within 2kb of the right genome end
stoperator_site_predictions88_2kbright <- subset(stoperator_site_predictions88,
                                                 stoperator_site_predictions88$site_right_end_dist >= -2000)

stoperator_site_predictions88_2kbleft <- subset(stoperator_site_predictions88,
                                                 stoperator_site_predictions88$site_right_end_dist < -2000)








stoperator_site_predictions88_2kbright_targetself <- subset(stoperator_site_predictions88_2kbright,
                                                            stoperator_site_predictions88_2kbright$tfbs88_site_self == TRUE)

stoperator_site_predictions88_2kbright_nontargetself <- subset(stoperator_site_predictions88_2kbright,
                                                            stoperator_site_predictions88_2kbright$tfbs88_site_self == FALSE)





stoperator_site_predictions88_2kbleft_targetself <- subset(stoperator_site_predictions88_2kbleft,
                                                           stoperator_site_predictions88_2kbleft$tfbs88_site_self == TRUE)

stoperator_site_predictions88_2kbleft_nontargetself <- subset(stoperator_site_predictions88_2kbleft,
                                                              stoperator_site_predictions88_2kbleft$tfbs88_site_self == FALSE)








stop88_tally <- as.data.frame(table(stoperator_site_predictions88$tfbs88_motif_target))
stop88_2kbright_tally <- as.data.frame(table(stoperator_site_predictions88_2kbright$tfbs88_motif_target))
stop88_2kbleft_tally <- as.data.frame(table(stoperator_site_predictions88_2kbleft$tfbs88_motif_target))


stop88_2kbright_targetself_tally <- as.data.frame(table(stoperator_site_predictions88_2kbright_targetself$tfbs88_motif_target))
stop88_2kbright_nontargetself_tally <- as.data.frame(table(stoperator_site_predictions88_2kbright_nontargetself$tfbs88_motif_target))

stop88_2kbleft_targetself_tally <- as.data.frame(table(stoperator_site_predictions88_2kbleft_targetself$tfbs88_motif_target))
stop88_2kbleft_nontargetself_tally <- as.data.frame(table(stoperator_site_predictions88_2kbleft_nontargetself$tfbs88_motif_target))


names(stop88_tally) <- c("tfbs88_motif_target","stoperator88_tally")
names(stop88_2kbright_tally) <- c("tfbs88_motif_target","stoperator88_2kbright_tally")
names(stop88_2kbleft_tally) <- c("tfbs88_motif_target","stoperator88_2kbleft_tally")

names(stop88_2kbright_targetself_tally) <- c("tfbs88_motif_target","stoperator88_2kbright_targetself_tally")
names(stop88_2kbright_nontargetself_tally) <- c("tfbs88_motif_target","stoperator88_2kbright_nontargetself_tally")
names(stop88_2kbleft_targetself_tally) <- c("tfbs88_motif_target","stoperator88_2kbleft_targetself_tally")
names(stop88_2kbleft_nontargetself_tally) <- c("tfbs88_motif_target","stoperator88_2kbleft_nontargetself_tally")


stoperator_tally <- merge(stop88_tally,stop88_2kbright_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")
stoperator_tally <- merge(stoperator_tally,stop88_2kbleft_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")
stoperator_tally <- merge(stoperator_tally,stop88_2kbright_targetself_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")
stoperator_tally <- merge(stoperator_tally,stop88_2kbright_nontargetself_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")
stoperator_tally <- merge(stoperator_tally,stop88_2kbleft_targetself_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")
stoperator_tally <- merge(stoperator_tally,stop88_2kbleft_nontargetself_tally,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")



#Remove all NA values in final table
stoperator_tally[is.na(stoperator_tally)] <- 0



#Add individual target and motif columns
stoperator_tally <- merge(stoperator_tally,all_levels,by.x="tfbs88_motif_target",by.y="tfbs88_motif_target")

#stoperator_tally now contains tally of all sites predicted from all possible pairwise combinations of PWM and target genome.
















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
main_immunity_data <- merge(main_immunity_data,repressor_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,portal_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,recb_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,repressor336_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,cas4311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,endovii306_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,dnapol311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)
main_immunity_data <- merge(main_immunity_data,portal311_distance_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)



#Match PWM distance data. Contains PWM data for 264 phages.
main_immunity_data <- merge(main_immunity_data,stoperator_pwm_data,by.x="defending_challenging",by.y="phage1_phage2",all.x=TRUE)



#TODO in progress
#Match stoperator site prediction data. Contains site tallies at five score thresholds (80, 85, 90, 95, 1), 
#actino1319 = between 248 PWMs and all 325 actino1319 Cluster A genomes
#actino1321 = between 264 PWMs and all 327 actino1321 Cluster A genomes
#TODO: update stoperator prediction data for phiTM1 and phiTM33, now in actino1321
#Need to keep all rows. phiTM41 and phiTM6 PWMs are not in the stoperator data, so this data would be lost otherwise. 
main_immunity_data <- merge(main_immunity_data,stoperator_tally,by.x="defending_challenging",by.y="tfbs88_motif_target",all.x=TRUE)



#Compute comparison fields
main_immunity_data <- compute_comparisons(main_immunity_data)





#Retain only confident data, and discard questionable data

main_immunity_data_unreduced <- main_immunity_data

main_immunity_data <- subset(main_immunity_data,
                          main_immunity_data$lawn_reliability != 1 &
                            main_immunity_data$phage_reliability != 1)




###At this point, all data in main_immunity_data is derived from phages present in actino1321 database AND only confident data###

























#Histogram of titers to assess the range of titers used
setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,15,4))
hist(log(main_immunity_data$tested_titer,10),xlim=c(0,10),col='black',ann=FALSE,main=NULL,las=1,breaks=20)
dev.copy(pdf,'titer_hist.pdf')
dev.off()


#Plot immunity phenotypes by titer to assess which phenotypes are observed at which range of titers
par(mar=c(4,8,4,4))
plot(log(main_immunity_data$tested_titer,10),
     as.numeric(as.character(main_immunity_data$rank6)),
     xlim=c(0,9),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'titer_vs_rank6.pdf')
dev.off()

###




































###Average data

convert_immunity_scores_to_numeric <- function(table){
  
  table$scored_infection_strength <- as.numeric(as.character(table$scored_infection_strength))
  table$scored_turbidity <- as.numeric(as.character(table$scored_turbidity))
  table$scored_plaque_size <- as.numeric(as.character(table$scored_plaque_size))
  table$scored_plaques <- as.numeric(as.character(table$scored_plaques))
  table$scored_four_factors <- as.numeric(as.character(table$scored_four_factors))
  table$rank6 <- as.numeric(as.character(table$rank6))
  return(table)
}

convert_immunity_scores_to_factor <- function(table){
  
  table$scored_infection_strength <- factor(as.character(table$scored_infection_strength))
  table$scored_turbidity <- factor(as.character(table$scored_turbidity))
  table$scored_plaque_size <- factor(as.character(table$scored_plaque_size))
  table$scored_plaques <- factor(as.character(table$scored_plaques))
  table$scored_four_factors <- factor(as.character(table$scored_four_factors))
  table$rank6 <- factor(as.character(table$rank6))
  return(table)
}


#For certain analyses, I need averaged non-duplicated immunity comparisons. 
#Averages can be computed by unique assay_strain_defending_challenging identifier, 
#or by unique strain_defending_challenging identifier (which merges multiple-titer and single-titer data)


#No need for this since main_immunity_data now retains only confident data
# conf_to_average <- subset(main_immunity_data,
#                           main_immunity_data$lawn_reliability != 1 &
#                             main_immunity_data$phage_reliability != 1)

conf_to_average <- main_immunity_data
conf_to_average$assay_strain_defending_challenging <- factor(conf_to_average$assay_strain_defending_challenging)
conf_to_average$strain_defending_challenging <- factor(conf_to_average$strain_defending_challenging)
conf_to_average <- convert_immunity_scores_to_numeric(conf_to_average)





#Average infection scores for each unique assay_strain_defending_challenging identifier
conf_assay_strain_def_chal_average <- aggregate(conf_to_average[,c('scored_infection_strength',
                                                                   'scored_turbidity',
                                                                   'scored_plaque_size',
                                                                   'scored_plaques',
                                                                   'scored_four_factors',
                                                                   'rank6')],
                                                list(conf_to_average$assay_strain_defending_challenging),mean)

names(conf_assay_strain_def_chal_average) <- c('assay_strain_defending_challenging',
                                               'averaged_infection_strength',
                                               'averaged_turbidity',
                                               'averaged_plaque_size',
                                               'averaged_plaques',
                                               'averaged_four_factors',
                                               'averaged_rank6') 

conf_assay_strain_def_chal_average$assay_strain_defending_challenging <- factor(conf_assay_strain_def_chal_average$assay_strain_defending_challenging)



#Compute the range of scores for each unique assay
#First compute the minimum score and maximum score
#Then compute the range
conf_assay_strain_def_chal_min <- aggregate(conf_to_average[,c('scored_infection_strength',
                                                               'scored_turbidity',
                                                               'scored_plaque_size',
                                                               'scored_plaques',
                                                               'scored_four_factors',
                                                               'rank6')],
                                            list(conf_to_average$assay_strain_defending_challenging),min)


names(conf_assay_strain_def_chal_min) <- c('assay_strain_defending_challenging',
                                               'min_infection_strength',
                                               'min_turbidity',
                                               'min_plaque_size',
                                               'min_plaques',
                                               'min_four_factors',
                                               'min_rank6') 

conf_assay_strain_def_chal_min$assay_strain_defending_challenging <- factor(conf_assay_strain_def_chal_min$assay_strain_defending_challenging)


conf_assay_strain_def_chal_max <- aggregate(conf_to_average[,c('scored_infection_strength',
                                                               'scored_turbidity',
                                                               'scored_plaque_size',
                                                               'scored_plaques',
                                                               'scored_four_factors',
                                                               'rank6')],
                                            list(conf_to_average$assay_strain_defending_challenging),max)

names(conf_assay_strain_def_chal_max) <- c('assay_strain_defending_challenging',
                                           'max_infection_strength',
                                           'max_turbidity',
                                           'max_plaque_size',
                                           'max_plaques',
                                           'max_four_factors',
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


conf_assay_strain_def_chal_average$range_infection_strength <- conf_assay_strain_def_chal_average$max_infection_strength - 
  conf_assay_strain_def_chal_average$min_infection_strength

conf_assay_strain_def_chal_average$range_turbidity <- conf_assay_strain_def_chal_average$max_turbidity - 
  conf_assay_strain_def_chal_average$min_turbidity

conf_assay_strain_def_chal_average$range_plaque_size <- conf_assay_strain_def_chal_average$max_plaque_size - 
  conf_assay_strain_def_chal_average$min_plaque_size

conf_assay_strain_def_chal_average$range_plaques <- conf_assay_strain_def_chal_average$max_plaques - 
  conf_assay_strain_def_chal_average$min_plaques

conf_assay_strain_def_chal_average$range_four_factors <- conf_assay_strain_def_chal_average$max_four_factors - 
  conf_assay_strain_def_chal_average$min_four_factors

conf_assay_strain_def_chal_average$range_rank6 <- conf_assay_strain_def_chal_average$max_rank6 - 
  conf_assay_strain_def_chal_average$min_rank6

conf_assay_strain_def_chal_average$range_infection_strength <- as.factor(conf_assay_strain_def_chal_average$range_infection_strength)
conf_assay_strain_def_chal_average$range_turbidity <- as.factor(conf_assay_strain_def_chal_average$range_turbidity)
conf_assay_strain_def_chal_average$range_plaque_size <- as.factor(conf_assay_strain_def_chal_average$range_plaque_size)
conf_assay_strain_def_chal_average$range_plaques <- as.factor(conf_assay_strain_def_chal_average$range_plaques)
conf_assay_strain_def_chal_average$range_four_factors <- as.factor(conf_assay_strain_def_chal_average$range_four_factors)
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
                                            repressor_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            portal_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            recb_distance_data,
                                            by.x="defending_challenging",
                                            by.y="phage1_phage2",
                                            all.x=TRUE)
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
conf_assay_strain_def_chal_average <- merge(conf_assay_strain_def_chal_average,
                                            stoperator_tally,
                                            by.x="defending_challenging",
                                            by.y="tfbs88_motif_target",
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







#Export all averaged data
setwd("~/scratch/immunity_analysis/output/")
write.table(conf_assay_strain_def_chal_average,
            "conf_assay_strain_def_chal_average.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)






#Create average assay_strain_defending_challenging data above





















#Map averages back to main immunity table, subtract average from original values, then plot histogram of differences to show how reliable the dataset is
conf_assay_strain_def_chal_average_reduced <- subset(conf_assay_strain_def_chal_average,
                                                     select = c("assay_strain_defending_challenging",
                                                                "averaged_infection_strength",
                                                                "averaged_turbidity",
                                                                "averaged_plaque_size",
                                                                "averaged_plaques",
                                                                "averaged_four_factors",
                                                                "averaged_rank6",
                                                                "min_infection_strength",
                                                                "min_turbidity",
                                                                "min_plaque_size",
                                                                "min_plaques",
                                                                "min_four_factors",
                                                                "min_rank6",
                                                                "max_infection_strength",
                                                                "max_turbidity",
                                                                "max_plaque_size",
                                                                "max_plaques",
                                                                "max_four_factors",
                                                                "max_rank6",
                                                                "range_infection_strength",
                                                                "range_turbidity",
                                                                "range_plaque_size",
                                                                "range_plaques",
                                                                "range_four_factors",
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




main_immunity_data$assay_strain_averaged_infection_strength_diff <- as.numeric(as.character(main_immunity_data$scored_infection_strength)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_infection_strength))

main_immunity_data$assay_strain_averaged_turbidity_diff <- as.numeric(as.character(main_immunity_data$scored_turbidity)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_turbidity))

main_immunity_data$assay_strain_averaged_plaque_size_diff <- as.numeric(as.character(main_immunity_data$scored_plaque_size)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_plaque_size))

main_immunity_data$assay_strain_averaged_plaques_diff <- as.numeric(as.character(main_immunity_data$scored_plaques)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_plaques))

main_immunity_data$assay_strain_averaged_four_factors_diff <- as.numeric(as.character(main_immunity_data$scored_four_factors)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_four_factors))

main_immunity_data$assay_strain_averaged_rank6_diff <- as.numeric(as.character(main_immunity_data$rank6)) - 
  as.numeric(as.character(main_immunity_data$assay_strain_averaged_rank6))
#Note: if the diff is positive = assay's rank is higher than the average











# #TEMP
# 
# #Identify best candidates for figure of infection spectrum
# 
# temp_l5 <- subset(conf_assay_strain_def_chal_average,
#                   conf_assay_strain_def_chal_average$challenging_phage == 'l5' &
#                     conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
#                     conf_assay_strain_def_chal_average$strain_type == 'lysogen',
#                   select=c("defending_phage","challenging_phage","averaged_rank6","frequency","min_rank6","max_rank6"))
# 
# temp_gladiator <- subset(conf_assay_strain_def_chal_average,
#                          conf_assay_strain_def_chal_average$challenging_phage == 'gladiator' &
#                            conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
#                            conf_assay_strain_def_chal_average$strain_type == 'lysogen',
#                          select=c("defending_phage","challenging_phage","averaged_rank6","frequency","min_rank6","max_rank6"))
# 
# temp_darthphader <- subset(conf_assay_strain_def_chal_average,
#                            conf_assay_strain_def_chal_average$challenging_phage == 'darthphader' &
#                              conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
#                              conf_assay_strain_def_chal_average$strain_type == 'lysogen',
#                            select=c("defending_phage","challenging_phage","averaged_rank6","frequency","min_rank6","max_rank6"))
# 
# temp_archernm <- subset(conf_assay_strain_def_chal_average,
#                         conf_assay_strain_def_chal_average$challenging_phage == 'archernm' &
#                           conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
#                           conf_assay_strain_def_chal_average$strain_type == 'lysogen',
#                         select=c("defending_phage","challenging_phage","averaged_rank6","frequency","min_rank6","max_rank6"))
# 
# temp_et2brutus <- subset(conf_assay_strain_def_chal_average,
#                          conf_assay_strain_def_chal_average$challenging_phage == 'et2brutus' &
#                            conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
#                            conf_assay_strain_def_chal_average$strain_type == 'lysogen',
#                          select=c("defending_phage","challenging_phage","averaged_rank6","frequency","min_rank6","max_rank6"))
# 
# 
# #




#Assess reproducibility for all comparisons with frequency > 1
main_immunity_data_assay_strain_reduced <- subset(main_immunity_data,
                                                  as.numeric(as.character(main_immunity_data$assay_strain_frequency)) > 1)

# main_immunity_data_strain_reduced <- subset(main_immunity_data,
#                                             main_immunity_data$strain_frequency > 1)

















#Plots
setwd("~/scratch/immunity_analysis/output/")

conf_assay_strain_def_chal_average$frequency <- as.numeric(as.character(conf_assay_strain_def_chal_average$frequency))
conf_assay_strain_def_chal_average$frequency2 <- ifelse(conf_assay_strain_def_chal_average$frequency > 10,11,conf_assay_strain_def_chal_average$frequency)
conf_assay_strain_def_chal_average$frequency <- as.factor(conf_assay_strain_def_chal_average$frequency)
conf_assay_strain_def_chal_average$frequency2 <- as.factor(conf_assay_strain_def_chal_average$frequency2)


par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average$frequency),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,800))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_frequency.pdf')
dev.off()

par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average$frequency2),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,800))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_frequency_adjusted.pdf')
dev.off()



#Number of unique assays with >1 replicate. This includes all single-titer assays
1 - nrow(subset(conf_assay_strain_def_chal_average,conf_assay_strain_def_chal_average$frequency == "1"))/nrow(conf_assay_strain_def_chal_average)
#As of 1/17/2019, 48% of unique comparisons has > 1 replicate 


par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average$range_rank6),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,1200))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_rank6_range.pdf')
dev.off()
#As of 11/16/2018, 92% of unique comparisons with > 1 replicate have a score range of < 2.










#Only look at multiple_titer assays
conf_assay_strain_def_chal_average_multi <- subset(conf_assay_strain_def_chal_average,conf_assay_strain_def_chal_average$assay_type == "multiple_titer")

par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average_multi$frequency),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,400))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_multi_frequency.pdf')
dev.off()

par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average_multi$frequency2),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,400))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_multi_frequency_adjusted.pdf')
dev.off()




#Final Number of unique multi-titer assays
nrow(conf_assay_strain_def_chal_average_multi)

#Final Number of unique multi-titer assays with >1 replicate.
conf_assay_strain_def_chal_average_multi_reps <- subset(conf_assay_strain_def_chal_average_multi,conf_assay_strain_def_chal_average_multi$frequency != "1")
nrow(conf_assay_strain_def_chal_average_multi_reps)/nrow(conf_assay_strain_def_chal_average_multi)
#As of 1/17/2019, 67% of unique comparisons has > 1 replicate 



par(mar=c(4,8,8,4))
barplot(summary(conf_assay_strain_def_chal_average_multi_reps$range_rank6),col="black",ann=FALSE,main=NULL,las=1,ylim=c(0,400))
dev.copy(pdf,'conf_assay_strain_def_chal_ave_rank6_multi_reps_range.pdf')
dev.off()


#Final number of unique multi-titer assays with > 1 replicate and score range < 2
nrow(subset(conf_assay_strain_def_chal_average_multi_reps,conf_assay_strain_def_chal_average_multi_reps$range_rank6 == "0" | conf_assay_strain_def_chal_average_multi_reps$range_rank6 == "1"))/nrow(conf_assay_strain_def_chal_average_multi_reps)
#As of 1/17/2019, 82% of unique multi-titer (lysogen or rep clone) comparisons with > 1 replicate have a score range of < 2.










hist(main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff,
     xlim=c(-3,3),col='black',ann=FALSE,main=NULL,las=1,breaks=9)
dev.copy(pdf,'assay_strain_averaged_rank6_diff_hist.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(as.numeric(as.character(main_immunity_data_assay_strain_reduced$rank6)),
     main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff,
     pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(main_immunity_data_assay_strain_reduced$rank6)),
     main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff,
     pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(10,8,4,4))
boxplot(main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff ~ main_immunity_data_assay_strain_reduced$rank6,
        cex=1,cex.axis=1,ann=FALSE,main=NULL,las=2)


par(mar=c(10,8,4,4))
boxplot(as.numeric(as.character(main_immunity_data_assay_strain_reduced$assay_strain_range_rank6)) ~ main_immunity_data_assay_strain_reduced$rank6,
        cex=1,cex.axis=1,ann=FALSE,main=NULL,las=2)

par(mar=c(10,8,4,4))
boxplot(as.numeric(as.character(main_immunity_data_assay_strain_reduced$assay_strain_range_rank6)) ~ main_immunity_data_assay_strain_reduced$rank6,
        cex=1,cex.axis=1,ann=FALSE,main=NULL,las=2)

par(mar=c(16,8,4,4))
stripchart(main_immunity_data_assay_strain_reduced$assay_strain_min_rank6 ~ main_immunity_data_assay_strain_reduced$rank6,
           vertical=TRUE,las=2,cex.axis=1,pch=19,method="jitter",cex=0.5,ylab='')

par(mar=c(16,8,4,4))
stripchart(main_immunity_data_assay_strain_reduced$assay_strain_max_rank6 ~ main_immunity_data_assay_strain_reduced$rank6,
           vertical=TRUE,las=2,cex.axis=1,pch=19,method="jitter",cex=0.5,ylab='',col="blue")




par(mar=c(16,8,4,4))
stripchart(main_immunity_data_assay_strain_reduced$assay_strain_min_rank6 ~ main_immunity_data_assay_strain_reduced$rank6,
           ylim=c(0,6),vertical=TRUE,las=2,cex.axis=1,pch=19,method="jitter",cex=0.5,ylab='',col="green")
par(new=TRUE)
stripchart(main_immunity_data_assay_strain_reduced$assay_strain_max_rank6 ~ main_immunity_data_assay_strain_reduced$rank6,
           ylim=c(0,6),vertical=TRUE,las=2,cex.axis=1,pch=19,method="jitter",cex=0.5,ylab='',col="blue")

temp_range0 <- subset(main_immunity_data_assay_strain_reduced,main_immunity_data_assay_strain_reduced$assay_strain_range_rank6 == 0)
barplot(summary(temp_range0$rank6))

temp_range2 <- subset(main_immunity_data_assay_strain_reduced,main_immunity_data_assay_strain_reduced$assay_strain_range_rank6 == 2)
barplot(summary(temp_range2$rank6))







# hist(main_immunity_data_strain_reduced$strain_averaged_rank6_diff,
#      xlim=c(-3,3),col='black',ann=FALSE,main=NULL,las=1,breaks=101)
# dev.copy(pdf,'strain_averaged_rank6_diff_hist.pdf')
# dev.off()



#Export the immunity assays that are most divergent from their average
assay_outliers <- subset(main_immunity_data_assay_strain_reduced,
                         main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff < -1 |
                           main_immunity_data_assay_strain_reduced$assay_strain_averaged_rank6_diff > 1)
write.table(assay_outliers,
            "assay_outliers.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

#Conclusion: of all confident, averaged, assay_strain_def_chal data, there are 27 comparisons with a rank6 that exceeds the average by 1.
#The largest difference is 2.
#None of them appear to be egregious or unrealistic. 










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









#Immunity correlation analysis
setwd("~/scratch/immunity_analysis/output/")

# par(mar=c(4,8,4,4))
# plot(immunity_correlation_data$defending_cor_reduced,
#      immunity_correlation_data$challenging_cor_reduced,
#      xlim=c(-1,1),ylim=c(-1,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,1)
# dev.copy(pdf,'immunity_correlation_reduced.pdf')
# dev.off()
# 
# 
# lm_immunity_cor_reduced <- lm(defending_cor_reduced ~
#                                 challenging_cor_reduced,
#                               data = immunity_correlation_data)
# summary(lm_immunity_cor_reduced)









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



#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_gcd.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$modified_mash_distance,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$modified_mash_distance,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$modified_mash_distance,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_mash.pdf')
dev.off()
#




#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_stopEuc.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repCterm.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_full_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$repressor_full_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_full_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repFull.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$cas4_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$cas4_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$cas4_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_cas4.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$endovii_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$endovii_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$endovii_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_endovii.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$dnapol_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$dnapol_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$dnapol_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_dnapol.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$portal_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$portal_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$portal_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_portal.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_nterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$repressor_nterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_nterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repNterm.pdf')
dev.off()



#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_hth_compare,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,10),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$repressor_hth_compare,
     conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,10),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_hth_compare,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$averaged_rank6,
     xlim=c(0,10),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_subtypes_by_repHTH.pdf')
dev.off()




#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic) +
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic) +
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade)

nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade)





# #evaluate specific assays with similar hth
# temp_hth0 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_hth_compare < 1)
# temp_hth1 <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic,conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic$repressor_hth_compare == 1)


#FinalCorrelations
lm_gcd <- lm(averaged_rank6 ~
               repressor_full_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               cas4_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               endovii_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               dnapol_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               portal_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               stoperator_pwd_dist_euc,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               repressor_cterm_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)


lm_gcd <- lm(averaged_rank6 ~
               pham_pham_dissimilarity,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)


lm_gcd <- lm(averaged_rank6 ~
               modified_mash_distance,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               repressor_nterm_mafft_dist_uncorrected,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)

lm_gcd <- lm(averaged_rank6 ~
               repressor_hth_compare,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2_heterotypic)
summary(lm_gcd)















par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_gcd.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$pham_pham_dissimilarity,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$averaged_rank6)),
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_TempDiff_RepDiff_GCD_vs_ave_rank6.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$modified_mash_distance,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     xlim=c(0,0.5),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$modified_mash_distance,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$averaged_rank6)),
     xlim=c(0,0.5),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_GCD.pdf')
dev.off()

#Conclusions:
#These are the best examples of mutants or non-empirically-temperate phages infecting lysogens of closely-related phages
#L5_phitm38, L5_phiTM4, LadyBird_Piro94, StarStuf_phiTM38, Trixie_phiTM42, Gladiator_phiTM46
escape_examples <- subset(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff,
                          conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$pham_pham_dissimilarity < 0.37 &
                            conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$averaged_rank6 > 1.4)

#These are the best examples of prophages that complete defend against distantly-related phages
#It seems that the majority involve at least one extrachromosomal phage
defense_examples <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                           conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity < 0.62 &
                              conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity > 0.5 &
                             conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6 < 1.1)







#misc phages 


conf_assay_strain_ave_lys_multi_misc <- subset(conf_assay_strain_def_chal_average,
                                               conf_assay_strain_def_chal_average$strain_type == 'lysogen' &
                                                 conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                                                 conf_assay_strain_def_chal_average$defending_source == 'environment' &
                                                 (conf_assay_strain_def_chal_average$challenging_phage == "jeffabunny" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "misswhite" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "d29" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "phitm33" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "piro94" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "journey13" |
                                                    conf_assay_strain_def_chal_average$challenging_phage == "echild"))



conf_assay_strain_ave_lys_multi_misc_interclade <- subset(conf_assay_strain_ave_lys_multi_misc,
                                                          conf_assay_strain_ave_lys_multi_misc$gene_content_clade_compare == "different")

conf_assay_strain_ave_lys_multi_misc_intraclade2 <- subset(conf_assay_strain_ave_lys_multi_misc,
                                                           conf_assay_strain_ave_lys_multi_misc$gene_content_clade_compare == "clade2")





par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_misc$cas4_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_misc$averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_misc_cas4_vs_averank6.pdf')
dev.off()



par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_misc$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_misc$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'misc_phages_ave_rank6_vs_stopEuc.pdf')
dev.off()








par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_misc_intraclade2$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_misc_intraclade2$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_misc_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_misc_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'misc_phages_ave_rank6_vs_stopEuc_subtypes.pdf')
dev.off()


#






par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$modified_mash_distance,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_mash.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator_pwd_dist_euc,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_stop_euc.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator_pwd_dist_pearson,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator80_tally,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),xlim=c(0,90),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)



















#Compare WT phages and all mutant phages
par(mar=c(4,8,4,4))
plot(conf_assay_strain_def_chal_average$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_def_chal_average$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_TempDiff_RepDiff$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_TempDiff_RepDiff_repCtermMafUn_vs_ave_rank6.pdf')
dev.off()
#





par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_full_mafft_phyml_dist,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_repFullMafft.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_nterm_mafft_phyml_dist,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_repNterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_cterm_mafft_phyml_dist,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_repCterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_muscle_bionj_distances,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_repFullMuscle.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_hth_compare,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_repHTH.pdf')
dev.off()

#These are the best examples of strong defense despite much different hth sequences
hth_defense_examples <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                               conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_hth_compare == 10 &
                                 conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6 < 1.5)



par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$recb_muscle_bionj_distances,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_recb_vs_ave_rank6.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$cas4_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_by_cas4MafftUn_vs_ave_rank6.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$portal_muscle_bionj_distances,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_portal_vs_ave_rank6.pdf')
dev.off()






par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$portal_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_portalMafftUn_vs_ave_rank6.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$endovii_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_endoviiMafftUn_vs_ave_rank6.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$dnapol_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_dnapolMafftUn_vs_ave_rank6.pdf')
dev.off()






par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_gcd_vs_repMafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_gcd_vs_portalMafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_gcd_vs_endoviiMafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_gcd_vs_dnapolMafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$pham_pham_dissimilarity,
     conf_assay_strain_ave_lys_multi_env_temp_rep$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_gcd_vs_cas4MafftUn.pdf')
dev.off()









#Correlation with genes
#use pearson for continuous variables
# cor(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_full_mafft_dist_uncorrected,
#     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
#     use="complete.obs",
#     method="pearson")

lm_gcd <- lm(averaged_rank6 ~
               pham_pham_dissimilarity,
             data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_gcd)

lm_mash <- lm(averaged_rank6 ~
                modified_mash_distance,
              data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_mash)


lm_repFull <- lm(averaged_rank6 ~
                   repressor_full_mafft_dist_uncorrected,
                 data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_repFull)

lm_cas4 <- lm(averaged_rank6 ~
                cas4_mafft_dist_uncorrected,
              data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_cas4)

lm_endovii <- lm(averaged_rank6 ~
                   endovii_mafft_dist_uncorrected,
                 data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_endovii)

lm_dnapol <- lm(averaged_rank6 ~
                  dnapol_mafft_dist_uncorrected,
                data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_dnapol)

lm_portal <- lm(averaged_rank6 ~
                  portal_mafft_dist_uncorrected,
                data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_portal)

lm_stoppwm <- lm(averaged_rank6 ~
                   stoperator_pwd_dist_euc,
                 data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_stoppwm)





# plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_full_mafft_dist_uncorrected,
#      conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator_pwd_dist_euc,
#      xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
# abline(lm_repFull)
















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


















par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)


#IntInt
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_repCtermMafUn.pdf')
dev.off()

#ExtraExtra
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_repCtermMafUn.pdf')
dev.off()

#IntExtra
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_repCtermMafUn.pdf')
dev.off()

#ExtraInt
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_repCtermMafUn.pdf')
dev.off()



#IntInt_IntPhamSame
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntSame_repCtermMafUn.pdf')
dev.off()

#IntInt_IntPhamDiff
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntDiff_repCtermMafUn.pdf')
dev.off()


#ExtraExtra_ParBPhamSame
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraSame_repCtermMafUn.pdf')
dev.off()

#ExtraExtra_ParBPhamDiff
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraDiff_repCtermMafUn.pdf')
dev.off()






#
#IntExtra
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_stopEuc.pdf')
dev.off()

#ExtraInt
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_stopEuc.pdf')
dev.off()



#IntInt_IntPhamSame
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntSame_stopEuc.pdf')
dev.off()

#IntInt_IntPhamDiff
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntDiff_stopEuc.pdf')
dev.off()


#ExtraExtra_ParBPhamSame
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraSame_stopEuc.pdf')
dev.off()

#ExtraExtra_ParBPhamDiff
par(mar=c(4,8,8,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraDiff_stopEuc.pdf')
dev.off()



















#Revised Int/Extra by intra-Clade2 and interclade
#IntExtra
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_stopEuc_subtypes.pdf')
dev.off()


#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade)

nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntExtra_interclade)







#ExtraInt
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade)

nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraInt_interclade)






#IntInt_IntPhamSame
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntSame_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_intraclade2_heterotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamSame_interclade)





#IntInt_IntPhamDiff
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntIntDiff_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt_IntPhamDiff_interclade)






#ExtraExtra_ParBPhamSame
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraSame_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_homotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_intraclade2_heterotypic)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBSame_interclade)







#ExtraExtra_ParBPhamDiff
#FinalFig
par(mar=c(4,8,16,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade$stoperator_pwd_dist_euc,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtraDiff_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays:
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2) + 
  nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_intraclade2)
nrow(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra_ParBDiff_interclade)


#











#int
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntInt$averaged_rank6)),
     ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntSame$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntSame$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_IntDiff$repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(conf_assay_strain_ave_lys_multi_env_temp_rep_IntDiff$averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_IntDiffandSame_repCtermMafftUn_vs_aveRank6.pdf')
dev.off()



#partitioning
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraExtra$averaged_rank6,
     ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraDiff$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraDiff$averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraSame$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraSame$averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_ExtraDiffandSame_repCtermMafftUn_vs_aveRank6.pdf')
dev.off()
#











#Bin data.
#Choose either all data or just clade2
#data_to_bin <- conf_assay_strain_ave_lys_multi_env_temp_rep

data_to_bin <- conf_assay_strain_ave_lys_multi_env_temp_rep_intraclade2



#Choose a sub-dataset:
rank6_0 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 < 0.5)

rank6_1 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 0.5 &
                    data_to_bin$averaged_rank6 < 1.5)

rank6_2 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 1.5 &
                    data_to_bin$averaged_rank6 < 2.5)

rank6_3 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 2.5 &
                    data_to_bin$averaged_rank6 < 3.5)

rank6_4 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 3.5 &
                    data_to_bin$averaged_rank6 < 4.5)

rank6_5 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 4.5 &
                    data_to_bin$averaged_rank6 < 5.5)

rank6_6 <- subset(data_to_bin,
                  data_to_bin$averaged_rank6 >= 5.5)

#QC Check
nrow(data_to_bin) - 
  (nrow(rank6_0) + 
     nrow(rank6_1) + 
     nrow(rank6_2) + 
     nrow(rank6_3) + 
     nrow(rank6_4) + 
     nrow(rank6_5) + 
     nrow(rank6_6))




#



rank6_subset_columns <- c("repressor_cterm_mafft_dist_uncorrected",
                          "repressor_full_mafft_dist_uncorrected",
                          "stoperator_pwd_dist_euc",
                          "stoperator85_tally",
                          "portal_mafft_dist_uncorrected")
  
  
  
rank6_0_subset <- subset(rank6_0,select=rank6_subset_columns)
rank6_1_subset <- subset(rank6_1,select=rank6_subset_columns)
rank6_2_subset <- subset(rank6_2,select=rank6_subset_columns)
rank6_3_subset <- subset(rank6_3,select=rank6_subset_columns)
rank6_4_subset <- subset(rank6_4,select=rank6_subset_columns)
rank6_5_subset <- subset(rank6_5,select=rank6_subset_columns)
rank6_6_subset <- subset(rank6_6,select=rank6_subset_columns)
rank6_interclade_subset <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade,select=rank6_subset_columns)

rank6_0_subset$bin <- "bin0"
rank6_1_subset$bin <- "bin1"
rank6_2_subset$bin <- "bin2"
rank6_3_subset$bin <- "bin3"
rank6_4_subset$bin <- "bin4"
rank6_5_subset$bin <- "bin5"
rank6_6_subset$bin <- "bin6"
rank6_interclade_subset$bin <- "binInterClade"


rank6_binned_subset <- rbind(rank6_0_subset,
                      rank6_1_subset,
                      rank6_2_subset,
                      rank6_3_subset,
                      rank6_4_subset,
                      rank6_5_subset,
                      rank6_6_subset,
                      rank6_interclade_subset)



rank6_binned_subset$bin <- factor(rank6_binned_subset$bin,
                           c("binInterClade",
                             "bin6",
                             "bin5",
                             "bin4",
                             "bin3",
                             "bin2",
                             "bin1",
                             "bin0"))





plot_colors <- c("black",
                 "dark green",
                 "green",
                 "light green",
                 "light grey",
                 "grey",
                 "dark grey",
                 "black")

errorCircles("repressor_cterm_mafft_dist_uncorrected",
             "stoperator85_tally",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,70),
             ylim=c(0,40),
             colors=plot_colors)
dev.copy(pdf,'averaged_repCtermMafftUn_vs_stop85tally.pdf')
dev.off()



errorCircles("repressor_cterm_mafft_dist_uncorrected",
             "stoperator_pwd_dist_euc",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,70),
             ylim=c(0,5),
             colors=plot_colors)
dev.copy(pdf,'averaged_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()






errorCircles("repressor_full_mafft_dist_uncorrected",
             "stoperator85_tally",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,70),
             ylim=c(0,40),
             colors=plot_colors)


errorCircles("repressor_full_mafft_dist_uncorrected",
             "stoperator_pwd_dist_euc",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,70),
             ylim=c(0,5),
             colors=plot_colors)


errorCircles("repressor_full_mafft_dist_uncorrected",
             "stoperator_pwd_dist_euc",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,45),
             ylim=c(0,3),
             colors=plot_colors)

errorCircles("portal_mafft_dist_uncorrected",
             "stoperator_pwd_dist_euc",
             data=rank6_binned_subset,
             group="bin",
             circles=FALSE,
             sd=TRUE,
             pch=16,
             cex=2,
             labels="",
             xlim=c(0,25),
             ylim=c(0,3),
             colors=plot_colors)
#













plot(data_to_bin$modified_mash_distance,
     data_to_bin$repressor_full_mafft_phyml_dist,xlim=c(0,0.5),ylim=c(0,2.5))
     


par(mar=c(4,8,4,4))
plot(rank6_0$modified_mash_distance,
     rank6_0$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_1$modified_mash_distance,
     rank6_1$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_2$modified_mash_distance,
     rank6_2$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_3$modified_mash_distance,
     rank6_3$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_4$modified_mash_distance,
     rank6_4$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_5$modified_mash_distance,
     rank6_5$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,0.5),ylim=c(0,2.5))
par(new=TRUE)
plot(rank6_6$modified_mash_distance,
     rank6_6$repressor_full_mafft_phyml_dist,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,0.5),ylim=c(0,2.5))




plot(data_to_bin$repressor_nterm_mafft_dist_uncorrected,
     data_to_bin$repressor_cterm_mafft_dist_uncorrected,xlim=c(0,70),ylim=c(0,70))


par(mar=c(4,8,4,4))
plot(rank6_0$repressor_nterm_mafft_dist_uncorrected,
     rank6_0$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_1$repressor_nterm_mafft_dist_uncorrected,
     rank6_1$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_2$repressor_nterm_mafft_dist_uncorrected,
     rank6_2$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_3$repressor_nterm_mafft_dist_uncorrected,
     rank6_3$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_4$repressor_nterm_mafft_dist_uncorrected,
     rank6_4$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_5$repressor_nterm_mafft_dist_uncorrected,
     rank6_5$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,70))
par(new=TRUE)
plot(rank6_6$repressor_nterm_mafft_dist_uncorrected,
     rank6_6$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,70))



plot(data_to_bin$repressor_hth_compare,
     data_to_bin$repressor_cterm_mafft_dist_uncorrected,xlim=c(0,10),ylim=c(0,70))

par(mar=c(4,8,4,4))
plot(rank6_0$repressor_hth_compare,
     rank6_0$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_1$repressor_hth_compare,
     rank6_1$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_2$repressor_hth_compare,
     rank6_2$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_3$repressor_hth_compare,
     rank6_3$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_4$repressor_hth_compare,
     rank6_4$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_5$repressor_hth_compare,
     rank6_5$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,10),ylim=c(0,70))
par(new=TRUE)
plot(rank6_6$repressor_hth_compare,
     rank6_6$repressor_cterm_mafft_dist_uncorrected,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,10),ylim=c(0,70))






#The multi-colored plot below is an alternative to the same plot further below that uses a color-spectrum
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator_pwd_dist_euc,xlim=c(0,70),ylim=c(0,5))


par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,5))
par(new=TRUE)
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,5))




#Now split up the plots individually
par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_0_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_1_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_2_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_3_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_4_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_5_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_6_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()
#






#Same as above, but no color
#Now split up the plots individually
par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_0_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_1_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_2_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_3_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_4_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_5_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'rank6_6_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()




#inter-clade data here
par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'interclade_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()



#































#stoperator_tally85
plot(conf_assay_strain_ave_lys_multi_env_temp_rep$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep$stoperator85_tally,xlim=c(0,70),ylim=c(0,40))



par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,40))
par(new=TRUE)
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_binned_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()






#Split into individual plots
par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_0_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_1_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_2_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_3_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_4_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_5_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_6_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()




#





#Same as above, but no color
#Split into individual plots


par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_0_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_1_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_2_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_3_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_4_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_5_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'rank6_6_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_interclade$stoperator85_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,40))
dev.copy(pdf,'interclade_repCtermMafftUn_vs_stop85_tally.pdf')
dev.off()





#































lm_r6_cterm_stoppwm <- lm(averaged_rank6 ~
                            repressor_cterm_mafft_dist_uncorrected + stoperator_pwd_dist_euc,
                          data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_r6_cterm_stoppwm)





lm_r6_cterm_stop85 <- lm(averaged_rank6 ~
                            repressor_cterm_mafft_dist_uncorrected + stoperator85_tally,
                          data = conf_assay_strain_ave_lys_multi_env_temp_rep)
summary(lm_r6_cterm_stop85)








#stop85_up
#this data looks good
par(mar=c(4,8,4,4))
plot(rank6_0$repressor_cterm_mafft_dist_uncorrected,
     rank6_0$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_1$repressor_cterm_mafft_dist_uncorrected,
     rank6_1$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_2$repressor_cterm_mafft_dist_uncorrected,
     rank6_2$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_3$repressor_cterm_mafft_dist_uncorrected,
     rank6_3$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_4$repressor_cterm_mafft_dist_uncorrected,
     rank6_4$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_5$repressor_cterm_mafft_dist_uncorrected,
     rank6_5$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,70),ylim=c(0,6))
par(new=TRUE)
plot(rank6_6$repressor_cterm_mafft_dist_uncorrected,
     rank6_6$stoperator85_up_tally,
     pch=16,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,70),ylim=c(0,6))
#













































conf_assay_strain_ave_lys_multi_env_temp_rep_temp <- conf_assay_strain_ave_lys_multi_env_temp_rep
conf_assay_strain_ave_lys_multi_env_temp_rep_temp$averaged_rank6_scaled <- conf_assay_strain_ave_lys_multi_env_temp_rep_temp$averaged_rank6/6

conf_assay_strain_ave_lys_multi_env_temp_rep_reduced <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep,
                                                               select=(c("repressor_cterm_mafft_dist_uncorrected",
                                                                         "stoperator_pwd_dist_euc",
                                                                         "averaged_rank6")))

# conf_assay_strain_ave_lys_multi_env_temp_rep_reduced <- conf_assay_strain_ave_lys_multi_env_temp_rep_reduced[complete.cases(conf_assay_strain_ave_lys_multi_env_temp_rep_reduced),]
# plot(conf_assay_strain_ave_lys_multi_env_temp_rep_reduced$repressor_cterm_mafft_dist_uncorrected,
#      conf_assay_strain_ave_lys_multi_env_temp_rep_reduced$stoperator_pwd_dist_euc,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
#      col=rgb((colorRamp(c("black","blue","cyan","yellow","orange","red"))(conf_assay_strain_ave_lys_multi_env_temp_rep_reduced$averaged_rank6_scaled))/255),
#      xlim=c(0,70),ylim=c(0,5))

plot(conf_assay_strain_ave_lys_multi_env_temp_rep_temp$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_temp$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
     col=rgb((colorRamp(c("black","blue","cyan","yellow","orange","red"))(conf_assay_strain_ave_lys_multi_env_temp_rep_temp$averaged_rank6_scaled))/255),
     xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_env_temp_rep_repCtermMafftUn_vs_stopEuc_aveRank6Col.pdf')
dev.off()



conf_assay_strain_ave_lys_multi_env_temp_rep_temp_gcd <- subset(conf_assay_strain_ave_lys_multi_env_temp_rep_temp,
                                                                conf_assay_strain_ave_lys_multi_env_temp_rep_temp$pham_pham_dissimilarity > 0.4 &
                                                                  conf_assay_strain_ave_lys_multi_env_temp_rep_temp$pham_pham_dissimilarity < 0.7)


plot(conf_assay_strain_ave_lys_multi_env_temp_rep_temp_gcd$repressor_cterm_mafft_dist_uncorrected,
     conf_assay_strain_ave_lys_multi_env_temp_rep_temp_gcd$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
     col=rgb((colorRamp(c("black","blue","cyan","yellow","orange","red"))(conf_assay_strain_ave_lys_multi_env_temp_rep_temp_gcd$averaged_rank6_scaled))/255),
     xlim=c(0,70),ylim=c(0,5))












#Updated Binned Data Stats

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



par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$total_num_assays,
        names.arg=clade2_binned_frequency$bin,
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_freq_total.pdf')
dev.off()

par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$defending_freq,
        names.arg=clade2_binned_frequency$bin,
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_freq_defending.pdf')
dev.off()



par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$challenging_freq,
        names.arg=clade2_binned_frequency$bin,
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_freq_challenging.pdf')
dev.off()

par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$num_intra_subcluster_assays,
        names.arg=clade2_binned_frequency$bin,
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_freq_intrasubcluster.pdf')
dev.off()



par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$num_inter_subcluster_assays,
        names.arg=clade2_binned_frequency$bin,
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_freq_intersubcluster.pdf')
dev.off()









#percent of total

#FinalFig
par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$total_assays_percent,
        names.arg=clade2_binned_frequency$bin,
        col='black',ylim=c(0,0.3))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_total.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$defending_percent,
        names.arg=clade2_binned_frequency$bin,
        col='black',ylim=c(0,1))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_defending.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$challenging_percent,
        names.arg=clade2_binned_frequency$bin,
        col='black',ylim=c(0,1))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_challenging.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$intra_subcluster_percent,
        names.arg=clade2_binned_frequency$bin,
        col='black',ylim=c(0,0.50))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_intrasubcluster.pdf')
dev.off()

#FinalFig
par(mar=c(4,8,4,4))
barplot(clade2_binned_frequency$inter_subcluster_percent,
        names.arg=clade2_binned_frequency$bin,
        col='black',ylim=c(0,0.50))
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_clade2_env_binned_percent_intersubcluster.pdf')
dev.off()
#Clade 2 data above







#Summarize diversity of phenotypes



temp <- conf_assay_strain_ave_lys_multi

temp$averaged_rank6 <- as.factor(temp$averaged_rank6)


par(mar=c(4,8,4,4))
barplot(summary(as.factor(conf_assay_strain_ave_lys_multi$averaged_rank6)),
        col='black')
dev.copy(pdf,'conf_assay_strain_ave_lys_multi_binned_freq_defending.pdf')
dev.off()
#TODO: in progress
#





###Above: averaged data






























###Reciprocal analysis
#Using averaged data, compute differences in reciprocal immunity data


#conf_assay_strain_def_chal_average column names:
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
# "prophage"
# "repressor_clone"
# "strain_type"
# "assay_type"
# "defending_host"
# "defending_cluster"
# "defending_subcluster"
# "defending_size"
# "defending_status"
# "defending_author"
# "defending_datelastmodified"
# "defending_mode_approx_80_percent"
# "defending_network005_interaction_tally"
# "defending_network005_group"
# "defending_network005_group_tally"
# "defending_lysogen_type"
# "defending_pham_integrase"
# "defending_pham_para"
# "defending_source"
# "defending_parent"
# "defending_cluster_a_functional_repressor_predicted"
# "defending_cluster_a_temperate_empirical"
# "defending_maxgcdgap_all"
# "defending_repressor_hth_domain_sequence"
# "defending_repressor_length_full"
# "defending_repressor_length_nterm"
# "defending_repressor_length_cterm"
# "defending_pham_parb"
# "challenging_host"
# "challenging_cluster"
# "challenging_subcluster"
# "challenging_size"
# "challenging_status"
# "challenging_author"
# "challenging_datelastmodified"
# "challenging_mode_approx_80_percent"
# "challenging_network005_interaction_tally"
# "challenging_network005_group"
# "challenging_network005_group_tally"
# "challenging_lysogen_type"
# "challenging_pham_integrase"
# "challenging_pham_para"
# "challenging_source"
# "challenging_parent"
# "challenging_cluster_a_functional_repressor_predicted"
# "challenging_cluster_a_temperate_empirical"
# "challenging_maxgcdgap_all"
# "challenging_repressor_hth_domain_sequence"
# "challenging_repressor_length_full"
# "challenging_repressor_length_nterm"
# "challenging_repressor_length_cterm"
# "challenging_pham_parb"
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
# "parb_compare"
# "repressor_hth_compare"
# "repressor_length_full_compare"
# "repressor_length_nterm_compare"
# "repressor_length_cterm_compare"
# "frequency"








#Split conf_assay_strain_def_chal_average columns into groups



#Phage-specific data = metadata that is not impacted by immunity vector
# "prophage"
# "repressor_clone"
# "strain_type"
# "defending_host"
# "defending_cluster"
# "defending_subcluster"
# "defending_size"
# "defending_status"
# "defending_author"
# "defending_datelastmodified"
# "defending_mode_approx_80_percent"
# "defending_network005_interaction_tally"
# "defending_network005_group"
# "defending_network005_group_tally"
# "defending_lysogen_type"
# "defending_pham_integrase"
# "defending_pham_para"
# "defending_source"
# "defending_parent"
# "defending_cluster_a_functional_repressor_predicted"
# "defending_cluster_a_temperate_empirical"
# "defending_maxgcdgap_all"
# "defending_repressor_hth_domain_sequence"
# "defending_repressor_length_full"
# "defending_repressor_length_nterm"
# "defending_repressor_length_cterm"
# "defending_pham_parb"
# "challenging_host"
# "challenging_cluster"
# "challenging_subcluster"
# "challenging_size"
# "challenging_status"
# "challenging_author"
# "challenging_datelastmodified"
# "challenging_mode_approx_80_percent"
# "challenging_network005_interaction_tally"
# "challenging_network005_group"
# "challenging_network005_group_tally"
# "challenging_lysogen_type"
# "challenging_pham_integrase"
# "challenging_pham_para"
# "challenging_source"
# "challenging_parent"
# "challenging_cluster_a_functional_repressor_predicted"
# "challenging_cluster_a_temperate_empirical"
# "challenging_maxgcdgap_all"
# "challenging_repressor_hth_domain_sequence"
# "challenging_repressor_length_full"
# "challenging_repressor_length_nterm"
# "challenging_repressor_length_cterm"
# "challenging_pham_parb"






#Phage metadata comparisons = data specific to both phages used in immunity but not impacted by vector orientation
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
# "stoperator88_tally"                                  
# "stoperator88_2kbright_tally"                         
# "stoperator88_2kbleft_tally"                          
# "stoperator88_2kbright_targetself_tally"              
# "stoperator88_2kbright_nontargetself_tally"           
# "stoperator88_2kbleft_targetself_tally"               
# "stoperator88_2kbleft_nontargetself_tally"            
# "tfbs88_stoperator_target"                            
# "tfbs88_stoperator_motif"                             





vectored_column_names <- c("assay_strain_defending_challenging",
                           "defending_challenging",
                           "defending_phage",
                           "challenging_phage",
                           "averaged_infection_strength",
                           "averaged_turbidity",
                           "averaged_plaque_size",
                           "averaged_plaques",
                           "averaged_four_factors",
                           "averaged_rank6",
                           "strain_defending_challenging",
                           "assay_type",
                           "frequency",
                           "stoperator88_tally",
                           "stoperator88_2kbright_tally",
                           "stoperator88_2kbleft_tally",
                           "stoperator88_2kbright_targetself_tally",
                           "stoperator88_2kbright_nontargetself_tally",
                           "stoperator88_2kbleft_targetself_tally",
                           "stoperator88_2kbleft_nontargetself_tally",
                           "tfbs88_stoperator_target",
                           "tfbs88_stoperator_motif"
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
reciprocal_data$averaged_infection_strength_diff <- abs(reciprocal_data$vector1_averaged_infection_strength - reciprocal_data$vector2_averaged_infection_strength)
reciprocal_data$averaged_four_factors_diff <- abs(reciprocal_data$vector1_averaged_four_factors - reciprocal_data$vector2_averaged_four_factors)
reciprocal_data$averaged_rank6_diff <- abs(reciprocal_data$vector1_averaged_rank6 - reciprocal_data$vector2_averaged_rank6)






#Now that all data has been matched, the data is duplicated. (vector1 will contain l5_alma and alma_l5, so both vector2's will be present as well)
#To solve this, remove all vector1 rows in which the defending_phage and challenging_phage are not alphabetically ordered.
#This does NOT retain self-comparisons (e.g. alma_alma)
reciprocal_data$vector1_alpha_ordered <- as.character(reciprocal_data$vector1_defending_phage) < as.character(reciprocal_data$vector1_challenging_phage)



#This also accounts for self-comparisons (e.g. alma_alma)
# reciprocal_data$vector1_alpha_ordered <- (as.character(reciprocal_data$vector1_defending_phage) < as.character(reciprocal_data$vector1_challenging_phage)) |
#                                          (as.character(reciprocal_data$vector1_defending_phage) == as.character(reciprocal_data$vector1_challenging_phage))





#Now retain only the unique pairwise comparisons
reciprocal_data_alpha_ordered <- subset(reciprocal_data,reciprocal_data$vector1_alpha_ordered == TRUE)

#This is the inverse dataset
#reciprocal_data_alpha_reversed <- subset(reciprocal_data,reciprocal_data$vector1_alpha_ordered == FALSE)






#Output the reciprocal dataset
setwd("~/scratch/immunity_analysis/output/")
write.table(reciprocal_data_alpha_ordered,
            "reciprocal_data_alpha_ordered.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)



#Plot data


plot(as.numeric(as.character(reciprocal_data_alpha_ordered$vector1_averaged_rank6)),
     as.numeric(as.character(reciprocal_data_alpha_ordered$vector2_averaged_rank6)),
     pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_nonself_ordered_ave_rank6.pdf')
dev.off()






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


setwd("~/scratch/immunity_analysis/output/")




#rank6 plots
plot(as.numeric(as.character(reciprocal_unique_envY$vector1_averaged_rank6)),
     as.numeric(as.character(reciprocal_unique_envY$vector2_averaged_rank6)),
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_envY_ave_rank6.pdf')
dev.off()


plot(as.numeric(as.character(reciprocal_unique_envY_intY$vector1_averaged_rank6)),
     as.numeric(as.character(reciprocal_unique_envY_intY$vector2_averaged_rank6)),
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_envY_intY_ave_rank6.pdf')
dev.off()

plot(as.numeric(as.character(reciprocal_unique_envY_extraY$vector1_averaged_rank6)),
     as.numeric(as.character(reciprocal_unique_envY_extraY$vector2_averaged_rank6)),
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_envY_extraY_ave_rank6.pdf')
dev.off()









hist(reciprocal_unique_envY_intY$averaged_four_factors_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=8)
dev.copy(pdf,'reciprocal_unique_nonself_envY_intY_ave_four_factors_diff_hist.pdf')
dev.off()

hist(reciprocal_unique_envY_extraY$averaged_four_factors_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=8)
dev.copy(pdf,'reciprocal_unique_nonself_envY_extraY_ave_four_factors_diff_hist.pdf')
dev.off()



hist(reciprocal_unique_envY_intY$averaged_infection_strength_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10)
dev.copy(pdf,'reciprocal_unique_nonself_envY_intY_ave_inf_strength_diff_hist.pdf')
dev.off()


hist(reciprocal_unique_envY_extraY$averaged_infection_strength_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10)
dev.copy(pdf,'reciprocal_unique_nonself_envY_extraY_ave_inf_strength_diff_hist.pdf')
dev.off()





#rank6 hist
hist(reciprocal_unique_envY$averaged_rank6_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10,xlim=c(0,6))
dev.copy(pdf,'reciprocal_unique_envY_rank6_diff_hist.pdf')
dev.off()

hist(reciprocal_unique_envY_intY$averaged_rank6_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10,xlim=c(0,6))
dev.copy(pdf,'reciprocal_unique_envY_intY_rank6_diff_hist.pdf')
dev.off()

hist(reciprocal_unique_envY_extraY$averaged_rank6_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10,xlim=c(0,6))
dev.copy(pdf,'reciprocal_unique_envY_extraY_rank6_diff_hist.pdf')
dev.off()







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


par(mar=c(4,8,4,4))
barplot(reciprocal_binned_freq$freq,
        names.arg=reciprocal_binned_freq$bin,
        col='black')
dev.copy(pdf,'reciprocal_unique_env_binned_freq_total.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,4,4))
barplot(reciprocal_binned_freq$freq_percent,
        names.arg=reciprocal_binned_freq$bin,
        col='black',ylim=c(0,0.5))
dev.copy(pdf,'reciprocal_unique_env_binned_freq_percent.pdf')
dev.off()

#Binning above










#How does various genome/gene distances correlate with differences in reciprocal profiles?



#rank6

#reciprocal_unique_envY: subsetted from conf, assay_strain_ave, lys, multi, env data. Since it is reciprocal, it is by definition tempY and repY data.
reciprocal_unique_envY_interclade <- subset(reciprocal_unique_envY,
                                            reciprocal_unique_envY$gene_content_clade_compare == "different")

reciprocal_unique_envY_intraclade2 <- subset(reciprocal_unique_envY,
                                             reciprocal_unique_envY$gene_content_clade_compare == "clade2")

#There are no reciprocal homotypic data, since this was intentionally removed above.

reciprocal_unique_envY_intraclade2_heterotypic <- subset(reciprocal_unique_envY_intraclade2,
                                                         as.character(reciprocal_unique_envY_intraclade2$vector1_defending_phage) != 
                                                           as.character(reciprocal_unique_envY_intraclade2$vector1_challenging_phage))






#split by clade and homo/heterotypic

par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY_interclade$pham_pham_dissimilarity,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(reciprocal_unique_envY_intraclade2_heterotypic$pham_pham_dissimilarity,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_gcd.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY_interclade$modified_mash_distance,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(reciprocal_unique_envY_intraclade2_heterotypic$modified_mash_distance,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_mash.pdf')
dev.off()
#








#

#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$modified_mash_distance,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$modified_mash_distance,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,0.5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_mash.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$pham_pham_dissimilarity,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$pham_pham_dissimilarity,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,1),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_gcd.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$portal_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$portal_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_portal.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$dnapol_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$dnapol_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_dnapol.pdf')
dev.off()



#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$endovii_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$endovii_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_endovii.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$cas4_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$cas4_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_cas4.pdf')
dev.off()



#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$repressor_full_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$repressor_full_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_repFull.pdf')
dev.off()



#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_repCterm.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$repressor_nterm_mafft_dist_uncorrected,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$repressor_nterm_mafft_dist_uncorrected,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,70),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_repNterm.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$repressor_hth_compare,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,10),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$repressor_hth_compare,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,10),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_repHTH.pdf')
dev.off()



#FinalFig
par(mar=c(4,8,16,4))
plot(reciprocal_unique_envY_intraclade2_heterotypic$stoperator_pwd_dist_euc,
     reciprocal_unique_envY_intraclade2_heterotypic$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(reciprocal_unique_envY_interclade$stoperator_pwd_dist_euc,
     reciprocal_unique_envY_interclade$averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'reciprocal_unique_envY_subtypes_by_stopEuc.pdf')
dev.off()


#Final Number of assays used:
nrow(reciprocal_unique_envY_intraclade2_heterotypic) +
  nrow(reciprocal_unique_envY_interclade)
nrow(reciprocal_unique_envY_intraclade2_heterotypic)
nrow(reciprocal_unique_envY_interclade)

#











par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$modified_mash_distance)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_mash_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$pham_pham_dissimilarity)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_GCD_vs_ave_rank6_diff.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_full_mafft_phyml_dist)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_mash_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_nterm_mafft_phyml_dist)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_mash_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_cterm_mafft_phyml_dist)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_mash_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_full_mafft_dist_uncorrected)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_repFullMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_nterm_mafft_dist_uncorrected)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_repNtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$repressor_hth_compare)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_repHth_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$recb_muscle_bionj_distances)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_recbMuscle_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(as.numeric(as.character(reciprocal_unique_envY$stoperator_pwd_dist_euc)),
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_stopEuc_vs_ave_rank6_diff.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY$cas4_mafft_dist_uncorrected,
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_cas4MafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY$endovii_mafft_dist_uncorrected,
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_endoMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY$dnapol_mafft_dist_uncorrected,
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_dnapolMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(reciprocal_unique_envY$portal_mafft_dist_uncorrected,
     as.numeric(as.character(reciprocal_unique_envY$averaged_rank6_diff)),
     pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'reciprocal_unique_nonself_envY_portalMafftUn_vs_ave_rank6_diff.pdf')
dev.off()





#Plot with smoothed/averaged line
scatter.smooth(reciprocal_unique_envY$portal_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_portalMafftUn_vs_ave_rank6_diff.pdf')
dev.off()


scatter.smooth(reciprocal_unique_envY$dnapol_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_dnapolMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$endovii_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_endoviiMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$cas4_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_cas4MafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$repressor_full_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_repFullMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$repressor_nterm_mafft_dist_uncorrected,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_repNtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$repressor_hth_compare,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_repHtH_vs_ave_rank6_diff.pdf')
dev.off()

scatter.smooth(reciprocal_unique_envY$stoperator_pwd_dist_euc,
               reciprocal_unique_envY$averaged_rank6_diff,
               pch=1,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(0,6))
dev.copy(pdf,'reciprocal_unique_nonself_envY_stopEuc_vs_ave_rank6_diff.pdf')
dev.off()
#











#TODO: plot reciprocal by repressor and PWD dist, binned


# temp <- subset(reciprocal_unique_envY,reciprocal_unique_envY$vector1_defending_challenging == 'gladiator_trixie')
# plot(temp$repressor_cterm_mafft_dist_uncorrected,
#      temp$stoperator_pwd_dist_euc,xlim=c(0,60),ylim=c(0,4))

reciprocal_unique_envY$averaged_rank6_scaled <- reciprocal_unique_envY$averaged_rank6/6

plot(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
     col=rgb((colorRamp(c("black","blue","cyan","green","orange","red"))(reciprocal_unique_envY$averaged_rank6_scaled))/255))

plot(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,
     col=rgb((colorRamp(c("black","blue","green","orange","yellow"))(reciprocal_unique_envY$averaged_rank6_scaled))/255))


#Choose a sub-dataset:
#all data or only clade2
#recip_data_to_bin <- reciprocal_unique_envY

reciprocal_unique_envY_intraclade2 <- subset(reciprocal_unique_envY,
                                             reciprocal_unique_envY$gene_content_clade_compare == "clade2")

recip_data_to_bin <- reciprocal_unique_envY_intraclade2

recip_rank6diff_0 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff < 0.5)

recip_rank6diff_1 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff >= 0.5 &
                              recip_data_to_bin$averaged_rank6_diff < 1.5)

recip_rank6diff_2 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6 >= 1.5 &
                              recip_data_to_bin$averaged_rank6_diff < 2.5)

recip_rank6diff_3 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff >= 2.5 &
                              recip_data_to_bin$averaged_rank6_diff < 3.5)

recip_rank6diff_4 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff >= 3.5 &
                              recip_data_to_bin$averaged_rank6_diff < 4.5)

recip_rank6diff_5 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff >= 4.5 &
                              recip_data_to_bin$averaged_rank6_diff < 5.5)

recip_rank6diff_6 <- subset(recip_data_to_bin,
                            recip_data_to_bin$averaged_rank6_diff >= 5.5)








#Color spectrum of infection change by stoperator tally
plot(reciprocal_unique_envY$vector1_stoperator88_tally,
     reciprocal_unique_envY$vector2_stoperator88_tally,
     xlim=c(0,40),ylim=c(0,40))



par(mar=c(4,8,4,4))
plot(recip_rank6diff_0$vector1_stoperator88_tally,
     recip_rank6diff_0$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_1$vector1_stoperator88_tally,
     recip_rank6diff_1$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_2$vector1_stoperator88_tally,
     recip_rank6diff_2$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_3$vector1_stoperator88_tally,
     recip_rank6diff_3$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_4$vector1_stoperator88_tally,
     recip_rank6diff_4$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_5$vector1_stoperator88_tally,
     recip_rank6diff_5$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,40),ylim=c(0,40))
par(new=TRUE)
plot(recip_rank6diff_6$vector1_stoperator88_tally,
     recip_rank6diff_6$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,40),ylim=c(0,40))
abline(0,1)



















#Stop motif versus 

plot(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,60),ylim=c(0,4))








par(mar=c(4,8,4,4))
plot(recip_rank6diff_0$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_1$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_2$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_3$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_4$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_5$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,60),ylim=c(0,4))
par(new=TRUE)
plot(recip_rank6diff_6$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,60),ylim=c(0,4))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc.pdf')
dev.off()







plot(reciprocal_unique_envY$repressor_cterm_mafft_dist_uncorrected,
     reciprocal_unique_envY$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,60),ylim=c(0,5))

plot(recip_data_to_bin$repressor_cterm_mafft_dist_uncorrected,
     recip_data_to_bin$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,60),ylim=c(0,5))





plot(recip_rank6diff_0$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank0.pdf')
dev.off()


plot(recip_rank6diff_1$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank1.pdf')
dev.off()


plot(recip_rank6diff_2$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='cyan',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank2.pdf')
dev.off()


plot(recip_rank6diff_3$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank3.pdf')
dev.off()


plot(recip_rank6diff_4$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='yellow',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank4.pdf')
dev.off()


plot(recip_rank6diff_5$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank5.pdf')
dev.off()


plot(recip_rank6diff_6$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red',xlim=c(0,60),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank6.pdf')
dev.off()








#same plots with no color
setwd("~/scratch/immunity_analysis/output/")


plot(recip_rank6diff_0$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_0$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank0.pdf')
dev.off()


plot(recip_rank6diff_1$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_1$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank1.pdf')
dev.off()


plot(recip_rank6diff_2$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_2$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank2.pdf')
dev.off()


plot(recip_rank6diff_3$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_3$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank3.pdf')
dev.off()


plot(recip_rank6diff_4$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_4$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank4.pdf')
dev.off()


plot(recip_rank6diff_5$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_5$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank5.pdf')
dev.off()


plot(recip_rank6diff_6$repressor_cterm_mafft_dist_uncorrected,
     recip_rank6diff_6$stoperator_pwd_dist_euc,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,70),ylim=c(0,5))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_repCtermMafftUn_vs_stopEuc_rank6.pdf')
dev.off()







#reciprocal # stop sites


plot(reciprocal_unique_envY$vector1_stoperator88_tally,
     reciprocal_unique_envY$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,60),ylim=c(0,5))


plot(recip_data_to_bin$vector1_stoperator88_tally,
     recip_data_to_bin$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,35),ylim=c(0,35))
abline(0,1)



plot(recip_rank6diff_0$vector1_stoperator88_tally,
     recip_rank6diff_0$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank0.pdf')
dev.off()


plot(recip_rank6diff_1$vector1_stoperator88_tally,
     recip_rank6diff_1$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank1.pdf')
dev.off()


plot(recip_rank6diff_2$vector1_stoperator88_tally,
     recip_rank6diff_2$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank2.pdf')
dev.off()


plot(recip_rank6diff_3$vector1_stoperator88_tally,
     recip_rank6diff_3$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank3.pdf')
dev.off()


plot(recip_rank6diff_4$vector1_stoperator88_tally,
     recip_rank6diff_4$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank4.pdf')
dev.off()


plot(recip_rank6diff_5$vector1_stoperator88_tally,
     recip_rank6diff_5$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank5.pdf')
dev.off()


plot(recip_rank6diff_6$vector1_stoperator88_tally,
     recip_rank6diff_6$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,35),ylim=c(0,35))
dev.copy(pdf,'recip_unique_nonself_envY_rank6_diff_stop85_vs_stop85_rank6.pdf')
dev.off()








#2kb right stop sites
plot(reciprocal_unique_envY$vector1_stoperator88_2kbright_tally,
     reciprocal_unique_envY$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,60),ylim=c(0,5))


plot(recip_data_to_bin$vector1_stoperator88_2kbright_tally,
     recip_data_to_bin$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,10),ylim=c(0,10))
abline(0,1)







plot(recip_rank6diff_0$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_0$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_1$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_1$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_2$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_2$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_3$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_3$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_4$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_4$vector2_stoperator88_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_5$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_5$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


plot(recip_rank6diff_6$vector1_stoperator88_2kbright_tally,
     recip_rank6diff_6$vector2_stoperator88_2kbright_tally,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black',xlim=c(0,10),ylim=c(0,10))


#



#nontargetself stop sites

#







###Reciprocal analysis above








































###Compare lysogen vs repressor clone data








#Lysogen-Repressor analyses using averaged data
clone_match_columns <- c('defending_challenging',
                         'defending_phage',
                         'challenging_phage',
                         'averaged_infection_strength',
                         'averaged_turbidity',
                         'averaged_plaque_size',
                         'averaged_plaques',
                         'averaged_rank6',
                         'modified_mash_distance',
                         'pham_pham_dissimilarity',
                         'repressor_cterm_mafft_dist_uncorrected',
                         'stoperator_pwd_dist_euc',
                         'gene_content_clade_compare')

match_average_lysogen_and_clone_data <- function(lysogen_data,clone_data){
  
  
  lysY_reduced <- subset(lysogen_data,select = c(clone_match_columns))
  names(lysY_reduced) <- paste('lys_',names(lysY_reduced),sep="")

  cloneY_reduced <- subset(clone_data,select = c(clone_match_columns))
  names(cloneY_reduced) <- paste('clone_',names(cloneY_reduced),sep="")
  

  
  lys_clone_compare <- merge(lysY_reduced,cloneY_reduced,by.x='lys_defending_challenging',by.y='clone_defending_challenging')
  
  lys_clone_compare$lys_averaged_infection_strength <- as.numeric(as.character(lys_clone_compare$lys_averaged_infection_strength))
  lys_clone_compare$lys_averaged_turbidity <- as.numeric(as.character(lys_clone_compare$lys_averaged_turbidity))
  lys_clone_compare$lys_averaged_plaque_size <- as.numeric(as.character(lys_clone_compare$lys_averaged_plaque_size))
  lys_clone_compare$lys_averaged_plaques <- as.numeric(as.character(lys_clone_compare$lys_averaged_plaques))
  lys_clone_compare$lys_averaged_rank6 <- as.numeric(as.character(lys_clone_compare$lys_averaged_rank6))


  lys_clone_compare$clone_averaged_infection_strength <- as.numeric(as.character(lys_clone_compare$clone_averaged_infection_strength))
  lys_clone_compare$clone_averaged_turbidity <- as.numeric(as.character(lys_clone_compare$clone_averaged_turbidity))
  lys_clone_compare$clone_averaged_plaque_size <- as.numeric(as.character(lys_clone_compare$clone_averaged_plaque_size))
  lys_clone_compare$clone_averaged_plaques <- as.numeric(as.character(lys_clone_compare$clone_averaged_plaques))
  lys_clone_compare$clone_averaged_rank6 <- as.numeric(as.character(lys_clone_compare$clone_averaged_rank6))
  
  
  
  lys_clone_compare$infection_strength_diff <- lys_clone_compare$clone_averaged_infection_strength - lys_clone_compare$lys_averaged_infection_strength
  lys_clone_compare$turbidity_diff <- lys_clone_compare$clone_averaged_turbidity - lys_clone_compare$lys_averaged_turbidity
  lys_clone_compare$plaque_size_diff <- lys_clone_compare$clone_averaged_plaque_size - lys_clone_compare$lys_averaged_plaque_size
  lys_clone_compare$plaque_diff <- lys_clone_compare$clone_averaged_plaques - lys_clone_compare$lys_averaged_plaques
  lys_clone_compare$rank6_diff <- lys_clone_compare$clone_averaged_rank6 - lys_clone_compare$lys_averaged_rank6
  


  return(lys_clone_compare)
  
  
}






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

par(mar=c(4,8,4,4))
plot(ave_multi_conf_lys_clone_matched$lys_averaged_rank6,
     ave_multi_conf_lys_clone_matched$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6))
dev.copy (pdf,'.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_averaged_rank6,
     ave_multi_conf_env_rep_emp_lys_clone_matched$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6))
abline(0,1)
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6compare.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(ave_multi_conf_envN_lys_clone_matched$lys_averaged_rank6,
     ave_multi_conf_envN_lys_clone_matched$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6))
abline(0,1)
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6compare.pdf')
dev.off()


#Final correlation stats
lm_lys_clone_env_averaged_rank6 <- lm(clone_averaged_rank6 ~
                                        lys_averaged_rank6,
                                      data = ave_multi_conf_env_rep_emp_lys_clone_matched)
summary(lm_lys_clone_env_averaged_rank6)

lm_lys_clone_envN_averaged_rank6 <- lm(clone_averaged_rank6 ~
                                         lys_averaged_rank6,
                                       data = ave_multi_conf_envN_lys_clone_matched)
summary(lm_lys_clone_envN_averaged_rank6)
#






#

#FinalFig
par(mar=c(4,8,8,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic$lys_averaged_rank6,
     ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade$lys_averaged_rank6,
     ave_multi_conf_env_rep_emp_lys_clone_matched_interclade$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
par(new=TRUE)
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic$lys_averaged_rank6,
     ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6compare_subtypes.pdf')
dev.off()

#Final Number of assays used:
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic) +
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade) +
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic)



#FinalFig
par(mar=c(4,8,8,4))
plot(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic$lys_averaged_rank6,
     ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(ave_multi_conf_envN_lys_clone_matched_interclade$lys_averaged_rank6,
     ave_multi_conf_envN_lys_clone_matched_interclade$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
par(new=TRUE)
plot(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic$lys_averaged_rank6,
     ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6compare_subtypes.pdf')
dev.off()


#Final Number of assays used:
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic) +
  nrow(ave_multi_conf_envN_lys_clone_matched_interclade) +
  nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_envN_lys_clone_matched_interclade)
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)

#












par(mar=c(4,8,4,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_averaged_rank6,
     ave_multi_conf_env_rep_emp_lys_clone_matched$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
par(new=TRUE)
plot(ave_multi_conf_envN_lys_clone_matched$lys_averaged_rank6,
     ave_multi_conf_envN_lys_clone_matched$clone_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
dev.copy (pdf,'.pdf')
dev.off()







par(mar=c(4,8,4,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_modified_mash_distance,
     ave_multi_conf_env_rep_emp_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,0.4),ylim=c(-4,4))
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_mash.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(ave_multi_conf_envN_lys_clone_matched$lys_modified_mash_distance,
     ave_multi_conf_envN_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,0.4),ylim=c(-4,4))
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_mash.pdf')
dev.off()





par(mar=c(4,8,4,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_repressor_cterm_mafft_dist_uncorrected,
     ave_multi_conf_env_rep_emp_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(-4,4))
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_repCtermUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(ave_multi_conf_envN_lys_clone_matched$lys_repressor_cterm_mafft_dist_uncorrected,
     ave_multi_conf_envN_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,ylim=c(-4,4))
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_repCtermUn.pdf')
dev.off()







par(mar=c(4,8,4,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_pham_pham_dissimilarity,
     ave_multi_conf_env_rep_emp_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4))
abline(h=0)
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_GCD.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(ave_multi_conf_envN_lys_clone_matched$lys_pham_pham_dissimilarity,
     ave_multi_conf_envN_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4))
abline(h=0)
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_GCD.pdf')
dev.off()








par(mar=c(4,8,8,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_env_rep_emp_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4))
abline(h=0)
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(ave_multi_conf_envN_lys_clone_matched$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_envN_lys_clone_matched$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4))
abline(h=0)
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_stopEuc.pdf')
dev.off()









#FinalFig
par(mar=c(4,8,16,4))
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_env_rep_emp_lys_clone_matched_interclade$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="grey")
par(new=TRUE)
plot(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'ave_multi_conf_env_rep_emp_lys_clone_matched_rank6diff_vs_stopEuc_subtypes.pdf')
dev.off()

#Final Number of assays used:
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic) +
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade) +
  nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_interclade)
nrow(ave_multi_conf_env_rep_emp_lys_clone_matched_intraclade2_homotypic)




#FinalFig
par(mar=c(4,8,16,4))
plot(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(ave_multi_conf_envN_lys_clone_matched_interclade$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_envN_lys_clone_matched_interclade$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="grey")
par(new=TRUE)
plot(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic$lys_stoperator_pwd_dist_euc,
     ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic$rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'ave_multi_conf_envN_lys_clone_matched_rank6diff_vs_stopEuc_subtypes.pdf')
dev.off()


#Final Number of assays used:
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic) +
  nrow(ave_multi_conf_envN_lys_clone_matched_interclade) +
  nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_heterotypic)
nrow(ave_multi_conf_envN_lys_clone_matched_interclade)
nrow(ave_multi_conf_envN_lys_clone_matched_intraclade2_homotypic)


###Lysogen-clone comparison above





















###Compare L5 and phiTM41 defense profiles


#Analysis using averaged data
#This can't be mixed with code for non-averaged data because averaged dataframe columns are named differently
#Averaged data should only reflect confident data
lys_l5_ave <- subset(conf_assay_strain_def_chal_average,
                     conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                       conf_assay_strain_def_chal_average$defending_phage == 'l5' &
                       conf_assay_strain_def_chal_average$strain_type == 'lysogen')

lys_phitm41_ave <- subset(conf_assay_strain_def_chal_average,
                          conf_assay_strain_def_chal_average$assay_type == 'multiple_titer' &
                            conf_assay_strain_def_chal_average$defending_phage == 'phitm41' &
                            conf_assay_strain_def_chal_average$strain_type == 'lysogen')





#Choose which dataset to reduce
l5_ave_reduced <- subset(lys_l5_ave,
                         select = c(
                           'defending_challenging',
                           'challenging_phage',
                           'averaged_infection_strength',
                           'averaged_turbidity',
                           'averaged_plaque_size',
                           'averaged_plaques',
                           'averaged_four_factors',
                           'averaged_rank6'))

names(l5_ave_reduced) <- c('l5_defending_challenging',
                       'challenging_phage',
                       'l5_averaged_infection_strength',
                       'l5_averaged_turbidity',
                       'l5_averaged_plaque_size',
                       'l5_averaged_plaques',
                       'l5_averaged_four_factors',
                       'l5_averaged_rank6')



phitm41_ave_reduced <- subset(lys_phitm41_ave,
                          select = c(
                            'defending_challenging',
                            'challenging_phage',
                            'averaged_infection_strength',
                            'averaged_turbidity',
                            'averaged_plaque_size',
                            'averaged_plaques',
                            'averaged_four_factors',
                            'averaged_rank6'))

names(phitm41_ave_reduced) <- c('phitm41_defending_challenging',
                            'challenging_phage',
                            'phitm41_averaged_infection_strength',
                            'phitm41_averaged_turbidity',
                            'phitm41_averaged_plaque_size',
                            'phitm41_averaged_plaques',
                            'phitm41_averaged_four_factors',
                            'phitm41_averaged_rank6')



l5_phitm41_ave_compare <- merge(l5_ave_reduced,phitm41_ave_reduced,by.x='challenging_phage',by.y='challenging_phage')

l5_phitm41_ave_compare$l5_averaged_infection_strength <- as.numeric(as.character(l5_phitm41_ave_compare$l5_averaged_infection_strength))
l5_phitm41_ave_compare$phitm41_averaged_infection_strength <- as.numeric(as.character(l5_phitm41_ave_compare$phitm41_averaged_infection_strength))

l5_phitm41_ave_compare$l5_averaged_rank6 <- as.numeric(as.character(l5_phitm41_ave_compare$l5_averaged_rank6))
l5_phitm41_ave_compare$phitm41_averaged_rank6 <- as.numeric(as.character(l5_phitm41_ave_compare$phitm41_averaged_rank6))

l5_phitm41_ave_compare$averaged_rank6_diff <- l5_phitm41_ave_compare$l5_averaged_rank6 - l5_phitm41_ave_compare$phitm41_averaged_rank6

setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,4,4))
hist(l5_phitm41_ave_compare$averaged_rank6_diff,
     xlim=c(-1.5,1.5),col='black',ann=FALSE,main=NULL,las=1,breaks=15)
dev.copy(pdf,'l5_phitm41_ave_compare_averaged_rank6_diff.pdf')
dev.off()


###Compare L5 and phiTM41 defense profiles above


















###L5-mutant comparisons using all averaged data (look at all L5 mutants, not just phiTM41)



#All averaged multiple-titer lysogen and clone data
l5_columns <- c('defending_challenging',
                'defending_phage',
                'challenging_phage',
                'averaged_infection_strength',
                'averaged_turbidity',
                'averaged_plaque_size',
                'averaged_plaques',
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





#Plots
setwd("~/scratch/immunity_analysis/output/")


#
par(mar=c(4,8,8,4))
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
par(new=TRUE)
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm4_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="green")
par(new=TRUE)
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="cyan")
abline(0,1)
dev.copy (pdf,'chal_l5_assays_ave_rank6_l5_vs_mutants.pdf')
dev.off()

#Number of assays used:
nrow(chal_l5_assays)







#
par(mar=c(4,8,8,4))
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm41_rank6.pdf')
dev.off()

#number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm41_averaged_rank6)))
                            



#
par(mar=c(4,8,8,4))
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm1_rank6.pdf')
dev.off()

#number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm1_averaged_rank6)))





# par(mar=c(4,8,8,4))
# plot(chal_l5_assays$l5_averaged_rank6,
#      chal_l5_assays$phitm4_averaged_rank6,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
# abline(0,1)
# dev.copy (pdf,'chal_l5_phitm4_rank6.pdf')
# dev.off()






#
par(mar=c(4,8,8,4))
plot(chal_l5_assays$l5_averaged_rank6,
     chal_l5_assays$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm6_rank6.pdf')
dev.off()

#number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm6_averaged_rank6)))





# 
# par(mar=c(4,8,8,4))
# plot(chal_l5_assays$phitm1_averaged_rank6,
#      chal_l5_assays$phitm4_averaged_rank6,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
# abline(0,1)
# dev.copy (pdf,'chal_phitm1_phitm4_rank6.pdf')
# dev.off()
# 
# 
# 
# par(mar=c(4,8,8,4))
# plot(chal_l5_assays$phitm1_averaged_rank6,
#      chal_l5_assays$phitm6_averaged_rank6,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
# abline(0,1)
# dev.copy (pdf,'chal_phitm1_phitm6_rank6.pdf')
# dev.off()
# 















par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_pham_pham_dissimilarity,
     chal_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="black")
par(new=TRUE)
plot(chal_l5_assays$l5_pham_pham_dissimilarity,
     chal_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(chal_l5_assays$l5_pham_pham_dissimilarity,
     chal_l5_assays$l5_phitm4_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="green")
par(new=TRUE)
plot(chal_l5_assays$l5_pham_pham_dissimilarity,
     chal_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="cyan")
abline(h=0)
dev.copy (pdf,'chal_l5_assays_ave_rank6_vs_gcd_l5_vs_mutants.pdf')
dev.off()




#
par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
par(new=TRUE)
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm4_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="green")
par(new=TRUE)
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="cyan")
abline(h=0)
dev.copy (pdf,'chal_l5_assays_ave_rank6_vs_stopEuc_l5_vs_mutants.pdf')
dev.off()










par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'chal_l5_phitm41_rankDiff_vs_stopEuc.pdf')
dev.off()


par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'chal_l5_phitm1_rankDiff_vs_stopEuc.pdf')
dev.off()


par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm4_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'chal_l5_phitm4_rankDiff_vs_stopEuc.pdf')
dev.off()


par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'chal_l5_phitm6_rankDiff_vs_stopEuc.pdf')
dev.off()






par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$phitm1_phitm4_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'chal_phitm1_phitm4_rankDiff_vs_stopEuc.pdf')
dev.off()












#
par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="black")
dev.copy (pdf,'chal_phitm1_rank6_vs_stopEuc.pdf')
dev.off()

#number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm1_averaged_rank6)))


#
par(mar=c(4,8,16,4))
plot(chal_l5_assays$l5_stoperator_pwd_dist_euc,
     chal_l5_assays$phitm4_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="black")
dev.copy (pdf,'chal_phitm4_rank6_vs_stopEuc.pdf')
dev.off()


#number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm4_averaged_rank6)))











#Same analyses but splitting by clade2 phages
chal_l5_assays_nonclade2 <- subset(chal_l5_assays,chal_l5_assays$l5_defending_gene_content_clade != "clade2")
chal_l5_assays_clade2 <- subset(chal_l5_assays,chal_l5_assays$l5_defending_gene_content_clade == "clade2")





#Right now I only need to compute homotypic comparisons for phiTM1
chal_l5_assays_clade2_phitm1homotypic <- subset(chal_l5_assays_clade2,
                                                as.character(chal_l5_assays_clade2$phitm1_challenging_phage) == 
                                                  as.character(chal_l5_assays_clade2$l5_defending_phage))
chal_l5_assays_clade2_phitm1heterotypic <- subset(chal_l5_assays_clade2,
                                                as.character(chal_l5_assays_clade2$phitm1_challenging_phage) != 
                                                  as.character(chal_l5_assays_clade2$l5_defending_phage))

#QC
nrow(chal_l5_assays_clade2)
nrow(subset(chal_l5_assays_clade2,is.na(chal_l5_assays_clade2$phitm1_challenging_phage)))
nrow(chal_l5_assays_clade2_phitm1homotypic)
nrow(chal_l5_assays_clade2_phitm1heterotypic)
#








#FinalFig
par(mar=c(4,8,8,4))
plot(chal_l5_assays_clade2$l5_averaged_rank6,
     chal_l5_assays_clade2$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(chal_l5_assays_nonclade2$l5_averaged_rank6,
     chal_l5_assays_nonclade2$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm41_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm41_averaged_rank6)))
nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm41_averaged_rank6)))
nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm41_averaged_rank6)))



#FinalFig
par(mar=c(4,8,8,4))
plot(chal_l5_assays_clade2$l5_averaged_rank6,
     chal_l5_assays_clade2$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(chal_l5_assays_nonclade2$l5_averaged_rank6,
     chal_l5_assays_nonclade2$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm1_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm1_averaged_rank6)))
nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm1_averaged_rank6)))
nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm1_averaged_rank6)))






#FinalFig
par(mar=c(4,8,8,4))
plot(chal_l5_assays_clade2$l5_averaged_rank6,
     chal_l5_assays_clade2$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(chal_l5_assays_nonclade2$l5_averaged_rank6,
     chal_l5_assays_nonclade2$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'chal_l5_phitm6_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm6_averaged_rank6)))
nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm6_averaged_rank6)))
nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm6_averaged_rank6)))







# #
# par(mar=c(4,8,16,4))
# plot(chal_l5_assays_clade2$l5_stoperator_pwd_dist_euc,
#      chal_l5_assays_clade2$phitm1_averaged_rank6,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="red")
# par(new=TRUE)
# plot(chal_l5_assays_nonclade2$l5_stoperator_pwd_dist_euc,
#      chal_l5_assays_nonclade2$phitm1_averaged_rank6,
#      pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="grey")
# dev.copy (pdf,'chal_phitm1_rank6_vs_stopEuc_by_subtype.pdf')
# dev.off()
# 
# #number of assays
# nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm1_averaged_rank6)))
# nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm1_averaged_rank6)))
# nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm1_averaged_rank6)))







#FinalFig
par(mar=c(4,8,16,4))
plot(chal_l5_assays_clade2_phitm1heterotypic$l5_stoperator_pwd_dist_euc,
     chal_l5_assays_clade2_phitm1heterotypic$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="red")
par(new=TRUE)
par(mar=c(4,8,16,4))
plot(chal_l5_assays_clade2_phitm1homotypic$l5_stoperator_pwd_dist_euc,
     chal_l5_assays_clade2_phitm1homotypic$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="black")
par(new=TRUE)
plot(chal_l5_assays_nonclade2$l5_stoperator_pwd_dist_euc,
     chal_l5_assays_nonclade2$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="grey")
dev.copy (pdf,'chal_phitm1_rank6_vs_stopEuc_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm1_averaged_rank6)))
nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm1_averaged_rank6)))
nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm1_averaged_rank6)))
nrow(chal_l5_assays_clade2_phitm1homotypic)
nrow(chal_l5_assays_clade2_phitm1heterotypic)




#FinalFig
par(mar=c(4,8,16,4))
plot(chal_l5_assays_clade2$l5_stoperator_pwd_dist_euc,
     chal_l5_assays_clade2$phitm4_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="red")
par(new=TRUE)
plot(chal_l5_assays_nonclade2$l5_stoperator_pwd_dist_euc,
     chal_l5_assays_nonclade2$phitm4_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(0,6),col="grey")
dev.copy (pdf,'chal_phitm4_rank6_vs_stopEuc_by_subtype.pdf')
dev.off()


#Final number of assays
nrow(subset(chal_l5_assays,!is.na(chal_l5_assays$phitm4_averaged_rank6)))
nrow(subset(chal_l5_assays_clade2,!is.na(chal_l5_assays_clade2$phitm4_averaged_rank6)))
nrow(subset(chal_l5_assays_nonclade2,!is.na(chal_l5_assays_nonclade2$phitm4_averaged_rank6)))


























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


#Plots



#
par(mar=c(4,8,8,4))
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
par(new=TRUE)
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="cyan")
abline(0,1)
dev.copy (pdf,'def_l5_assays_ave_rank6_l5_vs_mutants.pdf')
dev.off()

#Number of assays used:
nrow(def_l5_assays)










#
par(mar=c(4,8,8,4))
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'def_l5_phitm41_rank6.pdf')
dev.off()

#number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm41_averaged_rank6)))




#
par(mar=c(4,8,8,4))
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'def_l5_phitm1_rank6.pdf')
dev.off()

#number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm1_averaged_rank6)))




#
par(mar=c(4,8,8,4))
plot(def_l5_assays$l5_averaged_rank6,
     def_l5_assays$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'def_l5_phitm6_rank6.pdf')
dev.off()

#number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm6_averaged_rank6)))









par(mar=c(4,8,8,4))
plot(def_l5_assays$phitm1_averaged_rank6,
     def_l5_assays$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="black")
abline(0,1)
dev.copy (pdf,'def_phitm1_phitm6_rank6.pdf')
dev.off()













par(mar=c(4,8,8,4))
plot(def_l5_assays$l5_pham_pham_dissimilarity,
     def_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="black")
par(new=TRUE)
plot(def_l5_assays$l5_pham_pham_dissimilarity,
     def_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(def_l5_assays$l5_pham_pham_dissimilarity,
     def_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,1),ylim=c(-4,4),col="cyan")
abline(h=0)
dev.copy (pdf,'def_l5_assays_ave_rank6_vs_gcd_l5_vs_mutants.pdf')
dev.off()


#
par(mar=c(4,8,16,4))
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
par(new=TRUE)
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="red")
par(new=TRUE)
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="cyan")
abline(h=0)
dev.copy (pdf,'def_l5_assays_ave_rank6_vs_stopEuc_l5_vs_mutants.pdf')
dev.off()










par(mar=c(4,8,16,4))
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm41_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'def_l5_phiTM41_rankDiff_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,16,4))
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm1_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'def_l5_phiTM1_rankDiff_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,16,4))
plot(def_l5_assays$l5_stoperator_pwd_dist_euc,
     def_l5_assays$l5_phitm6_rank6_diff,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,5),ylim=c(-4,4),col="black")
abline(h=0)
dev.copy (pdf,'def_l5_phiTM6_rankDiff_vs_stopEuc.pdf')
dev.off()












#Same analyses but splitting by clade2 phages
def_l5_assays_clade2 <- subset(def_l5_assays,def_l5_assays$l5_challenging_gene_content_clade == "clade2")
def_l5_assays_nonclade2 <- subset(def_l5_assays,def_l5_assays$l5_challenging_gene_content_clade != "clade2")

nrow(def_l5_assays)
nrow(def_l5_assays_clade2)
nrow(def_l5_assays_nonclade2)
summary(def_l5_assays$l5_challenging_phage)
summary(def_l5_assays_clade2$l5_challenging_phage)
summary(def_l5_assays_nonclade2$l5_challenging_phage)





#FinalFig
par(mar=c(4,8,8,4))
plot(def_l5_assays_clade2$l5_averaged_rank6,
     def_l5_assays_clade2$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(def_l5_assays_nonclade2$l5_averaged_rank6,
     def_l5_assays_nonclade2$phitm41_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'def_l5_phitm41_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm41_averaged_rank6)))
nrow(subset(def_l5_assays_clade2,!is.na(def_l5_assays_clade2$phitm41_averaged_rank6)))
nrow(subset(def_l5_assays_nonclade2,!is.na(def_l5_assays_nonclade2$phitm41_averaged_rank6)))


#FinalFig
par(mar=c(4,8,8,4))
plot(def_l5_assays_clade2$l5_averaged_rank6,
     def_l5_assays_clade2$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(def_l5_assays_nonclade2$l5_averaged_rank6,
     def_l5_assays_nonclade2$phitm1_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'def_l5_phitm1_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm1_averaged_rank6)))
nrow(subset(def_l5_assays_clade2,!is.na(def_l5_assays_clade2$phitm1_averaged_rank6)))
nrow(subset(def_l5_assays_nonclade2,!is.na(def_l5_assays_nonclade2$phitm1_averaged_rank6)))


#FinalFig
par(mar=c(4,8,8,4))
plot(def_l5_assays_clade2$l5_averaged_rank6,
     def_l5_assays_clade2$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="red")
par(new=TRUE)
plot(def_l5_assays_nonclade2$l5_averaged_rank6,
     def_l5_assays_nonclade2$phitm6_averaged_rank6,
     pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,xlim=c(0,6),ylim=c(0,6),col="grey")
abline(0,1)
dev.copy (pdf,'def_l5_phitm6_rank6_by_subtype.pdf')
dev.off()

#Final number of assays
nrow(subset(def_l5_assays,!is.na(def_l5_assays$phitm6_averaged_rank6)))
nrow(subset(def_l5_assays_clade2,!is.na(def_l5_assays_clade2$phitm6_averaged_rank6)))
nrow(subset(def_l5_assays_nonclade2,!is.na(def_l5_assays_nonclade2$phitm6_averaged_rank6)))










###L5-mutant comparisons using all averaged data (look at all L5 mutants, not just phiTM41) above












































###Compare known empirical temperate to isolated mutant to escape mutant

#Mutant-parent phage sets
#1
#pioneer
#phitm35 (escape mutant)

#2
#eagleeye
#phitm36 (escape mutant)

#3
#starstuff
#d29
#phitm43 (isolated mutant)
#phitm44 (isolated mutant)
#phitm38 (escape mutant)

#4
#et2brutus
#phitm39 (escape mutant)	
#phitm40 (escape mutant)


#5
#l5
#phitm1 (engineered mutant)
#phitm4 (isolated mutant)
#phitm6 (engineered mutant)
#phitm41 (escape mutant)


#6
#trixie
#phitm42 (escape mutant)

#7 (clone EM)
#bxb1
#phitm45 (escape mutant)

#8 (clone EM)
#davinci
#phitm46 (escape mutant)

#9 (clone EM)
#gladiator
#phitm47 (escape mutant)

###


















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
mutant_analysis$averaged_infection_strength_diff <- mutant_analysis$mutant_averaged_infection_strength - mutant_analysis$parent_averaged_infection_strength
mutant_analysis$averaged_four_factors_diff <- mutant_analysis$mutant_averaged_four_factors - mutant_analysis$parent_averaged_four_factors
mutant_analysis$averaged_rank6_diff <- mutant_analysis$mutant_averaged_rank6 - mutant_analysis$parent_averaged_rank6











setwd("~/scratch/immunity_analysis/output/")
write.table(mutant_analysis,
            "mutant_analysis.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)






par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     mutant_analysis$parent_averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'parent_ave_rank6_vs_cas4MafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_ave_rank6_vs_cas4MafftUn.pdf')
dev.off()









par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$parent_averaged_infection_strength)),
     xlim=c(0,0.4),ylim=c(0,4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue')
par(new=TRUE)
plot(mutant_analysis$mutant_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$mutant_averaged_infection_strength)),
     xlim=c(0,0.4),ylim=c(0,4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$parent_averaged_four_factors)),
     xlim=c(0,0.4),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue')
par(new=TRUE)
plot(mutant_analysis$mutant_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$mutant_averaged_four_factors)),
     xlim=c(0,0.4),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'mutant_parent_.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_repressor_muscle_bionj_distances,
     as.numeric(as.character(mutant_analysis$parent_averaged_four_factors)),
     xlim=c(0,1),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue')
par(new=TRUE)
plot(mutant_analysis$parent_repressor_muscle_bionj_distances,
     as.numeric(as.character(mutant_analysis$mutant_averaged_four_factors)),
     xlim=c(0,1),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'.pdf')
dev.off()






par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_repressor_full_mafft_dist_uncorrected,
     as.numeric(as.character(mutant_analysis$parent_averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='blue')
par(new=TRUE)
plot(mutant_analysis$parent_repressor_full_mafft_dist_uncorrected,
     as.numeric(as.character(mutant_analysis$mutant_averaged_rank6)),
     xlim=c(0,60),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'.pdf')
dev.off()







par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$averaged_infection_strength_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_parent_ave_inf_strength_diff_by_mash.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(mutant_analysis$parent_modified_mash_distance,
     as.numeric(as.character(mutant_analysis$averaged_four_factors_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_parent_ave_four_factors_diff_by_mash.pdf')
dev.off()




hist(mutant_analysis$averaged_rank6_diff,col='black',ann=FALSE,main=NULL,las=1,breaks=10)
dev.copy(pdf,'mutant_parent_ave_rank6_diff_hist.pdf')
dev.off()












#Bin data

mutant_data_to_bin <- mutant_analysis



#Choose a sub-dataset:
mutant_rank6diff_0 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff < 0.5)

mutant_rank6diff_1 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 0.5 &
                           mutant_data_to_bin$averaged_rank6_diff < 1.5)

mutant_rank6diff_2 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 1.5 &
                           mutant_data_to_bin$averaged_rank6_diff < 2.5)

mutant_rank6diff_3 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 2.5 &
                           mutant_data_to_bin$averaged_rank6_diff < 3.5)

mutant_rank6diff_4 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 3.5 &
                           mutant_data_to_bin$averaged_rank6_diff < 4.5)

mutant_rank6diff_5 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 4.5 &
                           mutant_data_to_bin$averaged_rank6_diff < 5.5)

mutant_rank6diff_6 <- subset(mutant_data_to_bin,
                         mutant_data_to_bin$averaged_rank6_diff >= 5.5)

#QC Check
nrow(mutant_data_to_bin) - 
  (nrow(mutant_rank6diff_0) + 
     nrow(mutant_rank6diff_1) + 
     nrow(mutant_rank6diff_2) + 
     nrow(mutant_rank6diff_3) + 
     nrow(mutant_rank6diff_4) + 
     nrow(mutant_rank6diff_5) + 
     nrow(mutant_rank6diff_6))


mutant_rank6diff_0_freq <- data.frame("bin0",nrow(mutant_rank6diff_0))
mutant_rank6diff_1_freq <- data.frame("bin1",nrow(mutant_rank6diff_1))
mutant_rank6diff_2_freq <- data.frame("bin2",nrow(mutant_rank6diff_2))
mutant_rank6diff_3_freq <- data.frame("bin3",nrow(mutant_rank6diff_3))
mutant_rank6diff_4_freq <- data.frame("bin4",nrow(mutant_rank6diff_4))
mutant_rank6diff_5_freq <- data.frame("bin5",nrow(mutant_rank6diff_5))
mutant_rank6diff_6_freq <- data.frame("bin6",nrow(mutant_rank6diff_6))


names(mutant_rank6diff_0_freq) <- c("bin","freq")
names(mutant_rank6diff_1_freq) <- c("bin","freq")
names(mutant_rank6diff_2_freq) <- c("bin","freq")
names(mutant_rank6diff_3_freq) <- c("bin","freq")
names(mutant_rank6diff_4_freq) <- c("bin","freq")
names(mutant_rank6diff_5_freq) <- c("bin","freq")
names(mutant_rank6diff_6_freq) <- c("bin","freq")


mutant_rank6diff_0_freq$freq_percent <- mutant_rank6diff_0_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_1_freq$freq_percent <- mutant_rank6diff_1_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_2_freq$freq_percent <- mutant_rank6diff_2_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_3_freq$freq_percent <- mutant_rank6diff_3_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_4_freq$freq_percent <- mutant_rank6diff_4_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_5_freq$freq_percent <- mutant_rank6diff_5_freq$freq/nrow(mutant_data_to_bin)
mutant_rank6diff_6_freq$freq_percent <- mutant_rank6diff_6_freq$freq/nrow(mutant_data_to_bin)


mutant_rank6diff_binned_freq <- rbind(mutant_rank6diff_0_freq,
                                      mutant_rank6diff_1_freq,
                                      mutant_rank6diff_2_freq,
                                      mutant_rank6diff_3_freq,
                                      mutant_rank6diff_4_freq,
                                      mutant_rank6diff_5_freq,
                                      mutant_rank6diff_6_freq)


mutant_rank6diff_binned_freq$bin <- factor(mutant_rank6diff_binned_freq$bin,
                                           c("bin6",
                                             "bin5",
                                             "bin4",
                                             "bin3",
                                             "bin2",
                                             "bin1",
                                             "bin0"))


#QC
sum(mutant_rank6diff_binned_freq$freq) - nrow(mutant_data_to_bin)






par(mar=c(4,8,4,4))
barplot(mutant_rank6diff_binned_freq$freq,
        names.arg=mutant_rank6diff_binned_freq$bin,
        col='black')
dev.copy(pdf,'mutant_rank6diff_binned_freq_total.pdf')
dev.off()


par(mar=c(4,8,4,4))
barplot(mutant_rank6diff_binned_freq$freq_percent,
        names.arg=mutant_rank6diff_binned_freq$bin,
        col='black',ylim=c(0,0.6))
dev.copy(pdf,'mutant_rank6diff_binned_freq_percent.pdf')
dev.off()












#Conclusions: in general, mutants infect with equal or greater strength than the parent.
#Some good examples:
#1. phiTM4 and phiTM6 infecting L5
#2. phiTM42 infecting Trixie and Et2Brutus
#3. phiTM38 infecting StarStuff
#4. phiTM41 infecting Trixie and DarthPhader (moderate but not strong)

















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
nrow(escape_mutant_analysis_intraclade2_parenthomotypic)
nrow(escape_mutant_analysis_intraclade2_parentheterotypic)
nrow(escape_mutant_analysis_intraclade2_mutanthomotypic)
nrow(escape_mutant_analysis_intraclade2_mutantheterotypic)

summary(escape_mutant_analysis_intraclade2_parenthomotypic$parent_defending_phage)
summary(escape_mutant_analysis_intraclade2_mutanthomotypic$mutant_defending_phage)
#




par(mar=c(4,8,4,4))
plot(escape_mutant_analysis$parent_pham_pham_dissimilarity,
     escape_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm35_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm35_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm36_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm36_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm38_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm38_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm39_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm40_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm40_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm41_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm41_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm42_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm42_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm46_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm46_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm47_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm47_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')









par(mar=c(4,8,4,4))
plot(phitm35_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm35_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm36_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm36_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')


par(mar=c(4,8,4,4))
plot(phitm38_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm38_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm39_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm40_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm40_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm41_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm41_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm42_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm42_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')

par(mar=c(4,8,4,4))
plot(phitm46_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm46_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')


par(mar=c(4,8,4,4))
plot(phitm47_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm47_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')









par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_pham_pham_dissimilarity,
     phitm39_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,1),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_cas4_mafft_dist_uncorrected,
     phitm39_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_modified_mash_distance,
     phitm39_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,0.5),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_repressor_cterm_mafft_dist_uncorrected,
     phitm39_mutant_analysis$averaged_rank6_diff,
     ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')



par(mar=c(4,8,8,4))
plot(escape_mutant_analysis$parent_averaged_rank6,
     escape_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
dev.copy(pdf,'ave_rank6_parent_vs_mutant.pdf')
dev.off()





#FinalFig
par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_intraclade2$parent_averaged_rank6,
     escape_mutant_analysis_intraclade2$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
par(new=TRUE)
plot(escape_mutant_analysis_interclade$parent_averaged_rank6,
     escape_mutant_analysis_interclade$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
abline(0,1)
dev.copy(pdf,'ave_rank6_parent_vs_mutant_subtypes.pdf')
dev.off()


# #
# par(mar=c(4,8,16,4))
# plot(escape_mutant_analysis_intraclade2$parent_stoperator_pwd_dist_euc,
#      escape_mutant_analysis_intraclade2$parent_averaged_rank6,
#      xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
# par(new=TRUE)
# plot(escape_mutant_analysis_interclade$parent_stoperator_pwd_dist_euc,
#      escape_mutant_analysis_interclade$parent_averaged_rank6,
#      xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
# dev.copy(pdf,'ave_rank6_parent_vs_parent_stopEuc.pdf')
# dev.off()
# 
# 
# #
# par(mar=c(4,8,16,4))
# plot(escape_mutant_analysis_intraclade2$mutant_stoperator_pwd_dist_euc,
#      escape_mutant_analysis_intraclade2$mutant_averaged_rank6,
#      xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
# par(new=TRUE)
# plot(escape_mutant_analysis_interclade$mutant_stoperator_pwd_dist_euc,
#      escape_mutant_analysis_interclade$mutant_averaged_rank6,
#      xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
# dev.copy(pdf,'ave_rank6_mutant_vs_mutant_stopEuc.pdf')
# dev.off()
# 
# 
# 
# 
# #Number of assays used:
# nrow(escape_mutant_analysis_intraclade2) +
#   nrow(escape_mutant_analysis_interclade)
# nrow(escape_mutant_analysis_intraclade2)
# nrow(escape_mutant_analysis_interclade)
# 
# 
# 
# 
# #correlations
# lm_stop_vs_averank6_parent <- lm(parent_averaged_rank6 ~
#                                    parent_stoperator_pwd_dist_euc,
#                                  data = escape_mutant_analysis_intraclade2)
# summary(lm_stop_vs_averank6_parent)
# 
# lm_stop_vs_averank6_mutant <- lm(mutant_averaged_rank6 ~
#                                    mutant_stoperator_pwd_dist_euc,
#                                  data = escape_mutant_analysis_intraclade2)
# summary(lm_stop_vs_averank6_mutant)









#FinalFig
par(mar=c(4,8,16,4))
plot(escape_mutant_analysis_intraclade2_parentheterotypic$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_intraclade2_parentheterotypic$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
par(new=TRUE)
plot(escape_mutant_analysis_intraclade2_parenthomotypic$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_intraclade2_parenthomotypic$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(escape_mutant_analysis_interclade$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_interclade$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'ave_rank6_parent_vs_parent_stopEuc.pdf')
dev.off()


#FinalFig
par(mar=c(4,8,16,4))
plot(escape_mutant_analysis_intraclade2_mutantheterotypic$mutant_stoperator_pwd_dist_euc,
     escape_mutant_analysis_intraclade2_mutantheterotypic$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
par(new=TRUE)
plot(escape_mutant_analysis_intraclade2_mutanthomotypic$mutant_stoperator_pwd_dist_euc,
     escape_mutant_analysis_intraclade2_mutanthomotypic$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(escape_mutant_analysis_interclade$mutant_stoperator_pwd_dist_euc,
     escape_mutant_analysis_interclade$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'ave_rank6_mutant_vs_mutant_stopEuc.pdf')
dev.off()




#Final Number of assays used:
nrow(escape_mutant_analysis_intraclade2) +
  nrow(escape_mutant_analysis_interclade)
nrow(escape_mutant_analysis_intraclade2)
nrow(escape_mutant_analysis_interclade)
nrow(escape_mutant_analysis_intraclade2_parentheterotypic)
nrow(escape_mutant_analysis_intraclade2_parenthomotypic)
nrow(escape_mutant_analysis_intraclade2_mutantheterotypic)
nrow(escape_mutant_analysis_intraclade2_mutanthomotypic)



#Final correlations
lm_stop_vs_averank6_parent <- lm(parent_averaged_rank6 ~
                                   parent_stoperator_pwd_dist_euc,
                                 data = escape_mutant_analysis_intraclade2_parentheterotypic)
summary(lm_stop_vs_averank6_parent)

lm_stop_vs_averank6_mutant <- lm(mutant_averaged_rank6 ~
                                   mutant_stoperator_pwd_dist_euc,
                                 data = escape_mutant_analysis_intraclade2_mutantheterotypic)
summary(lm_stop_vs_averank6_mutant)













par(mar=c(4,8,4,4))
plot(phitm35_mutant_analysis$parent_averaged_rank6,
     phitm35_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm36_mutant_analysis$parent_averaged_rank6,
     phitm36_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm38_mutant_analysis$parent_averaged_rank6,
     phitm38_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm39_mutant_analysis$parent_averaged_rank6,
     phitm39_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm40_mutant_analysis$parent_averaged_rank6,
     phitm40_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm41_mutant_analysis$parent_averaged_rank6,
     phitm41_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm42_mutant_analysis$parent_averaged_rank6,
     phitm42_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm46_mutant_analysis$parent_averaged_rank6,
     phitm46_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
par(mar=c(4,8,4,4))
plot(phitm47_mutant_analysis$parent_averaged_rank6,
     phitm47_mutant_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)






par(mar=c(4,8,4,4))
plot(escape_mutant_analysis$parent_averaged_rank6,
     escape_mutant_analysis$averaged_rank6_diff,
     xlim=c(0,6),ylim=c(-6,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')










escape_mutant_analysis_repN <- subset(escape_mutant_analysis,escape_mutant_analysis$mutant_challenging_phage != "phitm41")
escape_mutant_analysis_repY <- subset(escape_mutant_analysis,escape_mutant_analysis$mutant_challenging_phage == "phitm41")


escape_mutant_analysis_repN_intraclade2 <- subset(escape_mutant_analysis_repN,escape_mutant_analysis_repN$parent_gene_content_clade_compare == "clade2")
escape_mutant_analysis_repN_interclade <- subset(escape_mutant_analysis_repN,escape_mutant_analysis_repN$parent_gene_content_clade_compare == "different")


escape_mutant_analysis_repY_intraclade2 <- subset(escape_mutant_analysis_repY,escape_mutant_analysis_repY$parent_gene_content_clade_compare == "clade2")
escape_mutant_analysis_repY_interclade <- subset(escape_mutant_analysis_repY,escape_mutant_analysis_repY$parent_gene_content_clade_compare == "different")




par(mar=c(4,8,4,4))
plot(escape_mutant_analysis_repN$parent_cas4_mafft_dist_uncorrected,
     escape_mutant_analysis_repN$mutant_averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_repN_ave_rank6_vs_cas4MafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(escape_mutant_analysis_repY$parent_cas4_mafft_dist_uncorrected,
     escape_mutant_analysis_repY$mutant_averaged_rank6,
     xlim=c(0,50),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_repY_ave_rank6_vs_cas4MafftUn.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(escape_mutant_analysis_repN$parent_cas4_mafft_dist_uncorrected,
     escape_mutant_analysis_repN$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(h=0)
dev.copy(pdf,'mutant_repN_ave_rank6Diff_vs_cas4MafftUn.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(escape_mutant_analysis_repY$parent_cas4_mafft_dist_uncorrected,
     escape_mutant_analysis_repY$averaged_rank6_diff,
     xlim=c(0,50),ylim=c(-2,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(h=0)
dev.copy(pdf,'mutant_repY_ave_rank6Diff_vs_cas4MafftUn.pdf')
dev.off()











#
par(mar=c(4,8,8,4))
plot(escape_mutant_analysis$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'parent_ave_rank6_vs_stopEuc.pdf')
dev.off()




par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repN$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repN$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_repN_ave_rank6_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repY$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repY$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'mutant_repY_ave_rank6_vs_stopEuc.pdf')
dev.off()




par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repN$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repN$averaged_rank6_diff,
     xlim=c(0,5),ylim=c(-2,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(h=0)
dev.copy(pdf,'mutant_repN_ave_rank6Diff_vs_stopEuc.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repY$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repY$averaged_rank6_diff,
     xlim=c(0,5),ylim=c(-2,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(h=0)
dev.copy(pdf,'mutant_repY_ave_rank6Diff_vs_stopEuc.pdf')
dev.off()













par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_intraclade2$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_intraclade2$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
par(new=TRUE)
plot(escape_mutant_analysis_interclade$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_interclade$parent_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'parent_ave_rank6_vs_stopEuc_subtypes.pdf')
dev.off()



par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repN_intraclade2$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repN_intraclade2$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
par(new=TRUE)
plot(escape_mutant_analysis_repN_interclade$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repN_interclade$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'mutant_repN_ave_rank6_vs_stopEuc_subtypes.pdf')
dev.off()

par(mar=c(4,8,8,4))
plot(escape_mutant_analysis_repY_intraclade2$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repY_intraclade2$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
par(new=TRUE)
plot(escape_mutant_analysis_repY_interclade$parent_stoperator_pwd_dist_euc,
     escape_mutant_analysis_repY_interclade$mutant_averaged_rank6,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
dev.copy(pdf,'mutant_repY_ave_rank6_vs_stopEuc_subtypes.pdf')
dev.off()


#











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
mutant_rep_analysis$averaged_infection_strength_diff <- mutant_rep_analysis$mutant_averaged_infection_strength - mutant_rep_analysis$parent_averaged_infection_strength
mutant_rep_analysis$averaged_four_factors_diff <- mutant_rep_analysis$mutant_averaged_four_factors - mutant_rep_analysis$parent_averaged_four_factors
mutant_rep_analysis$averaged_rank6_diff <- mutant_rep_analysis$mutant_averaged_rank6 - mutant_rep_analysis$parent_averaged_rank6





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


# par(mar=c(4,8,4,4))
# plot(escape_mutant_rep_analysis$parent_averaged_rank6,
#      escape_mutant_rep_analysis$mutant_averaged_rank6,
#      xlim=c(0,6),ylim=c(0,6),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
# abline(0,1)
# dev.copy(pdf,'ave_rank6_parent_vs_mutant_rep.pdf')
# dev.off()




par(mar=c(4,8,4,4))
plot(escape_mutant_rep_analysis$parent_averaged_rank6,
     escape_mutant_rep_analysis$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
abline(0,1)
dev.copy(pdf,'ave_rank6_parent_vs_mutant_rep.pdf')
dev.off()





#FinalFig
par(mar=c(4,8,8,4))
plot(escape_mutant_rep_analysis_intraclade2$parent_averaged_rank6,
     escape_mutant_rep_analysis_intraclade2$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='red')
par(new=TRUE)
plot(escape_mutant_rep_analysis_interclade$parent_averaged_rank6,
     escape_mutant_rep_analysis_interclade$mutant_averaged_rank6,
     xlim=c(0,6),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
abline(0,1)
dev.copy(pdf,'ave_rank6_parent_vs_mutant_rep_subtypes.pdf')
dev.off()

#Number of assays used:
nrow(escape_mutant_rep_analysis_intraclade2) +
  nrow(escape_mutant_rep_analysis_interclade)
nrow(escape_mutant_rep_analysis_intraclade2)
nrow(escape_mutant_rep_analysis_interclade)

###Compare known empirical temperate to isolated mutant to escape mutant




























###Compute distance between challengers to determine minimum distance where you see a change in infection



#Create list of all unique phages used in the immunity dataset
# defending_phages <- as.character(levels(main_immunity_data$defending_phage))
# challenging_phages <- as.character(levels(main_immunity_data$challenging_phage))
# immunity_phages <- levels(factor(c(defending_phages,challenging_phages)))

immunity_phages <- union(main_immunity_data$defending_phage,main_immunity_data$challenging_phage)

#Create list of all phage triplets  

challenge_distance_min <- expand.grid(challenging_phage1 = immunity_phages,
                                      challenging_phage2 = immunity_phages,
                                      defending_phage = immunity_phages)

challenge_distance_min$challenging_phage1_challenging_phage2 <- paste(challenge_distance_min$challenging_phage1,
                                                                      "_",
                                                                      challenge_distance_min$challenging_phage2,
                                                                      sep="")

challenge_distance_min$defending_challenging_phage1 <- paste(challenge_distance_min$defending_phage,
                                                             "_",
                                                             challenge_distance_min$challenging_phage1,
                                                             sep="")

challenge_distance_min$defending_challenging_phage2 <- paste(challenge_distance_min$defending_phage,
                                                             "_",
                                                             challenge_distance_min$challenging_phage2,
                                                             sep="")

#Start adding averaged, non-redundant immunity data
#Most of the phage triplets will not be present so many rows should be dropped
#Note: if using assay_strain_averaged data, ensure that only multi-titer lysogen data is used.
#This should not be needed if using assay_averaged data is used.

temp_immunity_data <- subset(conf_assay_strain_def_chal_average,
                      conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                        conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

names(temp_immunity_data) <- paste('dc1','_',names(temp_immunity_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                temp_immunity_data,
                                by.x="defending_challenging_phage1",
                                by.y="dc1_defending_challenging")

temp_immunity_data <- subset(conf_assay_strain_def_chal_average,
                             conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                               conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')


names(temp_immunity_data) <- paste('dc2','_',names(temp_immunity_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                temp_immunity_data,
                                by.x="defending_challenging_phage2",
                                by.y="dc2_defending_challenging")



#Retrieve distance data for challenging phage pairs


#Match the genomic distance data. Phages not in actino1319 will not be matched
c1c2_genomic_distance_data <- genomic_distance_data
names(c1c2_genomic_distance_data) <- paste('c1c2','_',names(c1c2_genomic_distance_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_genomic_distance_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2")


#Match the gene distance data. Many comparisons will not be matched
c1c2_repressor_distance_data <- repressor_distance_data
names(c1c2_repressor_distance_data) <- paste('c1c2','_',names(c1c2_repressor_distance_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_repressor_distance_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)

c1c2_portal_distance_data <- portal_distance_data
names(c1c2_portal_distance_data) <- paste('c1c2','_',names(c1c2_portal_distance_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_portal_distance_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)

c1c2_recb_distance_data <- recb_distance_data
names(c1c2_recb_distance_data) <- paste('c1c2','_',names(c1c2_recb_distance_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_recb_distance_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)

c1c2_repressor336_distance_data <- repressor336_distance_data
names(c1c2_repressor336_distance_data) <- paste('c1c2','_',names(c1c2_repressor336_distance_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_repressor336_distance_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)


#TODO: merge with new portal, cas4, and endovii data


c1c2_stoperator_pwm_data <- stoperator_pwm_data
names(c1c2_stoperator_pwm_data) <- paste('c1c2','_',names(c1c2_stoperator_pwm_data),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_stoperator_pwm_data,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)




#merge defense and challenging correlations with paired challenging phages

#TODO: confirm that after merging, defending_profile_cor does not have any NA values
c1c2_defending_cor_df <- subset(defending_cor_df,select=c("phage1_phage2","defending_profile_cor"))
names(c1c2_defending_cor_df) <- paste('c1c2','_',names(c1c2_defending_cor_df),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_defending_cor_df,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)

#TODO: confirm that after merging, challenging_profile_cor does not have any NA values
c1c2_challenging_cor_df <- subset(challenging_cor_df,select=c("phage1_phage2","challenging_profile_cor"))
names(c1c2_challenging_cor_df) <- paste('c1c2','_',names(c1c2_challenging_cor_df),sep="")
challenge_distance_min <- merge(challenge_distance_min,
                                c1c2_challenging_cor_df,
                                by.x="challenging_phage1_challenging_phage2",
                                by.y="c1c2_phage1_phage2",
                                all.x=TRUE)



#Compute difference in infection profiles
challenge_distance_min$dc1dc2_averaged_infection_strength_diff <- abs(challenge_distance_min$dc1_averaged_infection_strength - challenge_distance_min$dc2_averaged_infection_strength)
challenge_distance_min$dc1dc2_averaged_four_factors_diff <- abs(challenge_distance_min$dc1_averaged_four_factors - challenge_distance_min$dc2_averaged_four_factors)
challenge_distance_min$dc1dc2_averaged_rank6_diff <- abs(challenge_distance_min$dc1_averaged_rank6 - challenge_distance_min$dc2_averaged_rank6)




#Compute comparison fields
#I can't run these through the comparison function because these comparisons involve both challenging phages and not the defending phage
challenge_distance_min$c1c2_subcluster_compare <- ifelse(challenge_distance_min$dc1_challenging_subcluster==challenge_distance_min$dc2_challenging_subcluster,
                                                     as.character(challenge_distance_min$dc1_challenging_subcluster),
                                                     "different")

challenge_distance_min$c1c2_source_compare <- ifelse(challenge_distance_min$dc1_challenging_source==challenge_distance_min$dc2_challenging_source,
                                                       as.character(challenge_distance_min$dc1_challenging_source),
                                                       "different")

challenge_distance_min$c1c2_temperate_empirical_compare <- ifelse(challenge_distance_min$dc1_challenging_cluster_a_temperate_empirical==challenge_distance_min$dc2_challenging_cluster_a_temperate_empirical,
                                                                    as.character(challenge_distance_min$dc1_challenging_cluster_a_temperate_empirical),
                                                                    "different")

challenge_distance_min$c1c2_functional_repressor_compare <- ifelse(challenge_distance_min$dc1_challenging_cluster_a_functional_repressor_predicted==challenge_distance_min$dc2_challenging_cluster_a_functional_repressor_predicted,
                                                                    as.character(challenge_distance_min$dc1_challenging_cluster_a_functional_repressor_predicted),
                                                                    "different")

challenge_distance_min$c1c2_lysogen_type_compare <- ifelse(challenge_distance_min$dc1_challenging_lysogen_type==challenge_distance_min$dc2_challenging_lysogen_type,
                                                             as.character(challenge_distance_min$dc1_challenging_lysogen_type),
                                                             "different")

challenge_distance_min$c1c2_integrase_compare <- ifelse(challenge_distance_min$dc1_challenging_pham_integrase==challenge_distance_min$dc2_challenging_pham_integrase,
                                                          as.character(challenge_distance_min$dc1_challenging_pham_integrase),
                                                          "different")

challenge_distance_min$c1c2_parb_compare <- ifelse(challenge_distance_min$dc1_challenging_pham_parb==challenge_distance_min$dc2_challenging_pham_parb,
                                                        as.character(challenge_distance_min$dc1_challenging_pham_parb),
                                                        "different")

challenge_distance_min$c1c2_repressor_hth_compare <- stringdist(as.character(challenge_distance_min$dc1_challenging_repressor_hth_domain_sequence),
                                                                 as.character(challenge_distance_min$dc2_challenging_repressor_hth_domain_sequence),
                                                                 method="hamming")

challenge_distance_min$c1c2_repressor_length_full_compare <- abs(challenge_distance_min$dc1_challenging_repressor_length_full - challenge_distance_min$dc2_challenging_repressor_length_full)
challenge_distance_min$c1c2_repressor_length_nterm_compare <- abs(challenge_distance_min$dc1_challenging_repressor_length_nterm - challenge_distance_min$dc2_challenging_repressor_length_nterm)
challenge_distance_min$c1c2_repressor_length_cterm_compare <- abs(challenge_distance_min$dc1_challenging_repressor_length_cterm - challenge_distance_min$dc2_challenging_repressor_length_cterm)
challenge_distance_min$c1c2_gene_content_clade_compare <- ifelse(challenge_distance_min$dc1_challenging_gene_content_clade==challenge_distance_min$dc2_challenging_gene_content_clade,
                                                                 as.character(challenge_distance_min$dc1_challenging_gene_content_clade),
                                                                 "different")

challenge_distance_min$c1c2_subcluster_compare <- as.factor(challenge_distance_min$c1c2_subcluster_compare)
challenge_distance_min$c1c2_source_compare <- as.factor(challenge_distance_min$c1c2_source_compare)
challenge_distance_min$c1c2_temperate_empirical_compare <- as.factor(challenge_distance_min$c1c2_temperate_empirical_compare)
challenge_distance_min$c1c2_functional_repressor_compare <- as.factor(challenge_distance_min$c1c2_functional_repressor_compare)
challenge_distance_min$c1c2_lysogen_type_compare <- as.factor(challenge_distance_min$c1c2_lysogen_type_compare)
challenge_distance_min$c1c2_integrase_compare <- as.factor(challenge_distance_min$c1c2_integrase_compare)
challenge_distance_min$c1c2_parb_compare <- as.factor(challenge_distance_min$c1c2_parb_compare)
challenge_distance_min$c1c2_gene_content_clade_compare <- as.factor(challenge_distance_min$c1c2_gene_content_clade_compare)








#Now that all data has been matched, the data is duplicated (e.g. challenging_phage1_challenging_phage2 contains l5_alma and alma_l5)
#To solve this, remove all rows in which the challenging_phage1 and challenging_phage2 are not alphabetically ordered.
#This does NOT retain self-comparisons (e.g. alma_alma)
challenge_distance_min$c1c2_alpha_ordered <- as.character(challenge_distance_min$challenging_phage1) < as.character(challenge_distance_min$challenging_phage2)

#Now retain only the unique pairwise comparisons
challenge_distance_min_alpha_ordered <- subset(challenge_distance_min,challenge_distance_min$c1c2_alpha_ordered == TRUE)




#Subset by envY and lysogen type to investigate
challenge_distance_min_ordered_envY_intY <- subset(challenge_distance_min_alpha_ordered,
                                                   challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                     challenge_distance_min_alpha_ordered$c1c2_lysogen_type_compare == 'integration')

challenge_distance_min_ordered_envY_extraY <- subset(challenge_distance_min_alpha_ordered,
                                                     challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                       challenge_distance_min_alpha_ordered$c1c2_lysogen_type_compare == 'extrachromosomal')

challenge_distance_min_ordered_envY <- subset(challenge_distance_min_alpha_ordered,
                                              challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment')


challenge_distance_min_ordered_envY_repY_empY <- subset(challenge_distance_min_alpha_ordered,
                                                        challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                          challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare == 'yes' &
                                                          challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare == 'yes')


challenge_distance_min_ordered_envDiff_repDiff_empDiff <- subset(challenge_distance_min_alpha_ordered,
                                                        challenge_distance_min_alpha_ordered$c1c2_source_compare != 'environment' |
                                                          challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare != 'yes' |
                                                          challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare != 'yes')


challenge_distance_min_ordered_envDiff <- subset(challenge_distance_min_alpha_ordered,
                                                        challenge_distance_min_alpha_ordered$c1c2_source_compare != 'environment')

challenge_distance_min_ordered_repDiff <- subset(challenge_distance_min_alpha_ordered,
                                                 challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare != 'yes')

challenge_distance_min_ordered_empDiff <- subset(challenge_distance_min_alpha_ordered,
                                                 challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare != 'yes')


challenge_distance_min_ordered_envY_repY_empY_clade2 <- subset(challenge_distance_min_alpha_ordered,
                                                        challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                          challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare == 'yes' &
                                                          challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare == 'yes' &
                                                          challenge_distance_min_alpha_ordered$dc1_challenging_gene_content_clade == 'clade2' &
                                                          challenge_distance_min_alpha_ordered$dc2_challenging_gene_content_clade == 'clade2' &
                                                          challenge_distance_min_alpha_ordered$dc2_defending_gene_content_clade == 'clade2')

challenge_distance_min_ordered_envY_repY_empY_interclade <- subset(challenge_distance_min_alpha_ordered,
                                                                   challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                                     challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare == 'yes' &
                                                                     challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare == 'yes' &
                                                                     challenge_distance_min_alpha_ordered$c1c2_gene_content_clade_compare == 'different' &
                                                                     challenge_distance_min_alpha_ordered$dc2_defending_gene_content_clade == 'clade2')


challenge_distance_min_ordered_envY_repY_empY_denvY_clade2 <- subset(challenge_distance_min_alpha_ordered,
                                                               challenge_distance_min_alpha_ordered$c1c2_source_compare == 'environment' &
                                                                 challenge_distance_min_alpha_ordered$c1c2_functional_repressor_compare == 'yes' &
                                                                 challenge_distance_min_alpha_ordered$c1c2_temperate_empirical_compare == 'yes' &
                                                                 challenge_distance_min_alpha_ordered$dc1_challenging_gene_content_clade == 'clade2' &
                                                                 challenge_distance_min_alpha_ordered$dc2_challenging_gene_content_clade == 'clade2' &
                                                                 challenge_distance_min_alpha_ordered$dc2_defending_gene_content_clade == 'clade2' &
                                                                 challenge_distance_min_alpha_ordered$dc2_defending_source == 'environment')


#Export the results
setwd("~/scratch/immunity_analysis/output/")

write.table(challenge_distance_min_ordered_envY_repY_empY,
            "challenge_distance_min_ordered_envY_repY_empY.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


#Plot the results
setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,4,4))
plot(challenge_distance_min$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min$averaged_infection_strength_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(challenge_distance_min$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min$averaged_four_factors_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(challenge_distance_min_alpha_ordered$c1c2_modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_alpha_ordered$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_ordered_envY$averaged_infection_strength_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_ordered_envY$averaged_four_factors_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envDiff_repDiff_empDiff$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_ordered_envDiff_repDiff_empDiff$averaged_four_factors_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
par(new=TRUE)
plot(challenge_distance_min_ordered_envY_repY_empY$modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$averaged_four_factors_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_min_ordered_wildtype_others_mash_vs_ave_ff_diff.pdf')
dev.off()


temp <- subset(challenge_distance_min_ordered_envY_repY_empY,challenge_distance_min_ordered_envY_repY_empY$modified_mash_distance < 0.13 &
                 challenge_distance_min_ordered_envY_repY_empY$averaged_four_factors_diff > 2.5)





par(mar=c(4,8,4,4))
plot(challenge_distance_min_alpha_ordered$c1c2_modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_alpha_ordered$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY$c1c2_modified_mash_distance,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY$c1c2_pham_pham_dissimilarity,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()





par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY$c1c2_repressor_full_mafft_dist_uncorrected,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY$c1c2_repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()





par(mar=c(4,8,4,4))
plot(challenge_distance_min_alpha_ordered$c1c2_repressor_cterm_mafft_dist_uncorrected,
     challenge_distance_min_alpha_ordered$dc1dc2_averaged_rank6_diff,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_min_alpha_ordered_c1c2_repCtermMafftUn_vs_dc1dc2_ave_rank6_diff.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY_clade2$c1c2_repressor_cterm_mafft_dist_uncorrected,
     challenge_distance_min_ordered_envY_repY_empY_clade2$dc1dc2_averaged_rank6_diff,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_envYrepYempYclade2Y_c1c2__repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY_denvY_clade2$c1c2_repressor_cterm_mafft_dist_uncorrected,
     challenge_distance_min_ordered_envY_repY_empY_denvY_clade2$dc1dc2_averaged_rank6_diff,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_envYrepYempYclade2YdenvY_repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()




par(mar=c(4,8,8,4))
plot(challenge_distance_min_ordered_envY_repY_empY_interclade$c1c2_stoperator_pwd_dist_euc,
     challenge_distance_min_ordered_envY_repY_empY_interclade$dc1dc2_averaged_rank6_diff,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
par(new=TRUE)
plot(challenge_distance_min_ordered_envY_repY_empY_clade2$c1c2_stoperator_pwd_dist_euc,
     challenge_distance_min_ordered_envY_repY_empY_clade2$dc1dc2_averaged_rank6_diff,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'challenge_distance_envYrepYempY_c1c2_stopEuc_vs_ave_rank6_diff_by_subtypes.pdf')
dev.off()





par(mar=c(4,8,4,4))
plot(challenge_distance_min_ordered_envY_repY_empY$c1c2_repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(challenge_distance_min_ordered_envY_repY_empY$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(challenge_distance_min_ordered_envDiff_repDiff_empDiff$c1c2_repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(challenge_distance_min_ordered_envDiff_repDiff_empDiff$dc1dc2_averaged_rank6_diff)),
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'.pdf')
dev.off()


#Conclusions: The best examples of large infection differences among closely-related phages:
#1. StarStuff and Jaan infecting Trixie, or StarStuff and Serenity infecting Trixie
#2. DaVinci and Gladiator infecting Et2Brutus, Mulciber, or Drake55



#TODO: save this plot once full correlation data is included
par(mar=c(4,8,4,4))
plot(challenge_distance_min_alpha_ordered$c1c2_repressor_cterm_mafft_dist_uncorrected,
     challenge_distance_min_alpha_ordered$c1c2_challenging_profile_cor,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_min_alpha_ordered_c1c2_repCtermMafftUn_vs_challenging_profile_cor.pdf')
dev.off()

#TODO: save this plot once full correlation data is included
par(mar=c(4,8,4,4))
plot(challenge_distance_min_alpha_ordered$c1c2_repressor_cterm_mafft_dist_uncorrected,
     challenge_distance_min_alpha_ordered$c1c2_defending_profile_cor,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'challenge_distance_min_alpha_ordered_c1c2_repCtermMafftUn_vs_defending_profile_cor.pdf')
dev.off()





###Compute distance between challengers to determine minimum distance where you see a change in infection above





























###Compute distance between defendings to determine minimum distance where you see a change in infection

#Create list of all unique phages used in the immunity dataset
immunity_phages <- union(main_immunity_data$defending_phage,main_immunity_data$challenging_phage)

#Create list of all phage triplets  


defend_distance_min <- expand.grid(defending_phage1 = immunity_phages,
                                   defending_phage2 = immunity_phages,
                                   challenging_phage = immunity_phages)

defend_distance_min$defending_phage1_defending_phage2 <- paste(defend_distance_min$defending_phage1,
                                                               "_",
                                                               defend_distance_min$defending_phage2,
                                                               sep="")

defend_distance_min$defending_phage1_challenging <- paste(defend_distance_min$defending_phage1,
                                                          "_",
                                                          defend_distance_min$challenging_phage,
                                                          sep="")

defend_distance_min$defending_phage2_challenging <- paste(defend_distance_min$defending_phage2,
                                                          "_",
                                                          defend_distance_min$challenging_phage,
                                                          sep="")

#Start adding averaged, non-redundant immunity data. Most of the phage triplets will not be present so many rows should be dropped
#Note: if using assay_strain_averaged data, ensure that only multi-titer lysogen data is used.
#This should not be needed if using assay_averaged data is used.

temp_immunity_data <- subset(conf_assay_strain_def_chal_average,
                             conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                               conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

names(temp_immunity_data) <- paste('d1c','_',names(temp_immunity_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             temp_immunity_data,
                             by.x="defending_phage1_challenging",
                             by.y="d1c_defending_challenging")

temp_immunity_data <- subset(conf_assay_strain_def_chal_average,
                             conf_assay_strain_def_chal_average$strain_type == 'lysogen' & 
                               conf_assay_strain_def_chal_average$assay_type == 'multiple_titer')

names(temp_immunity_data) <- paste('d2c','_',names(temp_immunity_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             temp_immunity_data,
                             by.x="defending_phage2_challenging",
                             by.y="d2c_defending_challenging")






#Retrieve distance data for defending phage pairs


#Match the genomic distance data. Phages not in actino1319 will not be matched
d1d2_genomic_distance_data <- genomic_distance_data
names(d1d2_genomic_distance_data) <- paste('d1d2','_',names(d1d2_genomic_distance_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_genomic_distance_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2")


#Match the gene distance data. Many comparisons will not be matched
d1d2_repressor_distance_data <- repressor_distance_data
names(d1d2_repressor_distance_data) <- paste('d1d2','_',names(d1d2_repressor_distance_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_repressor_distance_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2",
                             all.x=TRUE)

d1d2_portal_distance_data <- portal_distance_data
names(d1d2_portal_distance_data) <- paste('d1d2','_',names(d1d2_portal_distance_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_portal_distance_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2",
                             all.x=TRUE)

d1d2_recb_distance_data <- recb_distance_data
names(d1d2_recb_distance_data) <- paste('d1d2','_',names(d1d2_recb_distance_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_recb_distance_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2",
                             all.x=TRUE)

d1d2_repressor336_distance_data <- repressor336_distance_data
names(d1d2_repressor336_distance_data) <- paste('d1d2','_',names(d1d2_repressor336_distance_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_repressor336_distance_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2",
                             all.x=TRUE)

#TODO: merge with new portal, cas4, endovii data

d1d2_stoperator_pwm_data <- stoperator_pwm_data
names(d1d2_stoperator_pwm_data) <- paste('d1d2','_',names(d1d2_stoperator_pwm_data),sep="")
defend_distance_min <- merge(defend_distance_min,
                             d1d2_stoperator_pwm_data,
                             by.x="defending_phage1_defending_phage2",
                             by.y="d1d2_phage1_phage2",
                             all.x=TRUE)



#merge defense and challenging correlations with paired challenging phages

#TODO: confirm that after merging, defending_profile_cor does not have any NA values
d1d2_defending_cor_df <- subset(defending_cor_df,select=c("phage1_phage2","defending_profile_cor"))
names(d1d2_defending_cor_df) <- paste('d1d2','_',names(d1d2_defending_cor_df),sep="")
defend_distance_min <- merge(defend_distance_min,
                                d1d2_defending_cor_df,
                                by.x="defending_phage1_defending_phage2",
                                by.y="d1d2_phage1_phage2",
                                all.x=TRUE)


#TODO: confirm that after merging, challenging_profile_cor does not have any NA values
d1d2_challenging_cor_df <- subset(challenging_cor_df,select=c("phage1_phage2","challenging_profile_cor"))
names(d1d2_challenging_cor_df) <- paste('d1d2','_',names(d1d2_challenging_cor_df),sep="")
defend_distance_min <- merge(defend_distance_min,
                                d1d2_challenging_cor_df,
                                by.x="defending_phage1_defending_phage2",
                                by.y="d1d2_phage1_phage2",
                                all.x=TRUE)



#Compute difference in infection profiles
defend_distance_min$d1cd2c_averaged_infection_strength_diff <- abs(defend_distance_min$d1c_averaged_infection_strength - defend_distance_min$d2c_averaged_infection_strength)
defend_distance_min$d1cd2c_averaged_four_factors_diff <- abs(defend_distance_min$d1c_averaged_four_factors - defend_distance_min$d2c_averaged_four_factors)
defend_distance_min$d1cd2c_averaged_rank6_diff <- abs(defend_distance_min$d1c_averaged_rank6 - defend_distance_min$d2c_averaged_rank6)



#Compute comparison fields
defend_distance_min$d1d2_subcluster_compare <- ifelse(defend_distance_min$d1c_defending_subcluster==defend_distance_min$d2c_defending_subcluster,
                                                  as.character(defend_distance_min$d1c_defending_subcluster),
                                                  "different")

defend_distance_min$d1d2_source_compare <- ifelse(defend_distance_min$d1c_defending_source==defend_distance_min$d2c_defending_source,
                                                      as.character(defend_distance_min$d1c_defending_source),
                                                      "different")

defend_distance_min$d1d2_temperate_empirical_compare <- ifelse(defend_distance_min$d1c_defending_cluster_a_temperate_empirical==defend_distance_min$d2c_defending_cluster_a_temperate_empirical,
                                                                   as.character(defend_distance_min$d1c_defending_cluster_a_temperate_empirical),
                                                                   "different")

defend_distance_min$d1d2_functional_repressor_compare <- ifelse(defend_distance_min$d1c_defending_cluster_a_functional_repressor_predicted==defend_distance_min$d2c_defending_cluster_a_functional_repressor_predicted,
                                                                   as.character(defend_distance_min$d1c_defending_cluster_a_functional_repressor_predicted),
                                                                   "different")

defend_distance_min$d1d2_lysogen_type_compare <- ifelse(defend_distance_min$d1c_defending_lysogen_type==defend_distance_min$d2c_defending_lysogen_type,
                                                            as.character(defend_distance_min$d1c_defending_lysogen_type),
                                                            "different")

defend_distance_min$d1d2_integrase_compare <- ifelse(defend_distance_min$d1c_defending_pham_integrase==defend_distance_min$d2c_defending_pham_integrase,
                                                         as.character(defend_distance_min$d1c_defending_pham_integrase),
                                                         "different")

defend_distance_min$d1d2_parb_compare <- ifelse(defend_distance_min$d1c_defending_pham_parb==defend_distance_min$d2c_defending_pham_parb,
                                                     as.character(defend_distance_min$d1c_defending_pham_parb),
                                                     "different")

defend_distance_min$d1d2_repressor_hth_compare <- stringdist(as.character(defend_distance_min$d1c_defending_repressor_hth_domain_sequence),
                                                                as.character(defend_distance_min$d2c_defending_repressor_hth_domain_sequence),
                                                                method="hamming")

defend_distance_min$d1d2_repressor_length_full_compare <- abs(defend_distance_min$d1c_defending_repressor_length_full - defend_distance_min$d2c_defending_repressor_length_full)
defend_distance_min$d1d2_repressor_length_nterm_compare <- abs(defend_distance_min$d1c_defending_repressor_length_nterm - defend_distance_min$d2c_defending_repressor_length_nterm)
defend_distance_min$d1d2_repressor_length_cterm_compare <- abs(defend_distance_min$d1c_defending_repressor_length_cterm - defend_distance_min$d2c_defending_repressor_length_cterm)

defend_distance_min$d1d2_gene_content_clade_compare <- ifelse(defend_distance_min$d1c_defending_gene_content_clade==defend_distance_min$d2c_defending_gene_content_clade,
                                                                 as.character(defend_distance_min$d1c_defending_gene_content_clade),
                                                                 "different")


defend_distance_min$d1d2_subcluster_compare <- as.factor(defend_distance_min$d1d2_subcluster_compare)
defend_distance_min$d1d2_source_compare <- as.factor(defend_distance_min$d1d2_source_compare)
defend_distance_min$d1d2_temperate_empirical_compare <- as.factor(defend_distance_min$d1d2_temperate_empirical_compare)
defend_distance_min$d1d2_functional_repressor_compare <- as.factor(defend_distance_min$d1d2_functional_repressor_compare)
defend_distance_min$d1d2_lysogen_type_compare <- as.factor(defend_distance_min$d1d2_lysogen_type_compare)
defend_distance_min$d1d2_integrase_compare <- as.factor(defend_distance_min$d1d2_integrase_compare)
defend_distance_min$d1d2_parb_compare <- as.factor(defend_distance_min$d1d2_parb_compare)
defend_distance_min$d1d2_gene_content_clade_compare <- as.factor(defend_distance_min$d1d2_gene_content_clade_compare)






#Now that all data has been matched, the data is duplicated (e.g. defending_phage1_defending_phage2 contains l5_alma and alma_l5)
#To solve this, remove all rows in which the defending_phage1 and defending_phage2 are not alphabetically ordered.
#This does NOT retain self-comparisons (e.g. alma_alma)
defend_distance_min$d1d2_alpha_ordered <- as.character(defend_distance_min$defending_phage1) < as.character(defend_distance_min$defending_phage2)

#Now retain only the unique pairwise comparisons
defend_distance_min_alpha_ordered <- subset(defend_distance_min,defend_distance_min$d1d2_alpha_ordered == TRUE)



#Subset by envY and lysogen type to investigate

defend_distance_min_ordered_envY <- subset(defend_distance_min_alpha_ordered,
                                              defend_distance_min_alpha_ordered$d1d2_source_compare == 'environment')

defend_distance_min_ordered_envN <- subset(defend_distance_min_alpha_ordered,
                                           defend_distance_min_alpha_ordered$d1d2_source_compare != 'environment')




defend_distance_min_ordered_envY_clade2 <- subset(defend_distance_min_alpha_ordered,
                                           defend_distance_min_alpha_ordered$d1d2_source_compare == 'environment' &
                                             defend_distance_min_alpha_ordered$d1c_defending_gene_content_clade == 'clade2' &
                                             defend_distance_min_alpha_ordered$d2c_defending_gene_content_clade == 'clade2' &
                                             defend_distance_min_alpha_ordered$d2c_challenging_gene_content_clade == 'clade2')


defend_distance_min_ordered_envY_cenvYcrepYctempY_clade2 <- subset(defend_distance_min_alpha_ordered,
                                                                   defend_distance_min_alpha_ordered$d1d2_source_compare == 'environment' &
                                                                     defend_distance_min_alpha_ordered$d1c_defending_gene_content_clade == 'clade2' &
                                                                     defend_distance_min_alpha_ordered$d2c_defending_gene_content_clade == 'clade2' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_gene_content_clade == 'clade2' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_source == 'environment' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_cluster_a_functional_repressor_predicted == 'yes' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_cluster_a_temperate_empirical == 'yes')


defend_distance_min_ordered_envY_cenvYcrepYctempY_interclade <- subset(defend_distance_min_alpha_ordered,
                                                                   defend_distance_min_alpha_ordered$d1d2_source_compare == 'environment' &
                                                                     defend_distance_min_alpha_ordered$d1d2_gene_content_clade_compare == 'different' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_gene_content_clade == 'clade2' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_source == 'environment' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_cluster_a_functional_repressor_predicted == 'yes' &
                                                                     defend_distance_min_alpha_ordered$d2c_challenging_cluster_a_temperate_empirical == 'yes')







#Plot the results
setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY$modified_mash_distance,
     as.numeric(as.character(defend_distance_min_ordered_envY$averaged_infection_strength_diff)),
     xlim=c(0,0.4),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY$modified_mash_distance,
     as.numeric(as.character(defend_distance_min_ordered_envY$averaged_four_factors_diff)),
     xlim=c(0,0.4),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(defend_distance_min_ordered_envN$modified_mash_distance,
     as.numeric(as.character(defend_distance_min_ordered_envN$averaged_four_factors_diff)),
     xlim=c(0,0.4),ylim=c(0,8),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'defend_distance_min_envY_envN_ordered_envY_mash_vs_ave_ff_diff.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY$d1d2_repressor_cterm_mafft_dist_uncorrected,
     as.numeric(as.character(defend_distance_min_ordered_envY$d1cd2c_averaged_rank6_diff)),
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'defend_distance_min_ordered_envY_d1d2_repCtermMafftUn_vs_d1cd2c_ave_rank6_diff.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY_clade2$d1d2_repressor_cterm_mafft_dist_uncorrected,
     defend_distance_min_ordered_envY_clade2$d1cd2c_averaged_rank6_diff,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'defend_distance_min_envY_clade2_d1d2_repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()



par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY_cenvYcrepYctempY_clade2$d1d2_repressor_cterm_mafft_dist_uncorrected,
     defend_distance_min_ordered_envY_cenvYcrepYctempY_clade2$d1cd2c_averaged_rank6_diff,
     xlim=c(0,60),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'defend_distance_min_envY_cenvYcrepYctempY_clade2_repCtermMafftUn_vs_ave_rank6_diff.pdf')
dev.off()






par(mar=c(4,8,8,4))
plot(defend_distance_min_ordered_envY_cenvYcrepYctempY_interclade$d1d2_stoperator_pwd_dist_euc,
     defend_distance_min_ordered_envY_cenvYcrepYctempY_interclade$d1cd2c_averaged_rank6_diff,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='grey')
par(new=TRUE)
plot(defend_distance_min_ordered_envY_cenvYcrepYctempY_clade2$d1d2_stoperator_pwd_dist_euc,
     defend_distance_min_ordered_envY_cenvYcrepYctempY_clade2$d1cd2c_averaged_rank6_diff,
     xlim=c(0,5),ylim=c(0,6),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col='orange')
dev.copy(pdf,'defend_distance_min_envY_cenvYcrepYctempY_d1d2_stopEuc_vs_ave_rank6_diff_by_subtypes.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(defend_distance_min_ordered_envY$d1d2_pham_pham_dissimilarity,
     as.numeric(as.character(defend_distance_min_ordered_envY$d1cd2c_averaged_rank6_diff)),
     xlim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(defend_distance_min_ordered_envN$d1d2_pham_pham_dissimilarity,
     defend_distance_min_ordered_envN$d1cd2c_averaged_rank6_diff,
     xlim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'.pdf')
dev.off()








#TODO: save plot once all correlation data is included
par(mar=c(4,8,4,4))
plot(defend_distance_min_alpha_ordered$d1d2_repressor_cterm_mafft_dist_uncorrected,
     defend_distance_min_alpha_ordered$d1d2_defending_profile_cor,
     xlim=c(0,70),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'defend_distance_min_ordered_d1d2_repCtermMafftUn_vs_defending_profile_cor.pdf')
dev.off()

#TODO: save plot once all correlation data is included
par(mar=c(4,8,4,4))
plot(defend_distance_min_alpha_ordered$d1d2_repressor_cterm_mafft_dist_uncorrected,
     defend_distance_min_alpha_ordered$d1d2_challenging_profile_cor,
     xlim=c(0,70),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
dev.copy(pdf,'defend_distance_min_ordered_d1d2_repCtermMafftUn_vs_challenging_profile_cor.pdf')
dev.off()





###Compute distance between defendings to determine minimum distance where you see a change in infection above































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
distance_metrics <- merge(distance_metrics,repressor_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,portal_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
distance_metrics <- merge(distance_metrics,recb_distance_data,by.x="phage1_phage2",by.y="phage1_phage2",all.x=TRUE)
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


#TODO: merge with stoperator site prediction tally data?


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




# group1 <- factor(c("A1"))
# group2 <- factor(c("A5","A8","A18","A3","A19","A7","A10"))
# group3 <- factor(c("A4"))
# group4 <- factor(c("A13","A15","A12","A9","A14","A2","A11","A6","A16","A17"))


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






#Plot # Cluster A phages and their GCD diversity


#When computing genomic diversity, do not include the escape mutants
clusterA_data <- subset(distance_metrics,
                        distance_metrics$cluster_compare == "A" &
                          distance_metrics$source_compare == "environment")



#Subcluster factored list
clusterA_subcluster_list2 <- factor(c("A17",
                                      "A18",
                                      "A13",
                                      "A19",
                                      "A14",
                                      "A12",
                                      "A16",
                                      "A11",
                                      "A15",
                                      "A8",
                                      "A9",
                                      "A6",
                                      "A4",
                                      "A5",
                                      "A10",
                                      "A7",
                                      "A1",
                                      "A3",
                                      "A2"))

clusterA_data$phage1_subcluster <- factor(clusterA_data$phage1_subcluster,clusterA_subcluster_list2)
clusterA_data$phage2_subcluster <- factor(clusterA_data$phage2_subcluster,clusterA_subcluster_list2)



clusterA_subcluster_colorlist2 <- c("orange","grey",#"A17"
                                    "cyan","grey",#"A18"
                                    "orange","grey",#"A13"
                                    "green","grey",#"A19"
                                    "orange","grey",#"A14"
                                    "orange","grey",#"A12"
                                    "orange","grey",#"A16"
                                    "orange","grey",#"A11"
                                    "orange","grey",#"A15"
                                    "cyan","grey",#"A8"
                                    "orange","grey",#"A9"
                                    "orange","grey",#"A6"
                                    "green","grey",#"A4"
                                    "cyan","grey",#"A5"
                                    "green","grey",#"A10"
                                    "green","grey",#"A7"
                                    "red","grey",#"A1"
                                    "green","grey",#"A3"
                                    "orange","grey"#"A2"
)



# par(mar=c(10,8,4,4))
# stripchart(clusterA_data$pham_pham_dissimilarity ~ clusterA_data$subcluster_compare2*clusterA_data$phage1_subcluster,
#            vertical=TRUE,method="jitter",pch=19,cex=0.5,ylim=c(0,1),las=2,cex.axis=1,ann=FALSE,main=NULL,col=c("blue","green"))
# dev.copy(pdf,'gcd_distribution_by_subcluster.pdf')
# dev.off()

setwd("~/scratch/immunity_analysis/output/")
par(mar=c(10,8,4,4))
stripchart(clusterA_data$pham_pham_dissimilarity ~ clusterA_data$subcluster_compare2*clusterA_data$phage1_subcluster,
           vertical=TRUE,method="jitter",pch=19,cex=0.5,ylim=c(0,1),las=2,cex.axis=1,ann=FALSE,main=NULL,col=clusterA_subcluster_colorlist2)
dev.copy(pdf,'gcd_distribution_by_subcluster2.pdf')
dev.off()






#Gene content clade factored list
clusterA_gene_content_clade_list <- factor(c("clade1",
                                      "clade4",
                                      "clade3",
                                      "clade2"))

clusterA_data$phage1_gene_content_clade <- factor(clusterA_data$phage1_gene_content_clade,clusterA_gene_content_clade_list)
clusterA_data$phage2_gene_content_clade <- factor(clusterA_data$phage2_gene_content_clade,clusterA_gene_content_clade_list)



clusterA_gene_content_clade_colorlist <- c("red","grey",#"clade1"
                                           "cyan","grey",#"clade4"
                                           "green","grey",#"clade3"
                                           "orange","grey"#"clade2"
)


par(mar=c(10,8,4,4))
stripchart(clusterA_data$pham_pham_dissimilarity ~ clusterA_data$gene_content_clade_compare2*clusterA_data$phage1_gene_content_clade,
           vertical=TRUE,method="jitter",pch=19,cex=0.5,ylim=c(0,1),las=2,cex.axis=1,ann=FALSE,main=NULL,col=clusterA_gene_content_clade_colorlist)
dev.copy(pdf,'gcd_distribution_by_gene_content_clade.pdf')
dev.off()








#Phage frequency

#By subcluster
clusterA_metadata <- subset(phage_metadata,
                            phage_metadata$cluster == 'A' &
                              phage_metadata$source == 'environment')


clusterA_metadata$subcluster <- factor(clusterA_metadata$subcluster,clusterA_subcluster_list2)

par(mar=c(10,8,4,4))
barplot(summary(clusterA_metadata$subcluster),
        ylim=c(0,100),col="black",las=2)
dev.copy(pdf,'subcluster_abundance.pdf')
dev.off()




#By gene content clade
clusterA_metadata$gene_content_clade <- factor(clusterA_metadata$gene_content_clade,clusterA_gene_content_clade_list)

par(mar=c(10,8,4,4))
barplot(summary(clusterA_metadata$gene_content_clade),
        ylim=c(0,150),col="black",las=2)
dev.copy(pdf,'clade_abundance.pdf')
dev.off()


















#Plot relationship between repressors, stoperators, and genomes


par(mar=c(4,8,4,4))
plot(distance_metrics$modified_mash_distance,
     distance_metrics$pham_pham_dissimilarity,
     xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(4,8,4,4))
plot(distance_metrics$modified_mash_distance,
     distance_metrics$stoperator_pwd_dist_euc,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$stoperator_pwd_dist_euc,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_stopEuc.pdf')
dev.off()
#TODO: review outliers in gcd_vs_stoperator 









par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$repressor_full_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$repressor_nterm_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_repCtermMafUn.pdf')
dev.off()









par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_nterm_mafft_dist_uncorrected,
     distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1)





par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$cas4_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_cas4MafftUn.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$endovii_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_endoviiMafftUn.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$dnapol_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_dnapolMafftUn.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$portal_mafft_dist_uncorrected,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'gcd_vs_portalMafftUn.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$portal_muscle_bionj_distances,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$portal_prank_phyml_distances,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$recb_muscle_bionj_distances,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
















#Keep only Mycobacteriophage data and drop Gordonia phages
clusterA_data <- subset(distance_metrics,
                        distance_metrics$cluster_compare == "A" &
                          distance_metrics$source_compare == "environment" &
                          distance_metrics$phage1_host == "Mycobacterium" &
                          distance_metrics$phage2_host == "Mycobacterium")


clusterA_clade1 <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "clade1")
clusterA_clade2 <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "clade2")
clusterA_clade3 <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "clade3")
clusterA_clade4 <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "clade4")
clusterA_cladeDiff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different")

clusterA_clade1_diff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different" &
                                 (clusterA_data$phage1_gene_content_clade == 'clade1' |
                                    clusterA_data$phage2_gene_content_clade == 'clade1'))

clusterA_clade2_diff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different" &
                                 (clusterA_data$phage1_gene_content_clade == 'clade2' |
                                    clusterA_data$phage2_gene_content_clade == 'clade2'))

clusterA_clade3_diff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different" &
                                 (clusterA_data$phage1_gene_content_clade == 'clade3' |
                                    clusterA_data$phage2_gene_content_clade == 'clade3'))

clusterA_clade4_diff <- subset(clusterA_data,clusterA_data$gene_content_clade_compare == "different" &
                                 (clusterA_data$phage1_gene_content_clade == 'clade4' |
                                    clusterA_data$phage2_gene_content_clade == 'clade4'))









#Plots
setwd("~/scratch/immunity_analysis/output/")










#GCD vs portal distances

par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$portal_mafft_dist_uncorrected,
     xlim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

#Keep y axis at 70 to directly compare to repressor plots
#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_portal_clade1.pdf')
dev.off()



#clade2
par(mar=c(4,8,4,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
dev.copy(pdf,'gcd_vs_portal_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_portal_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$portal_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_portal_clade4.pdf')
dev.off()




#






#GCD vs dna pol distances

par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

#Keep y axis at 70 to directly compare to repressor plots
#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_dnapol_clade1.pdf')
dev.off()



#clade2
par(mar=c(4,8,4,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
dev.copy(pdf,'gcd_vs_dnapol_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_dnapol_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$dnapol_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_dnapol_clade4.pdf')
dev.off()




#









#GCD vs endovii distances

par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

#Keep y axis at 70 to directly compare to repressor plots
#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_endovii_clade1.pdf')
dev.off()



#clade2
par(mar=c(4,8,4,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
dev.copy(pdf,'gcd_vs_endovii_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_endovii_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$endovii_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_endovii_clade4.pdf')
dev.off()




#










#GCD vs cas4 distances

par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

#Keep y axis at 70 to directly compare to repressor plots
#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_cas4_clade1.pdf')
dev.off()



#clade2
par(mar=c(4,8,4,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
dev.copy(pdf,'gcd_vs_cas4_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_cas4_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$cas4_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_cas4_clade4.pdf')
dev.off()




#


























#GCD vs repressor distances
par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)


par(mar=c(4,8,4,4))
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
par(mar=c(4,8,4,4))
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(clusterA_cladeDiff$pham_pham_dissimilarity,
     clusterA_cladeDiff$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")




#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_repFullMafUn_clade1.pdf')
dev.off()



#clade2
#FinalFig
par(mar=c(4,8,8,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_repFullMafUn_clade2.pdf')
dev.off()

#Number of comparisons:
nrow(clusterA_clade2_diff) +
  nrow(clusterA_clade2)



#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_repFullMafUn_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_repFullMafUn_clade4.pdf')
dev.off()

#













#GCD vs stoperator motif distances
par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$pham_pham_dissimilarity,
     clusterA_clade1_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$pham_pham_dissimilarity,
     clusterA_clade1$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_stopEuc_clade1.pdf')
dev.off()



#clade2
#FinalFig
par(mar=c(4,8,8,4))
plot(clusterA_clade2_diff$pham_pham_dissimilarity,
     clusterA_clade2_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$pham_pham_dissimilarity,
     clusterA_clade2$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'gcd_vs_stopEuc_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$pham_pham_dissimilarity,
     clusterA_clade3_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$pham_pham_dissimilarity,
     clusterA_clade3$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'gcd_vs_stopEuc_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$pham_pham_dissimilarity,
     clusterA_clade4_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$pham_pham_dissimilarity,
     clusterA_clade4$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'gcd_vs_stopEuc_clade4.pdf')
dev.off()








#repressor vs stoperator motif distances
par(mar=c(4,8,4,4))
plot(clusterA_data$repressor_full_mafft_dist_uncorrected,
     clusterA_data$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




#clade1
par(mar=c(4,8,4,4))
plot(clusterA_clade1_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_clade1_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade1$repressor_full_mafft_dist_uncorrected,
     clusterA_clade1$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'repMafftUn_vs_stopEuc_clade1.pdf')
dev.off()



#clade2
#FinalFig
par(mar=c(4,8,8,4))
plot(clusterA_clade2_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_clade2_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$repressor_full_mafft_dist_uncorrected,
     clusterA_clade2$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'repMafftUn_vs_stopEuc_clade2.pdf')
dev.off()


#clade3
par(mar=c(4,8,4,4))
plot(clusterA_clade3_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_clade3_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade3$repressor_full_mafft_dist_uncorrected,
     clusterA_clade3$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
dev.copy(pdf,'repMafftUn_vs_stopEuc_clade3.pdf')
dev.off()


#clade4
par(mar=c(4,8,4,4))
plot(clusterA_clade4_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_clade4_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade4$repressor_full_mafft_dist_uncorrected,
     clusterA_clade4$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
dev.copy(pdf,'repMafftUn_vs_stopEuc_clade4.pdf')
dev.off()




















#Mash vs GCD
#FinalFig
par(mar=c(4,8,8,4))
plot(clusterA_clade2_diff$modified_mash_distance,
     clusterA_clade2_diff$pham_pham_dissimilarity,
     xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$modified_mash_distance,
     clusterA_clade2$pham_pham_dissimilarity,
     xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
dev.copy(pdf,'mash_vs_gcd_clade2.pdf')
dev.off()










#RepNterm vs RepCterm
#clade2
#FinalFig
par(mar=c(4,8,8,4))
plot(clusterA_clade2_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_clade2_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_clade2$repressor_nterm_mafft_dist_uncorrected,
     clusterA_clade2$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1)
dev.copy(pdf,'repNtermMafftUn_vs_repCtermMaffUn_clade2.pdf')
dev.off()

















par(mar=c(4,8,4,4))
plot(clusterA_data$repressor_full_mafft_dist_uncorrected,
     clusterA_data$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")

par(mar=c(4,8,4,4))
plot(clusterA_data$recb_muscle_bionj_distances,
     clusterA_data$stoperator_pwd_dist_euc,
     ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")

par(mar=c(4,8,4,4))
plot(clusterA_data$portal_muscle_bionj_distances,
     clusterA_data$stoperator_pwd_dist_euc,
     ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")


par(mar=c(4,8,4,4))
plot(clusterA_data$repressor_hth_compare,
     clusterA_data$stoperator_pwd_dist_euc,
     ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")



temp3 <- subset(clusterA_data,select=c("repressor_full_mafft_dist_uncorrected","stoperator_pwd_dist_euc"))
temp3 <- temp3[complete.cases(temp3),]
cor(temp3$repressor_full_mafft_dist_uncorrected,temp3$stoperator_pwd_dist_euc)

#use pearson for continuous variables
cor(clusterA_data$repressor_full_mafft_dist_uncorrected,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_repFull <- lm(stoperator_pwd_dist_euc ~
                   repressor_full_mafft_dist_uncorrected,
                 data = clusterA_data)
summary(lm_repFull)

plot(clusterA_data$repressor_full_mafft_dist_uncorrected,
     clusterA_data$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(lm_repFull)



cor(clusterA_data$repressor_nterm_mafft_dist_uncorrected,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_repNterm <- lm(repressor_nterm_mafft_dist_uncorrected ~
                    stoperator_pwd_dist_euc,
                  data = clusterA_data)
summary(lm_repNterm)


cor(clusterA_data$repressor_cterm_mafft_dist_uncorrected,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_repCterm <- lm(repressor_cterm_mafft_dist_uncorrected ~
                    stoperator_pwd_dist_euc,
                  data = clusterA_data)
summary(lm_repCterm)


cor(clusterA_data$repressor_hth_compare,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_repHth <- lm(repressor_hth_compare ~
                  stoperator_pwd_dist_euc,
                data = clusterA_data)
summary(lm_repHth)

cor(clusterA_data$recb_muscle_bionj_distances,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_Recb <- lm(recb_muscle_bionj_distances ~
                stoperator_pwd_dist_euc,
              data = clusterA_data)
summary(lm_Recb)

cor(clusterA_data$portal_muscle_bionj_distances,
    clusterA_data$stoperator_pwd_dist_euc,
    use="complete.obs",
    method="pearson")
lm_portal <- lm(portal_muscle_bionj_distances ~
                  stoperator_pwd_dist_euc,
                data = clusterA_data)
summary(lm_portal)



lm_group1_repFull <- lm(repressor_full_mafft_dist_uncorrected ~
                          stoperator_pwd_dist_euc,
                        data = clusterA_group1)
summary(lm_group1_repFull)

lm_group2_repFull <- lm(repressor_full_mafft_dist_uncorrected ~
                          stoperator_pwd_dist_euc,
                        data = clusterA_group2)
summary(lm_group2_repFull)

lm_group3_repFull <- lm(repressor_full_mafft_dist_uncorrected ~
                          stoperator_pwd_dist_euc,
                        data = clusterA_group3)
summary(lm_group3_repFull)

lm_group4_repFull <- lm(repressor_full_mafft_dist_uncorrected ~
                          stoperator_pwd_dist_euc,
                        data = clusterA_group4)
summary(lm_group4_repFull)


#stoperator pwm
par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")

temp <- subset(clusterA_data,clusterA_data$pham_pham_dissimilarity > 0.6 & clusterA_data$stoperator_pwd_dist_euc < 0.8)
temp2 <- subset(clusterA_data,clusterA_data$pham_pham_dissimilarity < 0.2 & clusterA_data$stoperator_pwd_dist_euc > 2.2)


#group1
par(mar=c(4,8,4,4))
plot(clusterA_group1_diff$pham_pham_dissimilarity,
     clusterA_group1_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group1$pham_pham_dissimilarity,
     clusterA_group1$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

#group2
par(mar=c(4,8,4,4))
plot(clusterA_group2_diff$pham_pham_dissimilarity,
     clusterA_group2_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group2$pham_pham_dissimilarity,
     clusterA_group2$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")

#group3
par(mar=c(4,8,4,4))
plot(clusterA_group3_diff$pham_pham_dissimilarity,
     clusterA_group3_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group3$pham_pham_dissimilarity,
     clusterA_group3$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

#group4
par(mar=c(4,8,4,4))
plot(clusterA_group4_diff$pham_pham_dissimilarity,
     clusterA_group4_diff$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group4$pham_pham_dissimilarity,
     clusterA_group4$stoperator_pwd_dist_euc,
     xlim=c(0,1),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")







#group1
par(mar=c(4,8,4,4))
plot(clusterA_group1_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_group1_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group1$repressor_full_mafft_dist_uncorrected,
     clusterA_group1$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

#group2
par(mar=c(4,8,4,4))
plot(clusterA_group2_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_group2_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group2$repressor_full_mafft_dist_uncorrected,
     clusterA_group2$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")

#group3
par(mar=c(4,8,4,4))
plot(clusterA_group3_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_group3_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group3$repressor_full_mafft_dist_uncorrected,
     clusterA_group3$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

#group4
par(mar=c(4,8,4,4))
plot(clusterA_group4_diff$repressor_full_mafft_dist_uncorrected,
     clusterA_group4_diff$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group4$repressor_full_mafft_dist_uncorrected,
     clusterA_group4$stoperator_pwd_dist_euc,
     xlim=c(0,70),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
















par(mar=c(10,8,4,4))
stripchart(clusterA_data$pham_pham_dissimilarity ~ clusterA_data$group_compare2*clusterA_data$phage1_group,
           vertical=TRUE,method="jitter",pch=19,cex=0.01,ylim=c(0,1),las=2,cex.axis=1,ann=FALSE,main=NULL,col=c("blue","green"))
dev.copy(pdf,'gcd_distribution_by_group.pdf')
dev.off()


par(mar=c(10,8,4,4))

hist(clusterA_group1$pham_pham_dissimilarity,
     xlim=c(0,1),col='black',ann=FALSE,main=NULL,las=1,breaks=100,ylim=c(0,400))
dev.copy(pdf,'clusterA_group1_gcd_hist.pdf')
dev.off()

hist(clusterA_group2$pham_pham_dissimilarity,
     xlim=c(0,1),col='black',ann=FALSE,main=NULL,las=1,breaks=100,ylim=c(0,400))

hist(clusterA_group3$pham_pham_dissimilarity,
     xlim=c(0,1),col='black',ann=FALSE,main=NULL,las=1,breaks=100,ylim=c(0,400))

hist(clusterA_group4$pham_pham_dissimilarity,
     xlim=c(0,1),col='black',ann=FALSE,main=NULL,las=1,breaks=100,ylim=c(0,400))















clusterA_data$subcluster1_subcluster2 <- paste(clusterA_data$phage1_subcluster,"_",clusterA_data$phage2_subcluster,sep="")
clusterA_data$subcluster1_subcluster2 <- as.factor(clusterA_data$subcluster1_subcluster2)

clusterA_subcluster_data <- subset(clusterA_data,
                                   clusterA_data$subcluster_compare != "different")

clusterA_subcluster_data$subcluster_compare <- factor(clusterA_subcluster_data$subcluster_compare)


subclusterA1_data <- subset(clusterA_data,
                            clusterA_data$subcluster_compare == "A1")


subclusterA1_data <- subset(clusterA_data,
                            clusterA_data$subcluster_compare == "A1")


A3_A4_data <- subset(clusterA_data,
                     clusterA_data$subcluster1_subcluster2 == 'A3_A4' |
                       clusterA_data$subcluster1_subcluster2 == 'A4_A3')

par(mar=c(4,8,4,4))
plot(A3_A4_data$pham_pham_dissimilarity,
     A3_A4_data$repressor_full_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(4,8,4,4))
plot(clusterA_data$pham_pham_dissimilarity,
     clusterA_data$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(clusterA_subcluster_data$pham_pham_dissimilarity,
     clusterA_subcluster_data$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(subclusterA1_data$pham_pham_dissimilarity,
     subclusterA1_data$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,1),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)






par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$stoperator_pwd_dist_pearson,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     distance_metrics$stoperator_pwd_dist_euc,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'repCtermMafUn_vs_stopEuc.pdf')
dev.off()


#TODO: review outliers in repressorCterm_vs_stoperator 

par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     distance_metrics$stoperator_pwd_dist_pearson,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)


par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_nterm_mafft_dist_uncorrected,
     distance_metrics$stoperator_pwd_dist_euc,
     pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)









#Repressor similarity
par(mar=c(4,8,4,4))
plot(clusterA_group1$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group1$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1)
dev.copy(pdf,'clusterA_group1_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group2$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group2$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1)
dev.copy(pdf,'clusterA_group2_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group3$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group3$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1)
dev.copy(pdf,'clusterA_group3_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group4$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group4$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1)
dev.copy(pdf,'clusterA_group4_rep_nterm_vs_cterm.pdf')
dev.off()





#Repressor similarity
par(mar=c(4,8,4,4))
plot(clusterA_group1_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group1_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group1$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group1$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1)
dev.copy(pdf,'clusterA_group1_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group2_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group2_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group2$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group2$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
abline(0,1)
dev.copy(pdf,'clusterA_group2_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group3_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group3_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group3$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group3$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,1)
dev.copy(pdf,'clusterA_group3_rep_nterm_vs_cterm.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(clusterA_group4_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group4_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group4$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group4$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,70),ylim=c(0,70),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,1)
dev.copy(pdf,'clusterA_group4_rep_nterm_vs_cterm.pdf')
dev.off()




#Zoom in
par(mar=c(4,8,4,4))
plot(clusterA_group1_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group1_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group1$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group1$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(clusterA_group2_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group2_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group2$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group2$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(clusterA_group3_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group3_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group3$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group3$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(clusterA_group4_diff$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group4_diff$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group4$repressor_nterm_mafft_dist_uncorrected,
     clusterA_group4$repressor_cterm_mafft_dist_uncorrected,
     xlim=c(0,40),ylim=c(0,40),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,1)
dev.copy(pdf,'clusterA_all_groups_zoom_in.pdf')
dev.off()

























par(mar=c(4,8,4,4))
plot(clusterA_group1_diff$repressor_nterm_mafft_phyml_dist,
     clusterA_group1_diff$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group1$repressor_nterm_mafft_phyml_dist,
     clusterA_group1$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1)

par(mar=c(4,8,4,4))
plot(clusterA_group2_diff$repressor_nterm_mafft_phyml_dist,
     clusterA_group2_diff$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group2$repressor_nterm_mafft_phyml_dist,
     clusterA_group2$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
abline(0,1)

par(mar=c(4,8,4,4))
plot(clusterA_group3_diff$repressor_nterm_mafft_phyml_dist,
     clusterA_group3_diff$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group3$repressor_nterm_mafft_phyml_dist,
     clusterA_group3$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,1)

par(mar=c(4,8,4,4))
plot(clusterA_group4_diff$repressor_nterm_mafft_phyml_dist,
     clusterA_group4_diff$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(clusterA_group4$repressor_nterm_mafft_phyml_dist,
     clusterA_group4$repressor_cterm_mafft_phyml_dist,
     xlim=c(0,4),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,1)







scatter.smooth(clusterA_group4$repressor_hth_compare,
               clusterA_group4$stoperator_pwd_dist_euc,
               pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")










#Immunity correlation

par(mar=c(4,8,4,4))
plot(distance_metrics$pham_pham_dissimilarity,
     distance_metrics$defending_cor_reduced,
     xlim=c(0,1),ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     distance_metrics$defending_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$stoperator_pwd_dist_euc,
     distance_metrics$defending_cor_reduced,
     ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(4,8,4,4))
plot(distance_metrics$stoperator_pwd_dist_euc,
     distance_metrics$challenging_cor_reduced,
     ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(distance_metrics$repressor_cterm_mafft_dist_uncorrected,
     distance_metrics$challenging_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




par(mar=c(4,8,4,4))
plot(clusterA_clade2$stoperator_pwd_dist_euc,
     clusterA_clade2$challenging_cor_reduced,
     ylim=c(-1,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)



par(mar=c(4,8,4,4))
plot(clusterA_clade2$repressor_cterm_mafft_dist_uncorrected,
     clusterA_clade2$challenging_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=19,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(clusterA_clade2_diff$repressor_cterm_mafft_dist_uncorrected,
     clusterA_clade2_diff$challenging_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=19,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")



par(mar=c(4,8,4,4))
plot(clusterA_clade2$repressor_cterm_mafft_dist_uncorrected,
     clusterA_clade2$defending_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=19,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(clusterA_clade2_diff$repressor_cterm_mafft_dist_uncorrected,
     clusterA_clade2_diff$defending_cor_reduced,
     xlim=c(0,70),ylim=c(-1,1),pch=19,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")

















#FinalFig
par(mar=c(4,8,16,4))
plot(clusterA_clade2$stoperator_pwd_dist_euc,
     clusterA_clade2$challenging_cor_reduced,
     xlim=c(0,5),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(clusterA_clade2_diff$stoperator_pwd_dist_euc,
     clusterA_clade2_diff$challenging_cor_reduced,
     xlim=c(0,5),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'clusterA_clade2_challCor_vs_stopEuc_subtypes.pdf')
dev.off()

#Number of comparisons:
nrow(subset(clusterA_clade2,!is.na(clusterA_clade2$challenging_cor_reduced))) +
  nrow(subset(clusterA_clade2_diff,!is.na(clusterA_clade2_diff$challenging_cor_reduced)))
nrow(subset(clusterA_clade2,!is.na(clusterA_clade2$challenging_cor_reduced)))
nrow(subset(clusterA_clade2_diff,!is.na(clusterA_clade2_diff$challenging_cor_reduced)))

            


#FinalFig
par(mar=c(4,8,16,4))
plot(clusterA_clade2$stoperator_pwd_dist_euc,
     clusterA_clade2$defending_cor_reduced,
     xlim=c(0,5),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(clusterA_clade2_diff$stoperator_pwd_dist_euc,
     clusterA_clade2_diff$defending_cor_reduced,
     xlim=c(0,5),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'clusterA_clade2_defCor_vs_stopEuc_subtypes.pdf')
dev.off()


#Number of comparisons:
nrow(subset(clusterA_clade2,!is.na(clusterA_clade2$defending_cor_reduced))) +
  nrow(subset(clusterA_clade2_diff,!is.na(clusterA_clade2_diff$defending_cor_reduced)))
nrow(subset(clusterA_clade2,!is.na(clusterA_clade2$defending_cor_reduced)))
nrow(subset(clusterA_clade2_diff,!is.na(clusterA_clade2_diff$defending_cor_reduced)))






#Final correlations
lm_immunity_def_cor_reduced <- lm(defending_cor_reduced ~
                                    stoperator_pwd_dist_euc,
                              data = clusterA_clade2)
summary(lm_immunity_def_cor_reduced)

lm_immunity_chal_cor_reduced <- lm(challenging_cor_reduced ~
                                    stoperator_pwd_dist_euc,
                                  data = clusterA_clade2)
summary(lm_immunity_chal_cor_reduced)


















par(mar=c(4,8,8,4))
plot(clusterA_clade2$defending_cor_reduced,
     clusterA_clade2$challenging_cor_reduced,
     xlim=c(-1,1),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(clusterA_clade2_diff$defending_cor_reduced,
     clusterA_clade2_diff$challenging_cor_reduced,
     xlim=c(-1,1),ylim=c(-1,1),pch=16,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
dev.copy(pdf,'clusterA_clade2_defCor_vs_chalCor_subtypes.pdf')
dev.off()



lm_immunity_cor_reduced <- lm(defending_cor_reduced ~
                                challenging_cor_reduced,
                              data = clusterA_clade2)
summary(lm_immunity_cor_reduced)



#Repressor size


#By subcluster
clusterA_clade2_myco_env_temp_rep <- subset(phage_metadata,
                                            phage_metadata$cluster == 'A' &
                                              phage_metadata$source == 'environment' &
                                              phage_metadata$host == 'Mycobacterium' &
                                              phage_metadata$cluster_a_functional_repressor_predicted == 'yes' &
                                              phage_metadata$gene_content_clade == 'clade2')


#FinalFig
par(mar=c(4,8,16,20))
boxplot(clusterA_clade2_myco_env_temp_rep$repressor_length_full,las=1,cex.axis=2,ann=FALSE,main=NULL,outline=FALSE,ylim=c(150,250),col="light grey")
par(new=TRUE)
stripchart(clusterA_clade2_myco_env_temp_rep$repressor_length_full,
           vertical=TRUE,las=1,cex.axis=2,pch=16,method="jitter",cex=1,ann=FALSE,main=NULL,ylim=c(150,250))
dev.copy(pdf,'clusterA_clade2_myco_env_temp_repFull_sizes.pdf')
dev.off()



### Whole genome metrics (regardless of immunity data) above
































###Cluster N data
setwd("~/scratch/immunity_analysis/input/")
cluster_n_immunity_data <- read.csv("cluster_n_immunity_data.csv",sep=",",header=TRUE)
names(cluster_n_immunity_data) <- c('defending_phage','challenging_phage','log10_eop')
cluster_n_immunity_data$log10_eop <- as.factor(cluster_n_immunity_data$log10_eop)


cluster_n_immunity_data$defending_challenging <- paste(cluster_n_immunity_data$defending_phage,
                                                       "_",
                                                       cluster_n_immunity_data$challenging_phage,
                                                       sep="")
cluster_n_immunity_data$defending_challenging <- as.factor(cluster_n_immunity_data$defending_challenging)

cluster_n_immunity_data <- merge(cluster_n_immunity_data,
                                 genomic_distance_data,
                                 by.x="defending_challenging",
                                 by.y="phage1_phage2")




#merge with metadata
cluster_n_phage_metadata <- read.csv("cluster_n_phage_metadata.csv",sep=",",header=TRUE)

cluster_n_phage_metadata_to_match <- cluster_n_phage_metadata
names(cluster_n_phage_metadata_to_match) <- paste('defending_',
                                                  names(cluster_n_phage_metadata_to_match),
                                                  sep="")
cluster_n_immunity_data <- merge(cluster_n_immunity_data,
                                 cluster_n_phage_metadata_to_match,
                                 by.x="defending_phage",
                                 by.y="defending_phage")

cluster_n_phage_metadata_to_match <- cluster_n_phage_metadata
names(cluster_n_phage_metadata_to_match) <- paste('challenging_',
                                                  names(cluster_n_phage_metadata_to_match),
                                                  sep="")
cluster_n_immunity_data <- merge(cluster_n_immunity_data,
                                 cluster_n_phage_metadata_to_match,
                                 by.x="challenging_phage",
                                 by.y="challenging_phage")

cluster_n_immunity_data$temperate_compare <- ifelse(cluster_n_immunity_data$defending_temperate==cluster_n_immunity_data$challenging_temperate,
                                                    as.character(cluster_n_immunity_data$defending_temperate),
                                                    "different")

cluster_n_immunity_data$temperate_compare <- factor(cluster_n_immunity_data$temperate_compare)




#Plot results
setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,4,4))
plot(cluster_n_immunity_data$modified_mash_distance,
     as.numeric(as.character(cluster_n_immunity_data$log10_eop)),
     xlim=c(0,0.5),ylim=c(-9,0),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'cluster_n_immunity_mash_vs_log10eop.pdf')
dev.off()

plot(cluster_n_immunity_data$pham_pham_dissimilarity,
     as.numeric(as.character(cluster_n_immunity_data$log10_eop)),
     pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'cluster_n_immunity_gcd_vs_log10eop.pdf')
dev.off()




cluster_n_immunity_data_temperate <- subset(cluster_n_immunity_data,
                                            cluster_n_immunity_data$temperate_compare == 'yes')

cluster_n_immunity_data_lytic <- subset(cluster_n_immunity_data,
                                        cluster_n_immunity_data$temperate_compare != 'yes')


par(mar=c(4,8,4,4))
plot(cluster_n_immunity_data_temperate$modified_mash_distance,
     as.numeric(as.character(cluster_n_immunity_data_temperate$log10_eop)),
     xlim=c(0,0.5),ylim=c(-9,0),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(cluster_n_immunity_data_lytic$modified_mash_distance,
     as.numeric(as.character(cluster_n_immunity_data_lytic$log10_eop)),
     xlim=c(0,0.5),ylim=c(-9,0),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')
dev.copy(pdf,'.pdf')
dev.off()


par(mar=c(4,8,4,4))
plot(cluster_n_immunity_data_temperate$pham_pham_dissimilarity,
     as.numeric(as.character(cluster_n_immunity_data_temperate$log10_eop)),
     xlim=c(0,1),ylim=c(-9,0),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='black')
par(new=TRUE)
plot(cluster_n_immunity_data_lytic$pham_pham_dissimilarity,
     as.numeric(as.character(cluster_n_immunity_data_lytic$log10_eop)),
     xlim=c(0,1),ylim=c(-9,0),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,col='green')



###End of Cluster N analysis


















###Coliphage immunity data



setwd("~/scratch/immunity_analysis/input/")


#Import immunity data
#1 = "assay_number"          
#2 = "defending_challenging"
#3 = "defending_phage"
#4 = "challenging_phage"
#5 = "infection"
#6 = "source"
coliphage_immunity_data <- read.csv("coliphage_immunity_data.csv",sep=",",header=TRUE)


#Import phage metadata
#Format:
#0 = phageid
#1 = description
#2 = accession
#3 = family
#4 = size
#5 = genes
#6 = lambdoid_immunity_group
#7 = p2_immunity_group
#8 = immunity_group_source
#9 = immunity_group
coliphage_metadata_table <- read.csv("coliphage_metadata.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
coliphage_metadata_table[coliphage_metadata_table == "Unspecified"] <- NA
coliphage_metadata_table[coliphage_metadata_table == "unspecified"] <- NA


#Modify column names of host data and merge with immunity table
#Data for all phages in processed mash table should be present in phage metadata table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Match metadata for both the query and reference phages in each pairwise comparison
coliphage_metadata_to_match <- coliphage_metadata_table
names(coliphage_metadata_to_match) <- paste('defending_',names(coliphage_metadata_to_match),sep="")
coliphage_immunity_data <- merge(coliphage_immunity_data,coliphage_metadata_to_match,by.x="defending_phage",by.y="defending_phageid")
coliphage_metadata_to_match <- coliphage_metadata_table
names(coliphage_metadata_to_match) <- paste('challenging_',names(coliphage_metadata_to_match),sep="")
coliphage_immunity_data <- merge(coliphage_immunity_data,coliphage_metadata_to_match,by.x="challenging_phage",by.y="challenging_phageid")






#Import mash dataset and match to gene content dissimilarity data
#Format:
#0 = reference genome
#1 = query genome
#2 = mash distance 
#3 = mash p-value
#4 = mash kmer count
#5 = reference_query
coliphage_mash_table <- read.csv("coliphage_processed_mash_output.csv",sep=",",header=TRUE)
names(coliphage_mash_table) <- c("mash_reference",
                       "mash_query",
                       "mash_distance",
                       "mash_pvalue",
                       "mash_count",
                       "mash_ref_query")




#Import pham data and merge with mash table
#Format
#0 = phage1_name
#1 = phage1_number_unshared_phams
#2 = phage1_shared_proportion
#3 = phage2_name
#4 = phage2_number_unshared_phams
#5 = phage2_shared_proportion
#6 = number_shared_phams
#7 = average_shared_proportion
#8 = jaccard_similarity
#9 = shared_pham_distribution_mean
#10 = shared_pham_distribution_median
#11 = shared_pham_distribution_max
#12 = unshared_pham_distribution_mean
#13 = unshared_pham_distribution_median
#14 = unshared_pham_distribution_max
#15 = unshared_orpham_count
coliphage_gcd_table <- read.csv("coliphage_pairwise_pham_proportions.csv",sep=",",header=TRUE)

#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
coliphage_gcd_table$pham_phage1_phage2 <- paste(coliphage_gcd_table$phage1_name,"_",coliphage_gcd_table$phage2_name,sep="")
coliphage_gcd_table$pham_phage1_phage2 <- as.factor(coliphage_gcd_table$pham_phage1_phage2)


#Compute gene content dissimilarity
coliphage_gcd_table$pham_dissimilarity <- 1 - coliphage_gcd_table$average_shared_proportion

coliphage_gcd_table <- subset(coliphage_gcd_table,select=c("pham_phage1_phage2",
                                                           "phage1_name",
                                                           "phage2_name",
                                                           "pham_dissimilarity"))

names(coliphage_gcd_table) <- c("pham_phage1_phage2",
                                "pham_phage1",
                                "pham_phage2",
                                "pham_pham_dissimilarity")

coliphage_genomic_distance_data <- merge(coliphage_mash_table,coliphage_gcd_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")


#Assign filter status and change mash distance if data is not significant
#Alternatively, the max percent parameter can be omitted with minimal effects on final analysis.
coliphage_genomic_distance_data$filter <- ifelse(coliphage_genomic_distance_data$mash_pvalue < 1e-10,TRUE,FALSE)



#Determine the max filtered distance
# temp_filtered_table <- subset(coliphage_genomic_distance_data,
#                               coliphage_genomic_distance_data$filter == TRUE)
#summary(temp_filtered_table$mash_distance)
#max mash distance = 0.5094640
#temp_subset <- subset(temp_filtered_table,temp_filtered_table$mash_distance > 0.499999)
#There is only one pairwise comparison > 0.4999 and it does not involve any coliphages of interest

#At this point, the max mash distance of all filtered comparisons < 0.5. So set the distance of all comparisons that did not pass the filter = 0.5
coliphage_genomic_distance_data$modified_mash_distance <- ifelse(coliphage_genomic_distance_data$filter == TRUE,
                                                                 coliphage_genomic_distance_data$mash_distance,
                                                                 0.5)


#Currently, the pairwise data does not contain self-comparisons (i.e. L5 compared to L5), 
#and it does not contain duplicate reciprocal pairs (i.e. L5-Trixie, if Trixie-L5 is already present)
#Analyzing immunity data requires all duplicate reciprocal data and self-comparisons
#First remove only the columns of interest
coliphage_partial_genomic_distance_data <- subset(coliphage_genomic_distance_data,
                                                  select=c("mash_reference",
                                                           "mash_query",
                                                           "modified_mash_distance",
                                                           "pham_pham_dissimilarity"))

#Create reciprocal identifiers
coliphage_partial_genomic_distance_data$ref_query <- paste(coliphage_partial_genomic_distance_data$mash_reference,
                                                           "_",
                                                           coliphage_partial_genomic_distance_data$mash_query,
                                                           sep="")
coliphage_partial_genomic_distance_data$query_ref <- paste(coliphage_partial_genomic_distance_data$mash_query,
                                                           "_",
                                                           coliphage_partial_genomic_distance_data$mash_reference,
                                                           sep="")
coliphage_partial_genomic_distance_data$ref_query <- as.factor(coliphage_partial_genomic_distance_data$ref_query)
coliphage_partial_genomic_distance_data$query_ref <- as.factor(coliphage_partial_genomic_distance_data$query_ref)

coliphage_duplicate_data_for_immunity1 <- subset(coliphage_partial_genomic_distance_data,
                                                 select=c("ref_query",
                                                          "modified_mash_distance",
                                                          "pham_pham_dissimilarity"))
coliphage_duplicate_data_for_immunity2 <- subset(coliphage_partial_genomic_distance_data,
                                                 select=c("query_ref",
                                                          "modified_mash_distance",
                                                          "pham_pham_dissimilarity"))

#Rename the columns so the two tables match
names(coliphage_duplicate_data_for_immunity1) <- c("phage1_phage2","modified_mash_distance","pham_pham_dissimilarity")
names(coliphage_duplicate_data_for_immunity2) <- c("phage1_phage2","modified_mash_distance","pham_pham_dissimilarity")



#Create a table of self-comparisons
coliphage_self_comparison_for_immunity <- subset(coliphage_metadata_table,
                                                 select=c("phageid"))
coliphage_self_comparison_for_immunity$phage1_phage2 <- paste(coliphage_self_comparison_for_immunity$phageid,
                                                              "_",
                                                              coliphage_self_comparison_for_immunity$phageid,
                                                              sep="")
coliphage_self_comparison_for_immunity$modified_mash_distance <- 0
coliphage_self_comparison_for_immunity$pham_pham_dissimilarity <- 0
coliphage_self_comparison_for_immunity <- subset(coliphage_self_comparison_for_immunity,
                                                 select=c("phage1_phage2",
                                                          "modified_mash_distance",
                                                          "pham_pham_dissimilarity"))

#Now merge the three tables = duplicate/reciprocal data and self-comparisons
coliphage_all_genomic_distance_data_for_immunity <- rbind(coliphage_duplicate_data_for_immunity1,
                                                          coliphage_duplicate_data_for_immunity2,
                                                          coliphage_self_comparison_for_immunity)
coliphage_all_genomic_distance_data_for_immunity$phage1_phage2 <- as.factor(coliphage_all_genomic_distance_data_for_immunity$phage1_phage2)

#merge immunity data with genomic distance data
coliphage_immunity_data <- merge(coliphage_immunity_data,
                                 coliphage_all_genomic_distance_data_for_immunity,
                                 by.x="defending_challenging",
                                 by.y="phage1_phage2")



coliphage_immunity_data$immunity_group_compare <- ifelse(coliphage_immunity_data$defending_immunity_group==
                                                           coliphage_immunity_data$challenging_immunity_group,
                                                         as.character(coliphage_immunity_data$defending_immunity_group),
                                                         "different")

coliphage_immunity_data$lambdoid_immunity_group_compare <- ifelse(coliphage_immunity_data$defending_lambdoid_immunity_group==
                                                                    coliphage_immunity_data$challenging_lambdoid_immunity_group,
                                                                  as.character(coliphage_immunity_data$defending_lambdoid_immunity_group),
                                                                  "different")

coliphage_immunity_data$p2_immunity_group_compare <- ifelse(coliphage_immunity_data$defending_p2_immunity_group==
                                                           coliphage_immunity_data$challenging_p2_immunity_group,
                                                         as.character(coliphage_immunity_data$defending_p2_immunity_group),
                                                         "different")




coliphage_immunity_data$immunity_group_compare <- factor(coliphage_immunity_data$immunity_group_compare)
coliphage_immunity_data$lambdoid_immunity_group_compare <- factor(coliphage_immunity_data$lambdoid_immunity_group_compare)
coliphage_immunity_data$p2_immunity_group_compare <- factor(coliphage_immunity_data$p2_immunity_group_compare)



#convert infection category to infection score
infection_score_table <- data.frame(c("yes","no","moderate"),c(1,0.5,0))
names(infection_score_table) <- c("infection_category","infection_score")

coliphage_immunity_data <- merge(coliphage_immunity_data,
                                 infection_score_table,
                                 by.x="infection",
                                 by.y="infection_category")


lambdoid_data <- subset(coliphage_immunity_data,coliphage_immunity_data$immunity_group_compare == 'lambdoid')
p2_data <- subset(coliphage_immunity_data,coliphage_immunity_data$immunity_group_compare == 'p2')

#Plot results
setwd("~/scratch/immunity_analysis/output/")

par(mar=c(4,8,4,4))
plot(lambdoid_data$modified_mash_distance,
     lambdoid_data$infection_score,
     xlim=c(0,0.5),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(p2_data$modified_mash_distance,
     p2_data$infection_score,
     xlim=c(0,0.5),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'.pdf')
dev.off()




par(mar=c(4,8,4,4))
plot(lambdoid_data$pham_pham_dissimilarity,
     lambdoid_data$infection_score,
     xlim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'.pdf')
dev.off()

par(mar=c(4,8,4,4))
plot(p2_data$pham_pham_dissimilarity,
     p2_data$infection_score,
     xlim=c(0,1),pch=1,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'.pdf')
dev.off()


###End of coliphage immunity analysis




























































###End of main immunity assay code

###################################################################################














###Old non-averaged misc analyses below

#Subset only the lysogen data
lysogen_immunity <- subset(main_immunity_data,main_immunity_data$strain_type == 'lysogen')
lysogen_immunity$defending_challenging <- as.factor(lysogen_immunity$defending_challenging)

lysogen_immunity$defending_phage_parent <- factor(lysogen_immunity$defending_phage_parent)
lysogen_immunity$challenging_phage_parent <- factor(lysogen_immunity$challenging_phage_parent)
lysogen_immunity$lawn_reliability <- factor(lysogen_immunity$lawn_reliability)
lysogen_immunity$phage_reliability <- factor(lysogen_immunity$phage_reliability)


#Match the genomic distance data. The double lysogen strain carrying RedRock and Bxb1 should have a lookup value of NA,
#and the challenging phage phiTM42 (RedRock-Trixie hybrid) should have lookup value of NA
lysogen_immunity <- merge(lysogen_immunity,genomic_distance_data,by.x="defending_challenging",by.y="ref_query")



lysogen_immunity_reliable <- subset(lysogen_immunity,lysogen_immunity$lawn_reliability != 1)
lysogen_immunity_reliable <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$phage_reliability == 2 |
                                      lysogen_immunity_reliable$phage_reliability == 3)


lysogen_immunity_reliable <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$defending_phage_source == 'seaphages' &
                                      lysogen_immunity_reliable$challenging_phage_source == 'seaphages')


plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_four_factors)
hist(lysogen_immunity_reliable$modified_mash_distance,breaks=100,col="black")

plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_infection_strength)

plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_turbidity)

plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_plaque_size)

plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_plaques)

integrating_immunity_reliable <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$defending_phage_lysogen_type == 'integration')
extrachromosomal_immunity_reliable <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$defending_phage_lysogen_type == 'extrachromosomal')

both_integrating_immunity_reliable <- subset(integrating_immunity_reliable,integrating_immunity_reliable$challenging_phage_lysogen_type == 'integration')
both_extrachromosomal_immunity_reliable <- subset(extrachromosomal_immunity_reliable,extrachromosomal_immunity_reliable$challenging_phage_lysogen_type == 'extrachromosomal')


plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_four_factors,xlim=c(0,0.4))


plot(integrating_immunity_reliable$modified_mash_distance,integrating_immunity_reliable$score_four_factors,xlim=c(0,0.4))
plot(both_integrating_immunity_reliable$modified_mash_distance,both_integrating_immunity_reliable$score_four_factors,xlim=c(0,0.4))


plot(extrachromosomal_immunity_reliable$modified_mash_distance,extrachromosomal_immunity_reliable$score_four_factors,xlim=c(0,0.4))
plot(both_extrachromosomal_immunity_reliable$modified_mash_distance,both_extrachromosomal_immunity_reliable$score_four_factors,xlim=c(0,0.4))






#Evaluate single titer repressor clone data

single_titer_all <- subset(main_immunity_data,main_immunity_data$assay_type == 'single_titer')
single_titer_all$defending_phage <- factor(single_titer_all$defending_phage)
single_titer_all$challenging_phage <- factor(single_titer_all$challenging_phage)
single_titer_all$defending_challenging <- factor(single_titer_all$defending_challenging)
single_titer_all$defending_challenging_parent <- factor(single_titer_all$defending_challenging_parent)
single_titer_all$lawn_reliability <- factor(single_titer_all$lawn_reliability)
single_titer_all$phage_reliability <- factor(single_titer_all$phage_reliability)





#Only keep reliable data
single_titer_reliable <- subset(single_titer_all,single_titer_all$lawn_reliability != 1)
single_titer_reliable <- subset(single_titer_reliable,single_titer_reliable$phage_reliability == 2 |
                                  single_titer_reliable$phage_reliability == 3)


single_titer_lysogen <- subset(single_titer_reliable,single_titer_reliable$strain_type == 'lysogen')
single_titer_clone <- subset(single_titer_reliable,single_titer_reliable$strain_type == 'repressor_clone')



single_titer_lysogen_subset <- subset(single_titer_lysogen, select = c("defending_challenging",
                                                                "challenging_phage",
                                                                "challenging_phage_source",
                                                                "defending_phage",
                                                                "defending_phage_source",
                                                                "score_infection_strength",
                                                                "score_turbidity",
                                                                "score_plaque_size",
                                                                "score_plaques",
                                                                "score_four_factors"))


single_titer_clone_subset <- subset(single_titer_clone, select = c("defending_challenging",
                                                                       "challenging_phage",
                                                                       "challenging_phage_source",
                                                                       "defending_phage",
                                                                       "defending_phage_source",
                                                                       "score_infection_strength",
                                                                       "score_turbidity",
                                                                       "score_plaque_size",
                                                                       "score_plaques",
                                                                       "score_four_factors"))


names(single_titer_lysogen_subset) <- c("defending_challenging",
                                 "lysogen_challenging_phage",
                                 "lysogen_challenging_phage_source",
                                 "lysogen_defending_phage",
                                 "lysogen_defending_phage_source",
                                 "lysogen_score_infection_strength",
                                 "lysogen_score_turbidity",
                                 "lysogen_score_plaque_size",
                                 "lysogen_score_plaques",
                                 "lysogen_score_four_factors")


names(single_titer_clone_subset) <- c("defending_challenging",
                                        "clone_challenging_phage",
                                        "clone_challenging_phage_source",
                                        "clone_defending_phage",
                                        "clone_defending_phage_source",
                                        "clone_score_infection_strength",
                                        "clone_score_turbidity",
                                        "clone_score_plaque_size",
                                        "clone_score_plaques",
                                        "clone_score_four_factors")



single_titer_matched <- merge(single_titer_lysogen_subset,single_titer_clone_subset,by.x = "defending_challenging",by.y = "defending_challenging")

single_titer_matched$score_four_factors_diff <- single_titer_matched$lysogen_score_four_factors - single_titer_matched$clone_score_four_factors
single_titer_matched$score_infection_strength_diff <- single_titer_matched$lysogen_score_infection_strength - single_titer_matched$clone_score_infection_strength
single_titer_matched$score_turbidity_diff <- single_titer_matched$lysogen_score_turbidity - single_titer_matched$clone_score_turbidity
              


single_titer_matched_wildtype <- subset(single_titer_matched,single_titer_matched$lysogen_challenging_phage_source == 'seaphages' &
                                          single_titer_matched$clone_challenging_phage_source == 'seaphages')

plot(single_titer_matched$lysogen_score_four_factors,single_titer_matched$clone_score_four_factors,xlim=c(0,5),ylim=c(0,5))
plot(single_titer_matched$lysogen_score_infection_strength,single_titer_matched$clone_score_infection_strength,xlim=c(0,3),ylim=c(0,3))
plot(single_titer_matched$lysogen_score_turbidity,single_titer_matched$clone_score_turbidity,xlim=c(0,2),ylim=c(0,2))
plot(single_titer_matched$score_four_factors_diff)



plot(single_titer_matched_wildtype$lysogen_score_four_factors,single_titer_matched_wildtype$clone_score_four_factors,xlim=c(0,5),ylim=c(0,5))
plot(single_titer_matched_wildtype$lysogen_score_infection_strength,single_titer_matched_wildtype$clone_score_infection_strength,xlim=c(0,3),ylim=c(0,3))
plot(single_titer_matched_wildtype$lysogen_score_turbidity,single_titer_matched_wildtype$clone_score_turbidity,xlim=c(0,2),ylim=c(0,2))
plot(single_titer_matched_wildtype$score_four_factors_diff)



#plot frequencies of each infection score to compare wildtype phage to mutant phage
barplot(summary(as.factor(single_titer_matched$score_four_factors_diff)),ylim=c(0,130))
barplot(summary(as.factor(single_titer_matched_wildtype$score_four_factors_diff)),ylim=c(0,130))

#plot(single_titer_matched$lysogen_score_plaque_size,single_titer_matched$clone_score_plaque_size,xlim=c(0,5),ylim=c(0,5))
#plot(single_titer_matched$lysogen_score_plaques,single_titer_matched$clone_score_plaques,xlim=c(0,5),ylim=c(0,5))
                   

"immunity_assay_id"                             
"immunity._set"                                 
"date"                                          
"notebook"                                      
"page"                                          
"strain"                                        
"prophage"                                      
"repressor_clone"                               
"strain_type"                                   
"assay_type"                                    
"lawn_notes"                                    
"lawn_reliability"                              
"tested_titer"                                  
"phage_reliability"                             
"observation_infection_strength"                
"observation_turbidity"                         
"observation_plaque_size"                       
"observation_plaques"                           
"defending_phage_source"                        
"defending_phage_alias"                         
"defending_phage_description"                   
"challenging_phage_source"                      
"challenging_phage_alias"                       
"challenging_phage_description"                 
"defending_challenging_parent"                  
"defending_challenging"                         
"defending_phage_host"                          
"defending_phage_cluster"                       
"defending_phage_subcluster"                    
"defending_phage_size"                          
"defending_phage_status"                        
"defending_phage_author"                        
"defending_phage_datelastmodified"              
"defending_phage_mode_approx_80_percent"        
"defending_phage_network005_interaction_tally"  
"defending_phage_network005_mash_group"         
"defending_phage_network005_mash_group_tally"   
"defending_phage_lysogen_type"                  
"defending_phage_pham_integrase"                
"defending_phage_pham_para"                     
"challenging_phage_host"                        
"challenging_phage_cluster"                     
"challenging_phage_subcluster"                  
"challenging_phage_size"                        
"challenging_phage_status"                      
"challenging_phage_author"                      
"challenging_phage_datelastmodified"            
"challenging_phage_mode_approx_80_percent"      
"challenging_phage_network005_interaction_tally"
"challenging_phage_network005_mash_group"       
"challenging_phage_network005_mash_group_tally" 
"challenging_phage_lysogen_type"                
"challenging_phage_pham_integrase"              
"challenging_phage_pham_para"                   
"challenging_phage_parent"                      
"defending_phage_parent"                        

















#Now compare prophage and repressor clone multi-titer assay data


multi_lysogen_clone_all <- subset(main_immunity_data,main_immunity_data$assay_type == 'multiple_titer' & main_immunity_data$date == '6/17/18')
multi_lysogen_clone_all$defending_phage <- factor(multi_lysogen_clone_all$defending_phage)
multi_lysogen_clone_all$challenging_phage <- factor(multi_lysogen_clone_all$challenging_phage)
multi_lysogen_clone_all$defending_challenging <- factor(multi_lysogen_clone_all$defending_challenging)
multi_lysogen_clone_all$defending_challenging_parent <- factor(multi_lysogen_clone_all$defending_challenging_parent)
multi_lysogen_clone_all$lawn_reliability <- factor(multi_lysogen_clone_all$lawn_reliability)
multi_lysogen_clone_all$phage_reliability <- factor(multi_lysogen_clone_all$phage_reliability)


#Only keep reliable data
multi_lysogen_clone_reliable <- subset(multi_lysogen_clone_all,multi_lysogen_clone_all$lawn_reliability != 1)
multi_lysogen_clone_reliable <- subset(multi_lysogen_clone_reliable,multi_lysogen_clone_reliable$phage_reliability == 2 |
                                         multi_lysogen_clone_reliable$phage_reliability == 3)


multi_lysogen <- subset(multi_lysogen_clone_reliable,multi_lysogen_clone_reliable$strain_type == 'lysogen')
multi_clone <- subset(multi_lysogen_clone_reliable,multi_lysogen_clone_reliable$strain_type == 'repressor_clone')



multi_lysogen_subset <- subset(multi_lysogen, select = c("defending_challenging",
                                                                       "challenging_phage",
                                                                       "challenging_phage_source",
                                                                       "defending_phage",
                                                                       "defending_phage_source",
                                                                       "score_infection_strength",
                                                                       "score_turbidity",
                                                                       "score_plaque_size",
                                                                       "score_plaques",
                                                                       "score_four_factors"))


multi_clone_subset <- subset(multi_clone, select = c("defending_challenging",
                                                                   "challenging_phage",
                                                                   "challenging_phage_source",
                                                                   "defending_phage",
                                                                   "defending_phage_source",
                                                                   "score_infection_strength",
                                                                   "score_turbidity",
                                                                   "score_plaque_size",
                                                                   "score_plaques",
                                                                   "score_four_factors"))


names(multi_lysogen_subset) <- c("defending_challenging",
                                        "lysogen_challenging_phage",
                                        "lysogen_challenging_phage_source",
                                        "lysogen_defending_phage",
                                        "lysogen_defending_phage_source",
                                        "lysogen_score_infection_strength",
                                        "lysogen_score_turbidity",
                                        "lysogen_score_plaque_size",
                                        "lysogen_score_plaques",
                                        "lysogen_score_four_factors")


names(multi_clone_subset) <- c("defending_challenging",
                                      "clone_challenging_phage",
                                      "clone_challenging_phage_source",
                                      "clone_defending_phage",
                                      "clone_defending_phage_source",
                                      "clone_score_infection_strength",
                                      "clone_score_turbidity",
                                      "clone_score_plaque_size",
                                      "clone_score_plaques",
                                      "clone_score_four_factors")



multi_lysogen_clone_matched <- merge(multi_lysogen_subset,multi_clone_subset,by.x = "defending_challenging",by.y = "defending_challenging")


multi_lysogen_clone_matched$score_four_factors_diff <- multi_lysogen_clone_matched$lysogen_score_four_factors - multi_lysogen_clone_matched$clone_score_four_factors
multi_lysogen_clone_matched$score_infection_strength_diff <- multi_lysogen_clone_matched$lysogen_score_infection_strength - multi_lysogen_clone_matched$clone_score_infection_strength
multi_lysogen_clone_matched$score_turbidity_diff <- multi_lysogen_clone_matched$lysogen_score_turbidity - multi_lysogen_clone_matched$clone_score_turbidity



multi_lysogen_clone_matched_wildtype <- subset(multi_lysogen_clone_matched,multi_lysogen_clone_matched$lysogen_challenging_phage_source == 'seaphages' &
                                          multi_lysogen_clone_matched$clone_challenging_phage_source == 'seaphages')

plot(multi_lysogen_clone_matched$lysogen_score_four_factors,multi_lysogen_clone_matched$clone_score_four_factors)
plot(multi_lysogen_clone_matched$lysogen_score_infection_strength,multi_lysogen_clone_matched$clone_score_infection_strength)
plot(multi_lysogen_clone_matched$lysogen_score_turbidity,multi_lysogen_clone_matched$clone_score_turbidity)
plot(multi_lysogen_clone_matched$score_four_factors_diff)



plot(multi_lysogen_clone_matched_wildtype$lysogen_score_four_factors,multi_lysogen_clone_matched_wildtype$clone_score_four_factors,xlim=c(0,5),ylim=c(0,5))
plot(multi_lysogen_clone_matched_wildtype$lysogen_score_infection_strength,multi_lysogen_clone_matched_wildtype$clone_score_infection_strength,xlim=c(0,3),ylim=c(0,3))
plot(multi_lysogen_clone_matched_wildtype$lysogen_score_turbidity,multi_lysogen_clone_matched_wildtype$clone_score_turbidity,xlim=c(0,2),ylim=c(0,2))
plot(multi_lysogen_clone_matched_wildtype$score_four_factors_diff)



#plot frequencies of each infection score to compare wildtype phage to mutant phage
barplot(summary(as.factor(multi_lysogen_clone_matched$score_four_factors_diff)),ylim=c(0,25))
barplot(summary(as.factor(multi_lysogen_clone_matched_wildtype$score_four_factors_diff)),ylim=c(0,25))









#Plot immunity profiles per phage
#assumes lawn and phage reliability and lysogen matched to genomic distance


table <- subset(main_immunity_data,main_immunity_data$defending_phage == 'trixie')
table <- subset(table,table$strain_type == 'lysogen')
table <- subset(table,table$lawn_reliability != 1)
table <- subset(table,table$phage_reliability == 2 |
                  table$phage_reliability == 3)

plot(lysogen_immunity_reliable$modified_mash_distance,lysogen_immunity_reliable$score_four_factors)

plot(table$)

plot_immunity_profiles_by_phage <- function(table){

    
  table <- subset(table,table$)
  
  
  
  
}
  
  
  
  








#Plot superinfection profiles, matching escape mutant with parent phage


immunity_reliable <- subset(immunity_data,immunity_data$lawn_reliability != 1)
immunity_reliable <- subset(immunity_reliable,immunity_reliable$phage_reliability == 2 |
                              immunity_reliable$phage_reliability == 3)


# phitm46_data <- subset(immunity_reliable,immunity_reliable$challenging_phage == 'phitm46' |
#                          immunity_reliable$challenging_phage == 'davinci')



phitm46_data <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$challenging_phage == 'phitm46')
davinci_data <- subset(lysogen_immunity_reliable,lysogen_immunity_reliable$challenging_phage == 'davinci')


plot(phitm46_data$modified_mash_distance,phitm46_data$score_four_factors)



###








































###Extra

lys_conf_envDY <- subset(lys_conf,
                         lys_conf$defending_source == 'environment')


lys_conf_envDYCY <- subset(lys_conf_envDY,
                           lys_conf_envDY$challenging_source == 'environment')

lys_conf_envDYCN <- subset(lys_conf_envDY,
                           lys_conf_envDY$challenging_source == 'lab')

lys_conf_envDYCY_intDY <- subset(lys_conf_envDYCY,
                                 lys_conf_envDYCY$defending_lysogen_type == 'integration')


lys_conf_envDY_intDN <- subset(lys_conf_envDY,
                               lys_conf_envDY$defending_lysogen_type == 'extrachromosomal')



lys_conf_envDY_intDY_intCY <- subset(lys_conf_envDY_intDY,
                                     lys_conf_envDY_intDY$challenging_lysogen_type == 'integration')


lys_conf_envDY_intDN_intCN <- subset(lys_conf_envDY,
                                     lys_conf_envDY$challenging_lysogen_type == 'extrachromosomal')




lysogen_immunity_genomic_distances_reliable_env_repYes <- subset(lysY_mash_conf_envY,
                                                                 lysY_mash_conf_envY$challenging_cluster_a_functional_repressor_predicted == 'yes')


lysogen_immunity_genomic_distances_reliable_env_repNo <- subset(lysogen_immunity_genomic_distances_reliable_env,
                                                                lysogen_immunity_genomic_distances_reliable_env$challenging_cluster_a_functional_repressor_predicted == 'no')

lysogen_immunity_genomic_distances_reliable_env_repYes_intYes <- subset(lysogen_immunity_genomic_distances_reliable_env_repYes,
                                                                        lysogen_immunity_genomic_distances_reliable_env_repYes$defending_lysogen_type == 'integration')

lysogen_immunity_genomic_distances_reliable_env_repYes_extraYes <- subset(lysogen_immunity_genomic_distances_reliable_env_repYes_intYes,
                                                                          lysogen_immunity_genomic_distances_reliable_env_repYes_intYes$challenging_lysogen_type == 'extrachromosomal')


lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo <- subset(lysogen_immunity_genomic_distances_reliable_env_repYes,
                                                                         lysogen_immunity_genomic_distances_reliable_env_repYes$challenging_lysogen_type == 'integration')

lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo_phamSame <- subset(lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo,
                                                                                  lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo$integrase_compare != 'different')

lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo_phamDiff <- subset(lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo,
                                                                                  lysogen_immunity_genomic_distances_reliable_env_repYes_extraNo$integrase_compare == 'different')














































###Extra

###Assess scored_four_factors frequencies to compute new scoring total

#Multi, conf, lysogen data
setwd("~/scratch/immunity_analysis/output/")
multi_conf_lys_unique_scores <- subset(main_immunity_data,
                                       main_immunity_data$assay_type == 'multiple_titer' &
                                         main_immunity_data$lawn_reliability != 1 &
                                         main_immunity_data$phage_reliability != 1 &
                                         main_immunity_data$strain_type == 'lysogen',
                                       select=c("scored_infection_strength",
                                                "scored_turbidity",
                                                "scored_plaque_size",
                                                "scored_plaques"))


multi_conf_lys_unique_scores$scores_grouped <- paste(multi_conf_lys_unique_scores$scored_infection_strength,"_",
                                                     multi_conf_lys_unique_scores$scored_turbidity,"_",
                                                     multi_conf_lys_unique_scores$scored_plaques,"_",
                                                     multi_conf_lys_unique_scores$scored_plaque_size,
                                                     sep="")

multi_conf_lys_unique_scores$scores_grouped <- factor(multi_conf_lys_unique_scores$scores_grouped)

# #I may need to change this to table function. Summary doesn't output all values (only first 100).But with table function, columns are different order, and no need for rownames field.
# multi_conf_lys_unique_scores_count <- as.data.frame(summary(multi_conf_lys_unique_scores$scores_grouped))
# multi_conf_lys_unique_scores_count$scores_grouped <- rownames(multi_conf_lys_unique_scores_count)
# names(multi_conf_lys_unique_scores_count) <- c('frequency','scores_grouped')
# #



multi_conf_lys_unique_scores <- multi_conf_lys_unique_scores[!duplicated(multi_conf_lys_unique_scores),]
multi_conf_lys_unique_scores <- merge(multi_conf_lys_unique_scores,multi_conf_lys_unique_scores_count,by.x='scores_grouped',by.y='scores_grouped')


write.table(multi_conf_lys_unique_scores,
            "multi_conf_lys_unique_scores.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


multi_conf_lys_unique_scores_ordered <- multi_conf_lys_unique_scores[order(multi_conf_lys_unique_scores$scored_infection_strength,
                                                                           multi_conf_lys_unique_scores$scored_turbidity,
                                                                           multi_conf_lys_unique_scores$scored_plaques,
                                                                           multi_conf_lys_unique_scores$scored_plaque_size),]
multi_conf_lys_unique_scores_ordered <- multi_conf_lys_unique_scores_ordered[,c(1,2,4,3)]
multi_conf_lys_unique_scores_ordered$scored_rank20 <- c(0:20)





#
#Multi, conf, repressor_clone and lysogen data
setwd("~/scratch/immunity_analysis/output/")
multi_conf_unique_scores <- subset(main_immunity_data,
                                   main_immunity_data$assay_type == 'multiple_titer' &
                                     main_immunity_data$lawn_reliability != 1 &
                                     main_immunity_data$phage_reliability != 1,
                                   select=c("scored_infection_strength",
                                            "scored_turbidity",
                                            "scored_plaque_size",
                                            "scored_plaques"))


multi_conf_unique_scores$scores_grouped <- paste(multi_conf_unique_scores$scored_infection_strength,"_",
                                                 multi_conf_unique_scores$scored_turbidity,"_",
                                                 multi_conf_unique_scores$scored_plaques,"_",
                                                 multi_conf_unique_scores$scored_plaque_size,
                                                 sep="")

multi_conf_unique_scores$scores_grouped <- factor(multi_conf_unique_scores$scores_grouped)

# #I may need to change this to table function. Summary doesn't output all values (only first 100).  But with table function, columns are different order, and no need for rownames field.
# multi_conf_unique_scores_count <- as.data.frame(summary(multi_conf_unique_scores$scores_grouped))
# multi_conf_unique_scores_count$scores_grouped <- rownames(multi_conf_unique_scores_count)
# names(multi_conf_unique_scores_count) <- c('frequency','scores_grouped')
# #

multi_conf_unique_scores <- multi_conf_unique_scores[!duplicated(multi_conf_unique_scores),]
multi_conf_unique_scores <- merge(multi_conf_unique_scores,multi_conf_unique_scores_count,by.x='scores_grouped',by.y='scores_grouped')


write.table(multi_conf_unique_scores,
            "multi_conf_unique_scores.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)


# multi_conf_unique_scores_ordered <- multi_conf_unique_scores[order(multi_conf_unique_scores$scored_infection_strength,
#                                                                            multi_conf_unique_scores$scored_turbidity,
#                                                                            multi_conf_unique_scores$scored_plaques,
#                                                                            multi_conf_unique_scores$scored_plaque_size),]
# multi_conf_unique_scores_ordered <- multi_conf_unique_scores_ordered[,c(1,2,4,3)]
# multi_conf_unique_scores_ordered$scored_rank20 <- c(0:20)








#
#Single, conf, repressor_clone and lysogen data
setwd("~/scratch/immunity_analysis/output/")
single_conf_unique_scores <- subset(main_immunity_data,
                                    main_immunity_data$assay_type == 'single_titer' &
                                      main_immunity_data$lawn_reliability != 1 &
                                      main_immunity_data$phage_reliability != 1,
                                    select=c("scored_infection_strength",
                                             "scored_turbidity",
                                             "scored_plaque_size",
                                             "scored_plaques"))


single_conf_unique_scores$scores_grouped <- paste(single_conf_unique_scores$scored_infection_strength,"_",
                                                  single_conf_unique_scores$scored_turbidity,"_",
                                                  single_conf_unique_scores$scored_plaques,"_",
                                                  single_conf_unique_scores$scored_plaque_size,
                                                  sep="")

single_conf_unique_scores$scores_grouped <- factor(single_conf_unique_scores$scores_grouped)

# #I may need to change this to table function. Summary doesn't output all values (only first 100). But with table function, columns are different order, and no need for rownames field.
# single_conf_unique_scores_count <- as.data.frame(summary(single_conf_unique_scores$scores_grouped))
# single_conf_unique_scores_count$scores_grouped <- rownames(single_conf_unique_scores_count)
# names(single_conf_unique_scores_count) <- c('frequency','scores_grouped')
# #




single_conf_unique_scores <- single_conf_unique_scores[!duplicated(single_conf_unique_scores),]
single_conf_unique_scores <- merge(single_conf_unique_scores,single_conf_unique_scores_count,by.x='scores_grouped',by.y='scores_grouped')


write.table(single_conf_unique_scores,
            "single_conf_unique_scores.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)




#
#all conf data
setwd("~/scratch/immunity_analysis/output/")

conf_unique_scores <- subset(main_immunity_data,
                             main_immunity_data$lawn_reliability != 1 &
                               main_immunity_data$phage_reliability != 1,
                             select=c("scored_infection_strength",
                                      "scored_turbidity",
                                      "scored_plaque_size",
                                      "scored_plaques"))


conf_unique_scores$scores_grouped <- paste(conf_unique_scores$scored_infection_strength,"_",
                                           conf_unique_scores$scored_turbidity,"_",
                                           conf_unique_scores$scored_plaques,"_",
                                           conf_unique_scores$scored_plaque_size,
                                           sep="")

conf_unique_scores$scores_grouped <- factor(conf_unique_scores$scores_grouped)


# #I may need to change this to table function. Summary doesn't output all values (only first 100). But with table function, columns are different order, and no need for rownames field.
# conf_unique_scores_count <- as.data.frame(summary(conf_unique_scores$scores_grouped))
# conf_unique_scores_count$scores_grouped <- rownames(conf_unique_scores_count)
# names(conf_unique_scores_count) <- c('frequency','scores_grouped')
# #




conf_unique_scores <- conf_unique_scores[!duplicated(conf_unique_scores),]
conf_unique_scores <- merge(conf_unique_scores,conf_unique_scores_count,by.x='scores_grouped',by.y='scores_grouped')


write.table(conf_unique_scores,
            "conf_unique_scores.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)



#



rank_table <- expand.grid(scored_infection_strength = c(0:3),
                          scored_turbidity = c(0:2),
                          scored_plaques = c(0,1),
                          scored_plaque_size = c(0:2))
rank_table <- rank_table[order(rank_table$scored_infection_strength,
                               rank_table$scored_turbidity,
                               rank_table$scored_plaques,
                               rank_table$scored_plaque_size),]


###







