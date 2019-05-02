# R script to process genomic distance data in preparation for immunity data
# analysis in the analyze_immunity_data.R script published in
# Mavrich & Hatfull, mBio, 2019. 
# Travis Mavrich.
# Note: this code merges and analyzes various input data sets 
# prepared by other tools including Python and Excel.
# Whole genome nucleotide distance was computed from Mash.
# Whole genome gene content dissimilarity was computed from Python scripts.
# Both data sets were computed as published in Mavrich & Hatfull, Nature
# Microbiology 2017.
# Distinct analyses and code blocks separated by "###".

### 1. Prepare environment.
### 2. Define functions.
### 3. Import mash and gene content dissimilarity data.
### 4. Output data for immunity analysis.
### 5. Output data for predicting evolutionary mode (if needed).
### 6. Output data for computing MaxGCDGap (if needed).
### 7. Output data for creating phage networks (if needed).



### 1. Prepare environment.

# Set working directory variables specific to local directory structure.
DIR_INPUT = "~/scratch/process_genome_distances/input/"
DIR_OUTPUT = "~/scratch/process_genome_distances/output/"

# Set paths for all input files.

PHAGE_METADATA_FILENAME = 
  paste(DIR_INPUT,
        "phage_metadata.csv",
        sep="")

MASH_DATA_FILENAME = 
  paste(DIR_INPUT,
        "processed_mash_output.csv",
        sep="")

PHAM_DATA_FILENAME = 
  paste(DIR_INPUT,
        "pairwise_pham_proportions.csv",
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

# The standard genomic similarity plot parameters.
plot_genomic_similarity_standard <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,
       table$pham_pham_dissimilarity,
       xlim=c(0,0.5),ylim=c(0,1),
       pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
  abline(0,2,lty=2,lwd=3,col="grey")
  
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
### 3. Import mash and gene content dissimilarity data.


# Import mash dataset.
# This data is generated from the process_mash_data.py script.
# Data structure:
# 1 = reference genome
# 2 = query genome
# 3 = mash distance 
# 4 = mash p-value
# 5 = mash kmer count
# 6 = reference_query
mash_table <- read.csv(MASH_DATA_FILENAME,sep=",",header=TRUE)
names(mash_table) <- c("mash_reference",
                       "mash_query",
                       "mash_distance",
                       "mash_pvalue",
                       "mash_count",
                       "mash_ref_query")



# Import phage metadata. Derived from same metadata input used for immunity
# data analysis. Contains data for all phages in Actino1321 database.
# This data is generated manually in Excel.
# Data structure:
# 1 = phageid
# 2 = host
# 3 = cluster
# 4 = subcluster
phage_metadata_table <- read.csv(PHAGE_METADATA_FILENAME,sep=",",header=TRUE)


# Convert all Unspecified fields to NA missing value
phage_metadata_table[phage_metadata_table == "Unspecified"] <- NA



# Match metadata for the query and reference phages in each pairwise comparison.
# Data for all phages in processed mash table should be present in the
# phage metadata table. As a result, do not select all.x=TRUE option. Any
# missing rows indicates an error in the table.
metadata_to_match <- phage_metadata_table

names(metadata_to_match) <- c("phageid",
                              "query_host_genus",
                              "query_cluster",
                              "query_subcluster")

main_data_table <- merge(mash_table,
                         metadata_to_match,
                         by.x="mash_query",
                         by.y="phageid")

names(metadata_to_match) <- c("phageid",
                              "ref_host_genus",
                              "ref_cluster",
                              "ref_subcluster")

main_data_table <- merge(main_data_table,
                         metadata_to_match,
                         by.x="mash_reference",
                         by.y="phageid")


# Import pairwise pham proportion data and merge with mash table.
# This data is generated from the analyze_pham_data.py script.
# Data structure:
# 1 = phage1
# 2 = phage1_number_unshared_phams
# 3 = phage1_shared_proportion
# 4 = phage2
# 5 = phage2_number_unshared_phams
# 6 = phage2_shared_proportion
# 7 = number_shared_phams
# 8 = average_shared_proportion
# 9 = jaccard_similarity
# 10 = shared_pham_distribution_mean
# 11 = shared_pham_distribution_median
# 12 = shared_pham_distribution_max
# 13 = unshared_pham_distribution_mean
# 14 = unshared_pham_distribution_median
# 15 = unshared_pham_distribution_max
# 16 = unshared_orpham_count
pham_table <- read.csv(PHAM_DATA_FILENAME,sep=",",header=TRUE)

# Compute gene content dissimilarity.
pham_table$pham_dissimilarity <- 1 - pham_table$average_shared_proportion
pham_table$jaccard_dissimilarity <- 1 - pham_table$jaccard_similarity

names(pham_table) <- c("pham_phage1",
                       "pham_phage1_number_unshared_phams",
                       "pham_phage1_shared_proportion",
                       "pham_phage2",
                       "pham_phage2_number_unshared_phams",
                       "pham_phage2_shared_proportion",
                       "pham_number_shared_phams",
                       "pham_average_shared_proportion",
                       "pham_jaccard_similarity",
                       "pham_shared_pham_distribution_mean",
                       "pham_shared_pham_distribution_median",
                       "pham_shared_pham_distribution_max",
                       "pham_unshared_pham_distribution_mean",
                       "pham_unshared_pham_distribution_median",
                       "pham_unshared_pham_distribution_max",
                       "pham_unshared_orpham_count",
                       "pham_pham_dissimilarity",
                       "pham_jaccard_dissimilarity")


# Since pham data contains pairwise duplicates, no need to worry about the
# phage order when creating ref_query match column.
pham_table$pham_phage1_phage2 <- paste(pham_table$pham_phage1,
                                       "_",
                                       pham_table$pham_phage2,
                                       sep="")

pham_table$pham_phage1_phage2 <- as.factor(pham_table$pham_phage1_phage2)

main_data_table <- merge(main_data_table,pham_table,
                         by.x="mash_ref_query",
                         by.y="pham_phage1_phage2")

# Assign filter status and change mash distance if data is not significant.
# main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10 & 
#                                    main_data_table$size_diff_max_percent < 1,
#                                  TRUE,
#                                  FALSE)

# Alternatively, the max percent parameter can be omitted with minimal effects
# on final analysis.
main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10,
                                 TRUE,
                                 FALSE)


#  QC:Determine the max filtered distance.
filtered_table <- subset(main_data_table,
                         main_data_table$filter == TRUE)

summary(filtered_table$mash_distance)
# The max mash distance = 0.4991900.

# At this point, the max mash distance of all filtered comparisons < 0.5. So
# set the distance of all comparisons that did not pass the filter = 0.5.
main_data_table$modified_mash_distance <- 
  ifelse(main_data_table$filter == TRUE,
         main_data_table$mash_distance,0.5)

# Compare ref and query host data columns.
main_data_table$host_genus_compare <- 
  ifelse(main_data_table$ref_host_genus ==
           main_data_table$query_host_genus,
         as.character(main_data_table$ref_host_genus),
         "different")

main_data_table$cluster_compare <- 
  ifelse(main_data_table$ref_cluster ==
           main_data_table$query_cluster,
         as.character(main_data_table$ref_cluster),
         "different")

main_data_table$subcluster_compare <- 
  ifelse(main_data_table$ref_subcluster ==
           main_data_table$query_subcluster,
         as.character(main_data_table$ref_subcluster),
         "different")

main_data_table$host_genus_compare <- 
  as.factor(main_data_table$host_genus_compare)

main_data_table$cluster_compare <- 
  as.factor(main_data_table$cluster_compare)

main_data_table$subcluster_compare <- 
  as.factor(main_data_table$subcluster_compare)


# Convert all Unspecified fields to NA missing value.
main_data_table[main_data_table == "Unspecified"] <- NA


# QC: Check how all data looks.
plot_genomic_similarity_standard(main_data_table)




###
###
###
###
###
###
###
###
###
### 4. Output data for immunity analysis.
# Currently, the pairwise comparison data does not contain self-comparisons
# (i.e. L5 compared to L5), and it does not contain duplicate reciprocal pairs
# (i.e. L5-Trixie is absent if Trixie-L5 is already present). Analyzing
# immunity data requires all duplicate reciprocal data and self-comparisons.

# First subset only the columns of interest.
partial_data_for_immunity <- subset(main_data_table,
                                    select=c("mash_reference",
                                             "mash_query",
                                             "modified_mash_distance",
                                             "pham_pham_dissimilarity"))

# Create reciprocal identifiers.
partial_data_for_immunity$ref_query <- 
  paste(partial_data_for_immunity$mash_reference,
        "_",
        partial_data_for_immunity$mash_query,
        sep="")

partial_data_for_immunity$query_ref <- 
  paste(partial_data_for_immunity$mash_query,
        "_",
        partial_data_for_immunity$mash_reference,
        sep="")

partial_data_for_immunity$ref_query <- 
  as.factor(partial_data_for_immunity$ref_query)

partial_data_for_immunity$query_ref <- 
  as.factor(partial_data_for_immunity$query_ref)


duplicate_data_for_immunity1 <- subset(partial_data_for_immunity,
                                       select=c("ref_query",
                                                "modified_mash_distance",
                                                "pham_pham_dissimilarity"))
duplicate_data_for_immunity2 <- subset(partial_data_for_immunity,
                                       select=c("query_ref",
                                                "modified_mash_distance",
                                                "pham_pham_dissimilarity"))


# Rename the columns so the two tables match.
new_names <- c("phage1_phage2",
               "modified_mash_distance",
               "pham_pham_dissimilarity")

names(duplicate_data_for_immunity1) <- new_names
names(duplicate_data_for_immunity2) <- new_names



# Create a table of self-comparisons.
self_comparison_for_immunity <- subset(phage_metadata_table,
                                       select=c("phageid"))

self_comparison_for_immunity$phage1_phage2 <- 
  paste(self_comparison_for_immunity$phageid,
        "_",
        self_comparison_for_immunity$phageid,
        sep="")

self_comparison_for_immunity$modified_mash_distance <- 0
self_comparison_for_immunity$pham_pham_dissimilarity <- 0

self_comparison_for_immunity <- subset(self_comparison_for_immunity,
                                       select=c("phage1_phage2",
                                                "modified_mash_distance",
                                                "pham_pham_dissimilarity"))

# Now merge the three tables = duplicate/reciprocal data and self-comparisons.
all_data_for_immunity <- rbind(duplicate_data_for_immunity1,
                               duplicate_data_for_immunity2,
                               self_comparison_for_immunity)

all_data_for_immunity$phage1_phage2 <- 
  as.factor(all_data_for_immunity$phage1_phage2)

write.table(all_data_for_immunity,
            paste(DIR_OUTPUT,"genomic_distance_data.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

###
###
###
###
###
###
###
###
###
### 5. Output data for predicting evolutionary mode (if needed).
# Evolutionary mode is computed using the analyze_mash_network.py script as
# published in Mavrich & Hatfull 2017.
# Either use only intra-cluster boundary data or all data. If all data is used, 
# then mode will be predicted for all phages. If only subset of data is used,
# then mode will be predicted only for phages that contain at least one data 
# point within the subsetted data. Only data within intra-cluster boundaries 
# should be used to match what was done in Mavrich & Hatfull 2017.
data_for_mode_prediction <- 
  subset(main_data_table,
         main_data_table$modified_mash_distance < 0.42 &
           main_data_table$pham_pham_dissimilarity < 0.89)

data_for_mode_prediction <- subset(data_for_mode_prediction,
                                   select=c("mash_reference",
                                            "mash_query",
                                            "modified_mash_distance",
                                            "pham_pham_dissimilarity"))

write.table(data_for_mode_prediction,
            paste(DIR_OUTPUT,"data_for_mode_prediction.csv",sep=""),
            sep=",",
            row.names = FALSE,col.names = FALSE,quote=FALSE)


###
###
###
###
###
###
###
###
###
### 6. Output data for computing MaxGCDGap (if needed).
# MaxGCDGap is computed using the analyze_mash_network.py script as published
# in Pope et al. 2017.


# All pairwise comparisons can be used for this analysis.
data_for_maxgcdgap_all <- subset(main_data_table,
                                 select=c("mash_reference",
                                          "mash_query",
                                          "modified_mash_distance",
                                          "pham_pham_dissimilarity"))

write.table(data_for_maxgcdgap_all,
            paste(DIR_OUTPUT,"data_for_maxgcdgap_all.csv",sep=""),
            sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)

###
###
###
###
###
###
###
###
###
### 7. Output data for creating phage networks (if needed).
# The filter_process_mashed_data.py script filters pairwise data and only
# retains comparisons that fall within user-selected boundaries. 
# It also determines groups of phages that have comparisons that fall within
# those filter parameters (creating phage networks).
# This R code performs the filtering step first.
data_for_mash_network_nuc005 <- 
  subset(main_data_table,
         main_data_table$modified_mash_distance < 0.05,
         select=c("mash_reference",
                  "mash_query",
                  "modified_mash_distance",
                  "mash_pvalue",
                  "mash_count",
                  "mash_ref_query"))

write.table(data_for_mash_network_nuc005,
            paste(DIR_OUTPUT,"data_for_mash_network_nuc005.csv",sep=""),
            sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)


###

