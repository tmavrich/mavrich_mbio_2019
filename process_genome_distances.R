#Script to analyze Actinobacteriophage_1321 mash and gene content data
#20180926



#The standard genomic similarity plot parameters
plot_genomic_similarity_standard <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,table$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
  abline(0,2,lty=2,lwd=3,col="grey")
  
}



setwd("~/Desktop/Project_stoperator/Actinobacteriophage_1321/8_Actino1321_v2_data/5_r_analysis/")





###Import primary mash and gene content dissimilarity data

#Import mash dataset
#Format:
#0 = reference genome
#1 = query genome
#2 = mash distance 
#3 = mash p-value
#4 = mash kmer count
#5 = reference_query
mash_table <- read.csv("./input/processed_mash_output.csv",sep=",",header=TRUE)
names(mash_table) <- c("mash_reference",
                       "mash_query",
                       "mash_distance",
                       "mash_pvalue",
                       "mash_count",
                       "mash_ref_query")



#Import phage metadata
#Format:
#0 = phagename
#1 = host_genus
#2 = cluster
#3 = subcluster
phage_metadata_table <- read.csv("./input/phage_metadata.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
phage_metadata_table[phage_metadata_table == "Unspecified"] <- NA



#Modify column names of host data and merge with mash table
#Data for all phages in processed mash table should be present in phage metadata table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Match metadata for both the query and reference phages in each pairwise comparison
names(phage_metadata_table) <- c("phageid","query_host_genus","query_cluster","query_subcluster")
main_data_table <- merge(mash_table,phage_metadata_table,by.x="mash_query",by.y="phageid")

names(phage_metadata_table) <- c("phageid","ref_host_genus","ref_cluster","ref_subcluster")

main_data_table <- merge(main_data_table,phage_metadata_table,by.x="mash_reference",by.y="phageid")

names(phage_metadata_table) <- c("phageid","host_genus","cluster","subcluster")

 



#Import pham data and merge with mash table
#Format
#0 = phage1
#1 = phage1_number_unshared_phams
#2 = phage1_shared_proportion
#3 = phage2
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
pham_table <- read.csv("./input/pairwise_pham_proportions.csv",sep=",",header=TRUE)

#Compute gene content dissimilarity
pham_table$pham_dissimilarity <- 1 - pham_table$average_shared_proportion
pham_table$jaccard_dissimilarity <- 1 - pham_table$jaccard_similarity

names(pham_table) <- c("pham_phage1","pham_phage1_number_unshared_phams","pham_phage1_shared_proportion",
                       "pham_phage2","pham_phage2_number_unshared_phams","pham_phage2_shared_proportion",
                       "pham_number_shared_phams","pham_average_shared_proportion",
                       "pham_jaccard_similarity","pham_shared_pham_distribution_mean",
                       "pham_shared_pham_distribution_median","pham_shared_pham_distribution_max",
                       "pham_unshared_pham_distribution_mean","pham_unshared_pham_distribution_median",
                       "pham_unshared_pham_distribution_max","pham_unshared_orpham_count",
                       "pham_pham_dissimilarity","pham_jaccard_dissimilarity")


#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
pham_table$pham_phage1_phage2 <- paste(pham_table$pham_phage1,"_",pham_table$pham_phage2,sep="")
pham_table$pham_phage1_phage2 <- as.factor(pham_table$pham_phage1_phage2)

#To retain all rows, be sure to keep all.x=TRUE.
#But all rows don't need to be retained = making scatter plots or histograms can cause errors if some rows are missing data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
main_data_table <- merge(main_data_table,pham_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")




#Assign filter status and change mash distance if data is not significant
#main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10 & main_data_table$size_diff_max_percent < 1,TRUE,FALSE)

#Alternatively, the max percent parameter can be omitted with minimal effects on final analysis.
main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10,TRUE,FALSE)

#Determine the max filtered distance
#filtered_table <- subset(main_data_table,main_data_table$filter == TRUE)
#summary(filtered_table$mash_distance)
#max mash distance = 0.4991900

#At this point, the max mash distance of all filtered comparisons < 0.5. So set the distance of all comparisons that did not pass the filter = 0.5
main_data_table$modified_mash_distance <- ifelse(main_data_table$filter == TRUE,main_data_table$mash_distance,0.5)





#Compare ref and query host data columns and convert columns to factor class
main_data_table$host_genus_compare <- ifelse(main_data_table$ref_host_genus==main_data_table$query_host_genus,as.character(main_data_table$ref_host_genus),"different")
main_data_table$cluster_compare <- ifelse(main_data_table$ref_cluster==main_data_table$query_cluster,as.character(main_data_table$ref_cluster),"different")
main_data_table$subcluster_compare <- ifelse(main_data_table$ref_subcluster==main_data_table$query_subcluster,as.character(main_data_table$ref_subcluster),"different")

#Now set all new columns to factor class
main_data_table$host_genus_compare <- as.factor(main_data_table$host_genus_compare)
main_data_table$cluster_compare <- as.factor(main_data_table$cluster_compare)
main_data_table$subcluster_compare <- as.factor(main_data_table$subcluster_compare)


#Convert all Unspecified fields to NA missing value
main_data_table[main_data_table == "Unspecified"] <- NA




#Check how all data looks 
plot_genomic_similarity_standard(main_data_table)






#Output data for predicting evolutionary mode with analyze_mash_network.py script
#Either use only intra-cluster boundary data or all data
#If all data is used, then mode will be predicted for all phages
#If only subset of data is used, then mode will be predicted only for phages that contain at least one data point within the subsetted data
#However, only data within intra-cluster boundaries should be used to match what was done in Mavrich & Hatfull 2017.

#data_for_mode_prediction <- main_data_table

data_for_mode_prediction <- subset(main_data_table,
                                   main_data_table$modified_mash_distance < 0.42 &
                                     main_data_table$pham_pham_dissimilarity < 0.89)


data_for_mode_prediction <- subset(data_for_mode_prediction,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))
write.table(data_for_mode_prediction,"./output/data_for_mode_prediction.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)



#Output data for computing MaxGCDGap using all pairwise comparisons
data_for_maxgcdgap_all <- subset(main_data_table,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))
write.table(data_for_maxgcdgap_all,"./output/data_for_maxgcdgap_all.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)


#Output data for computing MaxGCDGap using only Cluster A pairwise comparisons
data_for_maxgcdgap_clustera <- subset(main_data_table,
                                      main_data_table$cluster_compare == 'A',
                                      select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))

write.table(data_for_maxgcdgap_clustera,"./output/data_for_maxgcdgap_clustera.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)





#Output data for creating phage networks with filter_process_mashed_data script
#The filter script filters the data and only retains comparisons that fall within user-selected boundaries
#It also determines groups of phages that have comparisons that fall within those filter parameters (creating phage networks)
#This R code performs the filtering step first.

data_for_mash_network_nuc005 <- subset(main_data_table,
                                       main_data_table$modified_mash_distance < 0.05,
                                       select=c("mash_reference","mash_query","modified_mash_distance","mash_pvalue","mash_count","mash_ref_query"))


write.table(data_for_mash_network_nuc005,"./output/data_for_mash_network_nuc005.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)







#Currently, the pairwise data does not contain self-comparisons (i.e. L5 compared to L5), 
#and it does not contain duplicate reciprocal pairs (i.e. L5-Trixie, if Trixie-L5 is already present)
#Analyzing immunity data requires all duplicate reciprocal data and self-comparisons
#First remove only the columns of interest
partial_data_for_immunity <- subset(main_data_table,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))

#Create reciprocal identifiers
partial_data_for_immunity$ref_query <- paste(partial_data_for_immunity$mash_reference,"_",partial_data_for_immunity$mash_query,sep="")
partial_data_for_immunity$query_ref <- paste(partial_data_for_immunity$mash_query,"_",partial_data_for_immunity$mash_reference,sep="")
partial_data_for_immunity$ref_query <- as.factor(partial_data_for_immunity$ref_query)
partial_data_for_immunity$query_ref <- as.factor(partial_data_for_immunity$query_ref)

duplicate_data_for_immunity1 <- subset(partial_data_for_immunity,select=c("ref_query","modified_mash_distance","pham_pham_dissimilarity"))
duplicate_data_for_immunity2 <- subset(partial_data_for_immunity,select=c("query_ref","modified_mash_distance","pham_pham_dissimilarity"))

#Rename the columns so the two tables match
names(duplicate_data_for_immunity1) <- c("phage1_phage2","modified_mash_distance","pham_pham_dissimilarity")
names(duplicate_data_for_immunity2) <- c("phage1_phage2","modified_mash_distance","pham_pham_dissimilarity")



#Create a table of self-comparisons
self_comparison_for_immunity <- subset(phage_metadata_table,select=c("phageid"))
self_comparison_for_immunity$phage1_phage2 <- paste(self_comparison_for_immunity$phageid,"_",self_comparison_for_immunity$phageid,sep="")
self_comparison_for_immunity$modified_mash_distance <- 0
self_comparison_for_immunity$pham_pham_dissimilarity <- 0
self_comparison_for_immunity <- subset(self_comparison_for_immunity,select=c("phage1_phage2","modified_mash_distance","pham_pham_dissimilarity"))

#Now merge the three tables = duplicate/reciprocal data and self-comparisons
all_data_for_immunity <- rbind(duplicate_data_for_immunity1,duplicate_data_for_immunity2,self_comparison_for_immunity)
all_data_for_immunity$phage1_phage2 <- as.factor(all_data_for_immunity$phage1_phage2)


write.table(all_data_for_immunity,"./output/all_data_for_immunity.csv",sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)






