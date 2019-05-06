# R script to perform analyze stoperator data for Mavrich & Hatfull, mBio, 2019.
# The scripts computes position weight matrices for stoperators and
# compares similarities between PWMs.
# Distinct analyses and code blocks separated by "###".

### 1. Prepare environment.
### 2. Define functions.
### 3. Import datasets.
### 4. Re-process sites using Biostrings.
### 5. Create PWMs using TFBSTools.
### 6. Predict sites using PWMs.


### 1. Prepare environment.

# Install dependencies
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("TFBSTools")
library("Biostrings")
library("TFBSTools")


# Set working directory variables specific to local directory structure.
DIR_INPUT_MISC = "~/scratch/process_stoperator_data/input_misc/"
DIR_INPUT_GENOMES = "~/scratch/process_stoperator_data/input_genomes/"
DIR_OUTPUT = "~/scratch/process_stoperator_data/output/"

#Set paths for all input files.
PHAGE_METADATA_FILENAME = paste(DIR_INPUT_MISC,"phage_metadata.csv",sep="")
STOPERATOR_DATA_FILENAME = paste(DIR_INPUT_MISC,"stoperator_data.csv",sep="")

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

compute_pwm_distances <- function(pwm_list){
  
  pwm_list1 <- pwm_list
  pwm_list2 <- pwm_list
  
  pwm_output_list <- 
    vector("list",length(pwm_list1)*length(pwm_list2))
  
  count <- 1
  for (pwm1 in pwm_list1){
    for (pwm2 in pwm_list2){
      
      dist_pearson <- 1 - PWMSimilarity(pwm1,pwm2,method="Pearson")
      dist_euc <- PWMSimilarity(pwm1,pwm2,method="Euclidean")
      
      # Append the distance data to the output list.
      pwm_output_list[[count]] <- c(ID(pwm1),
                                    ID(pwm2),
                                    round(dist_pearson,2),
                                    round(dist_euc,2))
      count <- count + 1
    }
  }
  
  pwm_output_df <- as.data.frame(t(as.data.frame(pwm_output_list)))
  rownames(pwm_output_df) <- NULL
  names(pwm_output_df) <- c("phage1","phage2","dist_pearson","dist_euc")
  
  pwm_output_df$phage1_phage2 <- paste(pwm_output_df$phage1,
                                       "_",
                                       pwm_output_df$phage2,
                                       sep="")
  pwm_output_df$dist_pearson <- 
    as.numeric(as.character(pwm_output_df$dist_pearson))
  pwm_output_df$dist_euc <- 
    as.numeric(as.character(pwm_output_df$dist_euc))
  
  return(pwm_output_df)
  
}

plot_scatter <- function(table,field1,field2,x_range,y_range,title){
  
  plot(table[,field1],
       table[,field2],
       pch=16,cex=2,cex.axis=2,las=1,
       xlim = x_range,ylim = y_range,
       xlab="",ylab="",main = title,col="black")
  abline(0,1)  

  correlation <- lm(table[,field1] ~ 
                      table[,field2],
                    data = freq_table)
  print(summary(correlation))

  print(paste("Total number of sites in Field 1: ",sum(table[,field1])))
  print(paste("Total number of sites in Field 2: ",sum(table[,field2])))
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

# Import all 327 Cluster A genomes from the Actino1321 database.
# Read individual fasta files by generating a list of files in the directory,
# and passing that list to the readDNAStringSet function.
setwd(DIR_INPUT_GENOMES)
actino1321_genomes <- readDNAStringSet(list.files())
setwd(DIR_OUTPUT)

# Import manually-compiled list of MEME-identified stoperators.
# Data structure:
# 1. phage name (imported as factor)
# 2. stoperator_id (imported as factor)
# 3. sequence (imported as factor)
meme_sites <- read.csv(STOPERATOR_DATA_FILENAME,sep=",",header=TRUE)

# Tally # of stoperators identified by MEME as input.
meme_freq <- as.data.frame(table(meme_sites$phage))
names(meme_freq) <- c("phage","meme_frequency")


# Import phage metadata.
# This data is derived from the same metadata input used for the
# analyze_immunity_data.R. It contains data for all phages in the Actino1321
# database.
# Data structure:
# 1. phageid (imported as factor)
# 2. size (imported as int)
# 3. coordinate_genome_center (imported as int)
phage_metadata <- read.csv(PHAGE_METADATA_FILENAME,sep=",",header=TRUE)

#Phages not in Cluster A do not have a specific genome center coordinate.
phage_metadata[phage_metadata == "Unspecified"] <- NA
phage_metadata$coordinate_genome_center <- 
  as.numeric(as.character(phage_metadata$coordinate_genome_center))

###
###
###
###
###
###
###
###
###
### 4. Re-process sites using Biostrings.
# Confirm presence of sites in each genomes in case genome sequence has
# changed since original MEME search.

# Create an empty dataframe to store all the search data
bios_sites <- data.frame(phage="empty",
                         forward_seq="empty",
                         reverse_seq="empty",
                         start=0,
                         end=0,
                         strand="empty")

# For each phage:
# 1. Create list of unique 13bp stoperators.
# 2. Search both strands of the corresponding parent genome.
# 3. Save list of sites and coordinates.
for (phageid in levels(meme_sites$phage)){
  
  reduced_table_seqs <- subset(meme_sites,
                               meme_sites$phage == phageid,
                               select = c("sequence"))
  
  reduced_table_seqs$sequence <- factor(reduced_table_seqs$sequence)
  
  for (stoperator in levels(reduced_table_seqs$sequence)){
    
    stop_forward <- DNAString(stoperator)
    stop_reverse <- reverseComplement(stop_forward)
    genome_to_search <- actino1321_genomes[[phageid]]
    
    hits_plus <- matchPattern(stop_forward,genome_to_search)
    hits_minus <- matchPattern(stop_reverse,genome_to_search)
    
    if (length(hits_plus) > 0){
      
      hits_plus_df <- 
        data.frame(phage=phageid,
                   forward_seq = as.character(stop_forward),
                   reverse_seq = as.character(stop_reverse),
                   start=start(hits_plus),
                   end=end(hits_plus),
                   strand="forward")
      bios_sites <- rbind(bios_sites,hits_plus_df)
    }

    if (length(hits_minus) > 0){
      
      hits_minus_df <- 
        data.frame(phage=phageid,
                   forward_seq=as.character(stop_forward),
                   reverse_seq=as.character(stop_reverse),
                   start=start(hits_minus),
                   end=end(hits_minus),
                   strand="reverse")
      bios_sites <- rbind(bios_sites,hits_minus_df)
    }
  }
}

# Remove the row of empty data used to initiate the table.
bios_sites <- bios_sites[bios_sites$phage != "empty",]
bios_sites$phage <- factor(bios_sites$phage)
bios_sites$strand <- factor(bios_sites$strand)
bios_sites$forward_seq <- factor(bios_sites$forward_seq)
bios_sites$reverse_seq <- factor(bios_sites$reverse_seq)

# Create unique identifier for each site.
bios_sites$stop_site_id <- paste(bios_sites$phage,
                                 bios_sites$strand,
                                 bios_sites$start,
                                 sep="_")

bios_sites$stop_site_id <- factor(bios_sites$stop_site_id)

# Add metadata.
bios_sites <- merge(bios_sites,phage_metadata,by.x="phage",by.y="phageid")
bios_sites$dist_from_center <- 
  bios_sites$start - bios_sites$coordinate_genome_center

# Tally # of stoperators confirmed.
bios_freq <- as.data.frame(table(bios_sites$phage))

names(bios_freq) <- c("phage","biostrings_freq")

# Tally # of stoperators on each side of genome center.
bios_sites_left <- subset(bios_sites,bios_sites$dist_from_center <= 0)
bios_sites_right <- subset(bios_sites,bios_sites$dist_from_center > 0)

bios_freq_left <- as.data.frame(table(bios_sites_left$phage))
bios_freq_right <- as.data.frame(table(bios_sites_right$phage))

names(bios_freq_left) <- c("phage","biostrings_freq_left")
names(bios_freq_right) <- c("phage","biostrings_freq_right")


# QC: should equal 0.
nrow(bios_sites) - nrow(bios_sites_left) - nrow(bios_sites_right)


# Tally # of stoperators on each strand of each side of genome center.
bios_sites_left_for <- subset(bios_sites_left,
                              bios_sites_left$strand == "forward")
bios_sites_left_rev <- subset(bios_sites_left,
                              bios_sites_left$strand == "reverse")
bios_sites_right_for <- subset(bios_sites_right,
                               bios_sites_right$strand == "forward")
bios_sites_right_rev <- subset(bios_sites_right,
                               bios_sites_right$strand == "reverse")


bios_freq_left_for <- as.data.frame(table(bios_sites_left_for$phage))
bios_freq_left_rev <- as.data.frame(table(bios_sites_left_rev$phage))
bios_freq_right_for <- as.data.frame(table(bios_sites_right_for$phage))
bios_freq_right_rev <- as.data.frame(table(bios_sites_right_rev$phage))


names(bios_freq_left_for) <- c("phage","biostrings_freq_left_forward")
names(bios_freq_left_rev) <- c("phage","biostrings_freq_left_reverse")
names(bios_freq_right_for) <- c("phage","biostrings_freq_right_forward")
names(bios_freq_right_rev) <- c("phage","biostrings_freq_right_reverse")


# Combine MEME and Biostrings data.
freq_table <- merge(meme_freq,bios_freq,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left_for,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left_rev,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right_for,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right_rev,by.x="phage",by.y="phage")



freq_table$meme_biostrings_diff <- 
  freq_table$meme_frequency - freq_table$biostrings_freq

freq_table$biostrings_percent_left <- 
  freq_table$biostrings_freq_left / freq_table$biostrings_freq

freq_table$biostrings_percent_right <- 
  freq_table$biostrings_freq_right / freq_table$biostrings_freq

freq_table$biostrings_percent_left_forward <- 
  freq_table$biostrings_freq_left_forward / freq_table$biostrings_freq_left

freq_table$biostrings_percent_left_reverse <- 
  freq_table$biostrings_freq_left_reverse / freq_table$biostrings_freq_left

freq_table$biostrings_percent_right_forward <- 
  freq_table$biostrings_freq_right_forward / freq_table$biostrings_freq_right

freq_table$biostrings_percent_right_reverse <- 
  freq_table$biostrings_freq_right_reverse / freq_table$biostrings_freq_right


# QC: Tally checks should equal 0.
freq_table$biostrings_check1 <- 
  freq_table$biostrings_freq - 
  freq_table$biostrings_freq_left - 
  freq_table$biostrings_freq_right

freq_table$biostrings_check2 <- 
  freq_table$biostrings_freq_left - 
  freq_table$biostrings_freq_left_forward - 
  freq_table$biostrings_freq_left_reverse


freq_table$biostrings_check3 <- 
  freq_table$biostrings_freq_right - 
  freq_table$biostrings_freq_right_forward - 
  freq_table$biostrings_freq_right_reverse

summary(freq_table$biostrings_check1)
summary(freq_table$biostrings_check2)
summary(freq_table$biostrings_check3)



# QC: Percent checks should equal 1.
freq_table$biostrings_check4 <- 
  freq_table$biostrings_percent_left +
  freq_table$biostrings_percent_right

freq_table$biostrings_check5 <- 
  freq_table$biostrings_percent_left_forward +
  freq_table$biostrings_percent_left_reverse

freq_table$biostrings_check6 <- 
  freq_table$biostrings_percent_right_forward +
  freq_table$biostrings_percent_right_reverse


summary(freq_table$biostrings_check4)
summary(freq_table$biostrings_check5)
summary(freq_table$biostrings_check6)



# QC:Compare # MEME sites versus # Biostrings sites per genome.
plot(freq_table$meme_frequency,
     freq_table$biostrings_freq)

# QC: If this does not equal 0, then there is a difference between the number
# of stoperators present in the two sets of stoperators.
nrow(meme_sites) - nrow(bios_sites)
nrow(meme_sites)
nrow(bios_sites)
hist(freq_table$meme_biostrings_diff)
# MEME sites tally = 9266. Biostrings sites tally = 9273.
# Biostrings contains 7 more stoperators. Overall,  only a few genomes contain
# small differences in stoperators than expected. Seven phages contain 1 fewer
# stoperator from MEME than Biostrings:
# BeesKnees (A1)
# Changeling (A2)
# Doom (A1)
# KBG (A1)
# Petruchio (A1)
# Saintus (A8)
# SarFire (A1)





# Plot distribution of number of stoperators.
par(mar=c(4,8,8,4))
hist(table(bios_sites$phage),col="black",
     xlim=c(0,50),ylim=c(0,80),
     breaks=10,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,paste(DIR_OUTPUT,
                   'biostrings_stoperators_per_genome.pdf',
                   sep=""))
dev.off()

# Output biostrings data for reference
# This list is more robust than the input meme-derived list, because it has 
# been generated using the most recent version of genome sequences.
write.table(bios_sites,
            paste(DIR_OUTPUT,"biostrings_stoperators.csv",sep=""),
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
### 5. Create PWMs using TFBSTools.
# Use the sites that were identified in each genome by Biostrings search.
# For each phage in the table, create a pwm.
# Create two lists of PWMs = one to store the linear probability and one to
# store the log2 probability ratio.
pwm_prob_list <- vector("list",nlevels(bios_sites$phage))
pwm_log2_list <- vector("list",nlevels(bios_sites$phage))


# For each phage genome, the list of stoperators must contain at least 
# one of each of the 4 nucleotides A,C,T,G. Otherwise, when PFMatrix is called,
# an error is encountered. The for loop will probably crash if one or more
# PWMs does not meet this criteria, but a QC step has been used to confirm it.
alphabet_check_sum <- 0

count <- 1
for (stoperator_phage in levels(bios_sites$phage)){
  
  reduced_table_seqs <- subset(bios_sites,
                               bios_sites$phage == 
                                 stoperator_phage,
                               select = c("forward_seq"))
  
  reduced_table_seqs <- as.character(reduced_table_seqs$forward_seq)
  
  consensus_matrix <- consensusMatrix(reduced_table_seqs)
  stoperator_pfm <- PFMatrix(ID=stoperator_phage,
                                  name=stoperator_phage,
                                  profileMatrix = consensus_matrix)
  
  pwm_log2 <- toPWM(stoperator_pfm,type="log2probratio") 
  pwm_prob <- toPWM(stoperator_pfm,type="prob") 

  # Append the pwm to the list of pwms.
  pwm_log2_list[[count]] <- pwm_log2
  pwm_prob_list[[count]] <- pwm_prob  
  
  
  # Append DNA letter frequency check.
  reduced_table_seqs_dnastringset <- DNAStringSet(reduced_table_seqs)
  all_4_nucleotides <- 
    ifelse(length(uniqueLetters(reduced_table_seqs_dnastringset)) == 4,1,0)
  
  alphabet_check_sum <- 
    alphabet_check_sum + all_4_nucleotides
  
  count <- count + 1
  
}  


# QC: should equal 0.
nlevels(bios_sites$phage) - alphabet_check_sum


### Compute PWM distances using the Prob PWMs or Log2ProbRatio PWMs.
# The log2-based score seems to be more commonly used than a linear score.
pwm_log2_distances <- compute_pwm_distances(pwm_log2_list)
pwm_prob_distances <- compute_pwm_distances(pwm_prob_list)




# Compare the euclidean distance of PWM linear and the log2probratio
# probabilities.
names(pwm_prob_distances) <- paste("prob_",names(pwm_prob_distances),sep = "")

pwm_distances <- merge(pwm_log2_distances,
                       pwm_prob_distances,
                       by.x="phage1_phage2",
                       by.y="prob_phage1_phage2")


#QC: Should equal 0.
nrow(pwm_log2_distances) - nrow(pwm_prob_distances)
nrow(pwm_distances) - nrow(pwm_log2_distances)
nrow(pwm_distances) - nrow(pwm_prob_distances)

plot(pwm_distances$dist_euc,
     pwm_distances$prob_dist_euc)

plot(pwm_distances$dist_euc,
     pwm_distances$dist_pearson)

plot(pwm_distances$prob_dist_euc,
     pwm_distances$prob_dist_pearson)


#TODO
#Summary of comparisons















#TODO will need to select columns
# Export data
write.table(pwm_log2_distances,
            paste(DIR_OUTPUT,"stoperator_pwm_data.csv",sep=""),
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
### 6. Predict sites using PWMs.
### After creating all stoperator PWMs, iterate through all genomes for each
# PWM and predict sites. 
time_start <- Sys.time()

# A PWMatrixList needs to be created first.
stoperator_pwm_matrixlist <- do.call(PWMatrixList,pwm_log2_list)

# Use default min.score of 80%, so that different cutoffs can be assessed.
stop_site_set80_list <- searchSeq(stoperator_pwm_matrixlist,actino1321_genomes)
stop_site_set80_df <- as(stop_site_set80_list,"data.frame")
time_stop <- Sys.time()


# QC: keep track of how long this step takes.
time_start
time_stop
# Note: it takes about an hour to run the two commands above using stop.

#TODO remove
# stop_site_set80_list_backup <- stop_site_set80_list
# stop_site_set80_df_backup <- stop_site_set80_df
#


# Dataframe output: each row is a predicted site in a genome from one PWM
# Column names:
# 1. seqnames = factor, 327 levels are individual phage genomes
# 2. source = factor, 1 level ("TFBS")
# 3. feature = factor, 1 level ("TFBS")
# 4. start = int
# 5. end = int
# 6. absScore = number (from ~4 to ~25)
# 7. relScore = number (from 0.8 to 1). It reflects the range of possible
#               scores that could be obtained with the motif and is
#               normalized between 0 and 1.
# 8. strand = factor, 2 levels ("+","-")
# 9. ID = factor, 327 levels are stoperator PWM names. Identical to "TF"
# 10. TF = factor, 327 levels are stoperator PWM names. Identical to "ID"
# 11. class = factor, 1 level ("Unknown")
# 12. siteSeqs = chr, 13bp stoperator DNA sequence


# Reduce data for analysis with immunity data
site_predictions_80 <- subset(stop_site_set80_df,
                              select=c("seqnames",
                                       "start",
                                       "end",
                                       "strand",
                                       "siteSeqs",
                                       "absScore",
                                       "relScore",
                                       "ID"))

names(site_predictions_80) <- c("stoperator_target",
                                "site_start",
                                "site_end",
                                "site_strand",
                                "site_seq",
                                "site_abs_score",
                                "site_rel_score",
                                "stoperator_motif")


site_predictions_80$motif_target <- 
  paste(site_predictions_80$stoperator_motif,
        "_",
        site_predictions_80$stoperator_target,
        sep="")

site_predictions_80$motif_target <- 
  as.factor(site_predictions_80$motif_target)


#TODO
# QC: Check number of unique comparisons.
nlevels(site_predictions_80$motif_target)
# This list has 106,300 levels, but in the Actino1321 database there are
# 106,929 unique possible pairwise combinations (327 PWMs * 327 target genomes).
# Re-factor motif_target using all possible pairwise combinations.
# This will enable quantification for all combinations, 
# even those that are missing.
all_levels <- 
  expand.grid(stoperator_target = 
                as.character(levels(site_predictions_80$stoperator_target)),
              stoperator_motif = 
                as.character(levels(site_predictions_80$stoperator_motif)))

all_levels$motif_target <- paste(all_levels$stoperator_motif,
                                 "_",
                                 all_levels$stoperator_target,
                                 sep="")

all_levels$motif_target <- as.factor(all_levels$motif_target)


site_predictions_80$motif_target <- 
  factor(site_predictions_80$motif_target,
         levels(all_levels$motif_target))

# Should equal 0.
(327*327) - nlevels(site_predictions_80$motif_target)
# Now there are 106,929 levels.




# Create unique site ids
site_predictions_80$site_strand2 <- 
  ifelse(site_predictions_80$site_strand == "+",
         "forward",
         "reverse")
site_predictions_80$site_strand2 <- 
  factor(site_predictions_80$site_strand2)

site_predictions_80$stop_site_id <- 
  paste(site_predictions_80$stoperator_target,
        site_predictions_80$site_strand2,
        site_predictions_80$site_start,
        sep="_")

site_predictions_80$stop_site_id <- 
  factor(site_predictions_80$stop_site_id)




site_predictions_80_self <- 
  subset(site_predictions_80,
         as.character(site_predictions_80$stoperator_target) == 
           as.character(site_predictions_80$stoperator_motif))


site_predictions_80_self$stoperator_target <- 
  factor(site_predictions_80_self$stoperator_target)
site_predictions_80_self$stoperator_motif <- 
  factor(site_predictions_80_self$stoperator_motif)
site_predictions_80_self$motif_target <- 
  factor(site_predictions_80_self$motif_target)
site_predictions_80_self$site_strand2 <- 
  factor(site_predictions_80_self$site_strand2)



site_predictions_85_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.85)
site_predictions_86_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.86)
site_predictions_87_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.87)
site_predictions_88_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.88)
site_predictions_89_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.89)
site_predictions_90_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.90)
site_predictions_91_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.91)
site_predictions_92_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.92)
site_predictions_93_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.93)
site_predictions_94_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.94)
site_predictions_95_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 0.95)
site_predictions_100_self <- 
  subset(site_predictions_80_self,
         site_predictions_80_self$site_rel_score >= 1)


# Table creates a two column table = the levels and the frequency (int)

site_predictions_80_self_freq <- 
  as.data.frame(table(site_predictions_80_self$stoperator_target))
names(site_predictions_80_self_freq) <- c("phage","tfbstools80_freq")

site_predictions_85_self_freq <- 
  as.data.frame(table(site_predictions_85_self$stoperator_target))
names(site_predictions_85_self_freq) <- c("phage","tfbstools85_freq")

site_predictions_86_self_freq <- 
  as.data.frame(table(site_predictions_86_self$stoperator_target))
names(site_predictions_86_self_freq) <- c("phage","tfbstools86_freq")

site_predictions_87_self_freq <- 
  as.data.frame(table(site_predictions_87_self$stoperator_target))
names(site_predictions_87_self_freq) <- c("phage","tfbstools87_freq")

site_predictions_88_self_freq <- 
  as.data.frame(table(site_predictions_88_self$stoperator_target))
names(site_predictions_88_self_freq) <- c("phage","tfbstools88_freq")

site_predictions_89_self_freq <- 
  as.data.frame(table(site_predictions_89_self$stoperator_target))
names(site_predictions_89_self_freq) <- c("phage","tfbstools89_freq")

site_predictions_90_self_freq <- 
  as.data.frame(table(site_predictions_90_self$stoperator_target))
names(site_predictions_90_self_freq) <- c("phage","tfbstools90_freq")

site_predictions_91_self_freq <- 
  as.data.frame(table(site_predictions_91_self$stoperator_target))
names(site_predictions_91_self_freq) <- c("phage","tfbstools91_freq")

site_predictions_92_self_freq <- 
  as.data.frame(table(site_predictions_92_self$stoperator_target))
names(site_predictions_92_self_freq) <- c("phage","tfbstools92_freq")

site_predictions_93_self_freq <- 
  as.data.frame(table(site_predictions_93_self$stoperator_target))
names(site_predictions_93_self_freq) <- c("phage","tfbstools93_freq")

site_predictions_94_self_freq <- 
  as.data.frame(table(site_predictions_94_self$stoperator_target))
names(site_predictions_94_self_freq) <- c("phage","tfbstools94_freq")

site_predictions_95_self_freq <- 
  as.data.frame(table(site_predictions_95_self$stoperator_target))
names(site_predictions_95_self_freq) <- c("phage","tfbstools95_freq")

site_predictions_100_self_freq <- 
  as.data.frame(table(site_predictions_100_self$stoperator_target))
names(site_predictions_100_self_freq) <- c("phage","tfbstools100_freq")


# Combine TFBSTools data with MEME and Biostrings data.
freq_table <- merge(freq_table,
                    site_predictions_80_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_85_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_86_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_87_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_88_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_89_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_90_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_91_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_92_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_93_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_94_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_95_self_freq,
                    by.x="phage",
                    by.y="phage")

freq_table <- merge(freq_table,
                    site_predictions_100_self_freq,
                    by.x="phage",
                    by.y="phage")





# QC: Check how well TFBSTools matches with Biostrings and determine which
# % cutoff is best to use.
x_coords <- c(0,100)
y_coords <- c(0,100)


plot_scatter(freq_table,"biostrings_freq","tfbstools80_freq",
             x_coords,y_coords,"Cutoff: 80")
plot_scatter(freq_table,"biostrings_freq","tfbstools85_freq",
             x_coords,y_coords,"Cutoff: 85")
plot_scatter(freq_table,"biostrings_freq","tfbstools86_freq",
             x_coords,y_coords,"Cutoff: 86")
plot_scatter(freq_table,"biostrings_freq","tfbstools87_freq",
             x_coords,y_coords,"Cutoff: 87")
plot_scatter(freq_table,"biostrings_freq","tfbstools88_freq",
             x_coords,y_coords,"Cutoff: 88")
plot_scatter(freq_table,"biostrings_freq","tfbstools89_freq",
             x_coords,y_coords,"Cutoff: 89")
plot_scatter(freq_table,"biostrings_freq","tfbstools90_freq",
             x_coords,y_coords,"Cutoff: 90")
plot_scatter(freq_table,"biostrings_freq","tfbstools91_freq",
             x_coords,y_coords,"Cutoff: 91")
plot_scatter(freq_table,"biostrings_freq","tfbstools92_freq",
             x_coords,y_coords,"Cutoff: 92")
plot_scatter(freq_table,"biostrings_freq","tfbstools93_freq",
             x_coords,y_coords,"Cutoff: 93")
plot_scatter(freq_table,"biostrings_freq","tfbstools94_freq",
             x_coords,y_coords,"Cutoff: 94")
plot_scatter(freq_table,"biostrings_freq","tfbstools95_freq",
             x_coords,y_coords,"Cutoff: 95")
plot_scatter(freq_table,"biostrings_freq","tfbstools100_freq",
             x_coords,y_coords,"Cutoff: 100")

# QC:
# CutofF____|__Multiple_R2__|__#_sites
# biostrings|               | 9273
# 80        |  0.1943       | 13193
# 85        |  0.8067       | 9686
# 86        |  0.876        | 9536
# 87        |  0.9224       | 9399
# 88        |  0.9532       | 9202
# 89        |  0.8966       | 8965
# 90        |  0.7866       | 8720
# 91        |  0.7791       | 8535
# 92        |  0.7465       | 8373
# 93        |  0.6517       | 8124
# 94        |  0.5796       | 7869
# 95        |  0.5526       | 7493
# 100       |  3.081e-05    | 2361

# TFBSTools cutoff of 88% identifies nearly the same number of
# stoperators in each genome (9202 compared to 9273), and it produces the most
# correlated data, with an R2 of 0.9532, So using 88% is probably the best 
# cutoff to use. 


names(site_predictions_88_self) <- 
  paste("tfbs88_",names(site_predictions_88_self),sep="")



# Compute the exact position of sites present in biostrings and 
# tfbstools88 datasets.
biostrings_tfbs88 <- merge(bios_sites,
                           site_predictions_88_self,
                           by.x="stop_site_id",
                           by.y="tfbs88_stop_site_id",
                           all.x=TRUE,
                           all.y=TRUE)


# QC: evaluate which sites are present in which dataset.
biostrings_tfbs88$source_biostrings <- 
  ifelse(is.na(biostrings_tfbs88$phage),
         FALSE,TRUE)

biostrings_tfbs88$source_tfbs88 <- 
  ifelse(is.na(biostrings_tfbs88$tfbs88_stoperator_target),
         FALSE,TRUE)

biostrings_tfbs88$source_both <- 
  ifelse(biostrings_tfbs88$source_biostrings & 
           biostrings_tfbs88$source_tfbs88,
         TRUE,FALSE)

summary(biostrings_tfbs88$source_biostrings)
summary(biostrings_tfbs88$source_tfbs88)
summary(biostrings_tfbs88$source_both)
# There are 9123 sites in both Biostrings and TFBS88 (98%).
# There are 150 sites in Biostrings that are not present in TFBS88. 
# There are 79 sites in TFBS88 that are not present in Biostrings.






#TODO
# The 88% cutoff dataset now represents endogenous stoperators in each genome
# that were identified by MEME (and re-processed and confirmed with biostrings)
# as well as matching stoperators in all other genomes that match the PWM
# of the endogenous stoperator set.




# Now that the best TFBS cutoff has been determined (88%), search all genomes
# with each PWM and determine how many sites of each motif are in each genome.


# Re-factor columns in the self-site dataset.
site_predictions_88_self$tfbs88_stoperator_target <- 
  factor(site_predictions_88_self$tfbs88_stoperator_target)

site_predictions_88_self$tfbs88_stoperator_motif <- 
  factor(site_predictions_88_self$tfbs88_stoperator_motif)

site_predictions_88_self$tfbs88_motif_target <- 
  factor(site_predictions_88_self$tfbs88_motif_target)

site_predictions_88_self$tfbs88_site_strand2 <- 
  factor(site_predictions_88_self$tfbs88_site_strand2)



# Subset all sites predicted in each genome from each PWM above 88% threshold.
site_predictions_80$site_self <- 
  ifelse(is.element(site_predictions_80$stop_site_id,
                    site_predictions_88_self$tfbs88_stop_site_id),
         TRUE,FALSE)

site_predictions_88 <- subset(site_predictions_80,
                              site_predictions_80$site_rel_score >= 0.88)



site_predictions_88$stoperator_target <- 
  factor(site_predictions_88$stoperator_target)

site_predictions_88$stoperator_motif <- 
  factor(site_predictions_88$stoperator_motif)

site_predictions_88$motif_target <- 
  factor(site_predictions_88$motif_target)

site_predictions_88$site_strand2 <- 
  factor(site_predictions_88$site_strand2)

names(site_predictions_88) <- 
  paste("tfbs88_",names(site_predictions_88),sep="")

site_predictions_88_reduced <- 
  subset(site_predictions_88,select=c("tfbs88_stop_site_id",
                                      "tfbs88_motif_target",
                                      "tfbs88_stoperator_target",
                                      "tfbs88_stoperator_motif",
                                      "tfbs88_site_start",
                                      "tfbs88_site_end",
                                      "tfbs88_site_strand2",
                                      "tfbs88_site_seq",
                                      "tfbs88_site_abs_score",
                                      "tfbs88_site_rel_score",
                                      "tfbs88_site_self"))


site_predictions_88_reduced$tfbs88_site_abs_score <- 
  round(site_predictions_88_reduced$tfbs88_site_abs_score,2)

site_predictions_88_reduced$tfbs88_site_rel_score <- 
  round(site_predictions_88_reduced$tfbs88_site_rel_score,2)



# Export data for analysis with immunity data.
write.table(site_predictions_88_reduced,
            paste(DIR_OUTPUT,"stoperator_site_predictions.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)




###
